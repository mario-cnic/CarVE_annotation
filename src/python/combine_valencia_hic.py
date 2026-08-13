#!/usr/bin/env python3
"""
combine_valencia_hic.py

Standardize cDNA names from raw Oracle (cases) and NDG (controls) datasets,
query GeneBe in batches to obtain GRCh38 coordinates, rename existing hg19 coordinate
columns to old_* prefixes, perform a genomic outer join on GRCh38 coords, and save
to integrated Parquet format.

Usage:
    mamba run -p /home/mruizp/data_lab_PGP/shared/utils/conda_envs/liftover python3 src/python/combine_valencia_hic.py
"""
import os
import sys
import glob
import time
import re
import json
import pandas as pd
import numpy as np
# pyrefly: ignore [missing-import]
import genebe as gnb
# pyrefly: ignore [missing-import]
from pyliftover import LiftOver

# Batch size for GeneBe queries
BATCH_SIZE = 100

# Initialize pyliftover LiftOver lazy-load variables
lo = None
lo_initialized = False

def save_with_retry(df, path, max_retries=5, delay=2):
    """
    Saves a DataFrame to Parquet with retries and exponential backoff to handle
    transient network filesystem (NFS) hiccups.
    Uses an atomic temp-file write and rename to prevent file corruption in case
    of interruptions (e.g., job cancellation or SIGTERM).
    """
    tmp_path = f"{path}.tmp"
    for attempt in range(max_retries):
        try:
            # Write to a temporary file first
            df.to_parquet(tmp_path, index=False)
            # Perform atomic replace
            os.replace(tmp_path, path)
            return
        except OSError as e:
            # Clean up temp file if it was created
            if os.path.exists(tmp_path):
                try:
                    os.remove(tmp_path)
                except Exception:
                    pass
            if attempt == max_retries - 1:
                raise e
            print(f"Warning: OSError writing to {path} (attempt {attempt+1}/{max_retries}): {e}. Retrying in {delay}s...")
            time.sleep(delay)
            delay *= 2

def check_manifest_match(gene, oracle_path, ndg_path, out_dir):
    """
    Checks if the size and modification time (mtime) of the input raw files
    match the records stored in the output directory's hidden manifest.
    Returns True if they match, False otherwise.
    """
    manifest_path = os.path.join(out_dir, ".combined_manifest.json")
    if not os.path.exists(manifest_path):
        return False
    try:
        with open(manifest_path, 'r') as f:
            manifest = json.load(f)
        entry = manifest.get(gene)
        if not entry:
            return False
            
        current_o_size = os.path.getsize(oracle_path)
        current_o_mtime = os.path.getmtime(oracle_path)
        current_n_size = os.path.getsize(ndg_path)
        current_n_mtime = os.path.getmtime(ndg_path)
        
        return (
            entry.get('oracle_size') == current_o_size and
            entry.get('oracle_mtime') == current_o_mtime and
            entry.get('ndg_size') == current_n_size and
            entry.get('ndg_mtime') == current_n_mtime
        )
    except Exception:
        return False

def update_manifest(gene, oracle_path, ndg_path, out_dir):
    """
    Updates the destination directory's hidden manifest with the current sizes
    and modification times (mtimes) of the input raw files.
    """
    manifest_path = os.path.join(out_dir, ".combined_manifest.json")
    os.makedirs(out_dir, exist_ok=True)
    manifest = {}
    if os.path.exists(manifest_path):
        try:
            with open(manifest_path, 'r') as f:
                manifest = json.load(f)
        except Exception:
            manifest = {}
            
    try:
        manifest[gene] = {
            "oracle_size": os.path.getsize(oracle_path),
            "oracle_mtime": os.path.getmtime(oracle_path),
            "ndg_size": os.path.getsize(ndg_path),
            "ndg_mtime": os.path.getmtime(ndg_path)
        }
        with open(manifest_path, 'w') as f:
            json.dump(manifest, f, indent=2)
    except Exception as e:
        print(f"Warning: Failed to update manifest at {manifest_path}: {e}")

def standardize_chrom(chrom_num):
    """
    Standardize chromosome identifiers to chrN format (e.g. 3 -> chr3, MT -> chrm, CHROMOSOME_3 -> chr3).
    """
    if pd.isna(chrom_num):
        return None
    val = str(chrom_num).strip().upper()
    # Strip prefix CHROMOSOME_ or CHROMOSOME or CHR
    if val.startswith("CHROMOSOME_"):
        val = val[11:]
    elif val.startswith("CHROMOSOME"):
        val = val[10:]
    elif val.startswith("CHR"):
        val = val[3:]
        
    if val in ("MT", "M", "CHRM"):
        return "chrm"
    return f"chr{val.lower()}"

def query_genebe_batch(cdna_list):
    """
    Queries genebe.parse_variants in batches and returns a dictionary mapping
    cDNA strings to (CHROM, POS, REF, ALT) tuples in GRCh38.
    """
    mapping = {}
    total = len(cdna_list)
    print(f"Querying GeneBe for {total} unique cDNA accessions in batches of {BATCH_SIZE}...")

    for i in range(0, total, BATCH_SIZE):
        batch = cdna_list[i : i + BATCH_SIZE]
        results = None
        max_retries = 3
        delay = 2
        for attempt in range(max_retries):
            try:
                results = gnb.parse_variants(batch)
                break
            except Exception as e:
                print(f"Warning: GeneBe query attempt {attempt+1}/{max_retries} failed: {e}")
                if attempt == max_retries - 1:
                    print(f"Error: Batch starting at index {i} failed after {max_retries} attempts.")
                else:
                    time.sleep(delay)
                    delay *= 2

        if results:
            for cdna, res in zip(batch, results):
                if not res:
                    mapping[cdna] = (None, None, None, None)
                    continue
                # Expected output: e.g. '3-183651103-C-T'
                fields = res.split("-")
                if len(fields) == 4:
                    chrom = standardize_chrom(fields[0])
                    try:
                        pos = int(fields[1])
                    except ValueError:
                        pos = None
                    ref = fields[2].strip().upper()
                    alt = fields[3].strip().upper()
                    mapping[cdna] = (chrom, pos, ref, alt)
                else:
                    mapping[cdna] = (None, None, None, None)
        else:
            for cdna in batch:
                mapping[cdna] = (None, None, None, None)

    return mapping

def liftover_coordinate_hg19_to_hg38(chrom, pos):
    """
    Lifts over a coordinate from hg19 to hg38 using pyliftover.
    Checks that the chromosome and position represent standard hg19 coordinates first.
    """
    global lo, lo_initialized
    if not lo_initialized:
        lo_initialized = True
        try:
            local_chain = "_RAW/data_Valencia_plus_hic/hg19ToHg38.over.chain.gz"
            if os.path.exists(local_chain):
                print(f"Initializing LiftOver with local chain: {local_chain}")
                lo = LiftOver(local_chain)
            else:
                print("Initializing LiftOver from UCSC (internet connection required if not cached)...")
                import socket
                socket.setdefaulttimeout(10)
                lo = LiftOver('hg19', 'hg38')
        except Exception as e:
            print(f"Warning: Could not initialize LiftOver fallback: {e}")
            print("If you are offline, please download the chain file from:")
            print("  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz")
            print("and place it at: _RAW/data_Valencia_plus_hic/hg19ToHg38.over.chain.gz")
            print("Coordinates-based liftover fallback will be skipped for failed cDNA mappings.")
            lo = None

    if not lo:
        return None, None
    if pd.isna(chrom) or pd.isna(pos):
        return None, None

    # Standardize chrom format to standard UCSC camel-case e.g. chr3, chrX
    chrom_str = standardize_chrom(chrom)
    if not chrom_str:
        return None, None
    
    # Map to proper case expected by pyliftover
    if chrom_str == "chrx":
        chrom_str = "chrX"
    elif chrom_str == "chry":
        chrom_str = "chrY"
    elif chrom_str == "chrm":
        chrom_str = "chrM"

    # Verify standard hg19 assembly chromosome list
    valid_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    if chrom_str not in valid_chroms:
        return None, None

    try:
        pos_int = int(float(pos))
        if pos_int <= 0:
            return None, None
    except ValueError:
        return None, None

    try:
        # UCSC coordinates are 0-based in pyliftover.
        # Convert input 1-based coordinate to 0-based, then convert output back to 1-based.
        res = lo.convert_coordinate(chrom_str, pos_int - 1)
        if res:
            target_chrom, target_pos, target_strand, target_size = res[0]
            # Convert target_chrom back to standard lowercase to match the rest of the script
            target_chrom_std = standardize_chrom(target_chrom)
            return target_chrom_std, int(target_pos) + 1
    except Exception as e:
        print(f"Error in liftover for {chrom_str}:{pos_int}: {e}")

    return None, None

def normalize_indel_hg38(chrom, pos, ref, alt, abs_desc, lo, fa):
    """
    Normalizes and left-aligns an indel using liftover and the hg38 reference genome FASTA.
    Returns (target_chrom, target_pos, target_ref, target_alt) in GRCh38.
    """
    if not lo or not fa:
        return None, None, None, None

    # Parse raw alleles
    raw_ref = str(ref).strip() if pd.notna(ref) else ""
    raw_alt = str(alt).strip() if pd.notna(alt) else ""
    
    # Clean up empty strings or hyphens
    if raw_ref == "-": raw_ref = ""
    if raw_alt == "-": raw_alt = ""

    # Parse hg19 coordinate and length from abs_desc (NC_...) if present
    start_pos = None
    deleted_len = 0
    ins_seq = ""
    is_del = False
    is_ins = False

    if pd.notna(abs_desc):
        abs_desc_str = str(abs_desc).strip()
        # 1. Check deletion range
        m_del_range = re.search(r'g\.(\d+)_(\d+)del([ACGTN]*)', abs_desc_str)
        if m_del_range:
            start_pos = int(m_del_range.group(1))
            end_pos = int(m_del_range.group(2))
            deleted_len = end_pos - start_pos + 1
            is_del = True
            m_ins = re.search(r'ins([ACGTN]+)', abs_desc_str)
            if m_ins:
                ins_seq = m_ins.group(1)
        else:
            # 2. Check single deletion
            m_del_single = re.search(r'g\.(\d+)del([ACGTN]*)', abs_desc_str)
            if m_del_single:
                start_pos = int(m_del_single.group(1))
                deleted_len = 1
                is_del = True
            else:
                # 3. Check insertion range
                m_ins_range = re.search(r'g\.(\d+)_(\d+)ins([ACGTN]+)', abs_desc_str)
                if m_ins_range:
                    start_pos = int(m_ins_range.group(1)) # preceding base position
                    ins_seq = m_ins_range.group(3)
                    is_ins = True

    # Fallback to parsing from raw_ref / raw_alt if abs_desc parsing failed
    if not is_del and not is_ins:
        if raw_alt == "" and raw_ref != "":
            is_del = True
            start_pos = int(pos)
            deleted_len = len(raw_ref)
        elif raw_ref == "" and raw_alt != "":
            is_ins = True
            start_pos = int(pos)
            ins_seq = raw_alt
        else:
            return None, None, None, None

    # Determine preceding base position in hg19 (1-based coordinate)
    if is_del:
        hg19_vcf_pos = start_pos - 1
    else: # is_ins
        hg19_vcf_pos = start_pos

    # Liftover the preceding base position
    target_chrom, target_pos = liftover_coordinate_hg19_to_hg38(chrom, hg19_vcf_pos)
    if not target_chrom or not target_pos:
        return None, None, None, None

    # Map to proper case expected by pysam FASTA (e.g. chrx -> chrX)
    fasta_chrom = target_chrom
    if fasta_chrom == "chrx":
        fasta_chrom = "chrX"
    elif fasta_chrom == "chry":
        fasta_chrom = "chrY"
    elif fasta_chrom == "chrm":
        fasta_chrom = "chrM"

    # Fetch reference base at target_pos in hg38 (0-based start)
    try:
        ref_base = fa.fetch(fasta_chrom, target_pos - 1, target_pos).upper()
    except Exception as e:
        print(f"Warning: Fasta fetch failed at {fasta_chrom}:{target_pos}: {e}")
        return None, None, None, None

    if is_del:
        # Fetch the entire reference sequence (preceding base + deleted length) in hg38
        try:
            ref_seq = fa.fetch(fasta_chrom, target_pos - 1, target_pos - 1 + deleted_len + 1).upper()
        except Exception:
            ref_seq = ref_base
        alt_seq = ref_base + ins_seq
    else: # is_ins
        ref_seq = ref_base
        alt_seq = ref_base + ins_seq

    return target_chrom, target_pos, ref_seq, alt_seq

def deduplicate_dataset(df, is_oracle=True):
    """
    Groups the dataset by ['CHROM', 'POS', 'REF', 'ALT'] and de-duplicates them on GRCh38.
    Aggregates:
    - Text and ID columns: unique semicolon-separated values
    - Count fields: max of the values (documenting this choice to avoid double-counting same variants)
    - Other columns: first non-null value
    """
    # Ensure all ID, coordinate, and text columns are consistently string typed to prevent pyarrow type errors
    df = df.copy()
    str_cols = [
        'mut_num_id', 'cdna', 'hgvs_c', 'NM', 'std_cdna',
        'old_chr', 'old_chromosome', 'old_ref', 'old_alt', 'old_reference', 'old_alternative',
        'old_pos', 'old_position'
    ]
    for col in str_cols:
        if col in df.columns:
            def to_clean_str(val):
                if pd.isna(val):
                    return None
                val_str = str(val).strip()
                if val_str.endswith(".0"):
                    val_str = val_str[:-2]
                if val_str in ("nan", "None", ""):
                    return None
                return val_str
            df[col] = df[col].apply(to_clean_str)

    coord_cols = ['CHROM', 'POS', 'REF', 'ALT']
    
    # Separate rows with valid coordinates from those with missing coordinates
    # Rows with missing coordinates cannot be grouped and should be preserved as-is.
    df_valid = df[df[coord_cols].notnull().all(axis=1)].copy()
    df_invalid = df[~df[coord_cols].notnull().all(axis=1)].copy()
    
    if df_valid.empty:
        return df

    dup_mask = df_valid.duplicated(subset=coord_cols, keep=False)
    num_dups = dup_mask.sum()
    if num_dups > 0:
        print(f"De-duplicating {num_dups} duplicate coordinate rows in {'Oracle' if is_oracle else 'NDG'} on GRCh38...")

    agg_dict = {}
    for col in df_valid.columns:
        if col in coord_cols:
            continue
        
        col_lower = col.lower()
        
        # 1. Text descriptions, HGVS, transcripts, IDs, and hg19 chromosomes/references
        if col in ['mut_num_id', 'cdna', 'hgvs_c', 'NM', 'std_cdna', 'old_chr', 'old_chromosome', 'old_ref', 'old_alt', 'old_reference', 'old_alternative', 'abs_desc', 'phenotypes']:
            def agg_str_id(s):
                vals = []
                for val in s.dropna():
                    val_str = str(val).strip()
                    if val_str.endswith(".0"):
                        val_str = val_str[:-2]
                    vals.append(val_str)
                # Keep sorted and unique semicolon-separated values
                return ";".join(sorted(set(vals)))
            agg_dict[col] = agg_str_id
        
        # 2. hg19 positions (old_pos, old_position)
        elif col in ['old_pos', 'old_position']:
            def agg_pos(s):
                vals = []
                for val in s.dropna():
                    try:
                        vals.append(str(int(float(val))))
                    except ValueError:
                        vals.append(str(val))
                return ";".join(sorted(set(vals)))
            agg_dict[col] = agg_pos
            
        # 3. Counts and allele frequencies (max to prevent duplication/double-counting)
        elif any(x in col_lower for x in ['_ac', '_an', 'n_pac', 'n_pac_hq', 'n_het', 'n_hom', 'n_hemi', 'n_mos', 'cardio_n', 'ctrl_n', 'cardio_n_exome', 'ctrl_n_exome', 'cardio_n_panel', 'ctrl_n_panel', 'cardio_n_het', 'ctrl_n_het', 'cardio_n_hom', 'ctrl_n_hom', '_freq']):
            agg_dict[col] = 'max'
            
        # 4. Fallback for all other columns
        else:
            agg_dict[col] = 'first'
            
    # Group and aggregate
    df_deduped = df_valid.groupby(coord_cols, as_index=False).agg(agg_dict)
    
    # Recombine with coordinate-invalid rows
    return pd.concat([df_deduped, df_invalid], ignore_index=True)

def process_gene(gene, oracle_path, ndg_path, out_dir, only_liftover=False):
    """
    Loads raw cases and controls for a gene, renames coordinates, lifts over via GeneBe,
    performs outer join, adds SOURCE_ORIGIN column, and saves to Parquet.
    """
    print(f"\n--- Processing gene: {gene} ---")

    # Load datasets
    df_oracle = pd.read_parquet(oracle_path)
    df_ndg = pd.read_parquet(ndg_path)

    # Open reference FASTA file using pysam (supporting both local and cluster environments)
    fasta_paths = [
        "/home/mruizp/data_references/genomes/Homo_sapiens/GATK_bundle/v0/Homo_sapiens_assembly38.fasta",
        "/references/genomes/Homo_sapiens/GATK_bundle/v0/Homo_sapiens_assembly38.fasta"
    ]
    fasta_path = None
    for p in fasta_paths:
        if os.path.exists(p):
            fasta_path = p
            break

    fa = None
    if fasta_path:
        try:
            import pysam
            fa = pysam.FastaFile(fasta_path)
            print(f"Loaded GRCh38 reference FASTA from: {fasta_path}")
        except Exception as e:
            print(f"Warning: Failed to load FASTA index for {fasta_path}: {e}")
    else:
        print("Warning: GRCh38 reference FASTA file not found. Indel normalization will be skipped.")

    print(f"Loaded Oracle: {len(df_oracle)} rows | NDG: {len(df_ndg)} rows")

    # 1. Rename existing coordinate columns to old_* prefixes to prevent mixing
    # Oracle coordinate columns are: chr, pos, ref, alt
    oracle_rename = {}
    for col in ['chr', 'pos', 'ref', 'alt']:
        if col in df_oracle.columns:
            oracle_rename[col] = f"old_{col}"
    if oracle_rename:
        df_oracle.rename(columns=oracle_rename, inplace=True)
        print(f"Renamed Oracle coordinates: {list(oracle_rename.values())}")

    # NDG coordinate columns are: chromosome, position, reference, alternative
    ndg_rename = {}
    for col in ['chromosome', 'position', 'reference', 'alternative']:
        if col in df_ndg.columns:
            ndg_rename[col] = f"old_{col}"
    if ndg_rename:
        df_ndg.rename(columns=ndg_rename, inplace=True)
        print(f"Renamed NDG coordinates: {list(ndg_rename.values())}")

    # 2. Standardize cDNA identifiers
    # Oracle has cdna (e.g. NM_017644.3:c.747C>T)
    # NDG has NM and hgvs_c. Combine them as NM:hgvs_c
    df_oracle['std_cdna'] = df_oracle['cdna'].astype(str).str.strip()
    
    def get_ndg_cdna(row):
        nm = row.get('NM')
        hgvs = row.get('hgvs_c')
        if pd.isna(nm) or pd.isna(hgvs):
            return None
        nm_str = str(nm).strip()
        hgvs_str = str(hgvs).strip()
        if nm_str and hgvs_str:
            return f"{nm_str}:{hgvs_str}"
        return None

    df_ndg['std_cdna'] = df_ndg.apply(get_ndg_cdna, axis=1)

    # 3. Gather unique cDNA strings to query GeneBe
    oracle_cdnas = df_oracle['std_cdna'].dropna().unique().tolist()
    ndg_cdnas = df_ndg['std_cdna'].dropna().unique().tolist()
    all_unique_cdnas = list(set(oracle_cdnas + ndg_cdnas))
    
    # Filter to valid looking cDNA strings (must contain ':')
    valid_cdnas = [c for c in all_unique_cdnas if ":" in c]
    print(f"Unique cDNA candidates: {len(all_unique_cdnas)} | Valid cDNA formats: {len(valid_cdnas)}")

    # Query GeneBe
    if only_liftover:
        print("Bypassing GeneBe queries (--only-liftover enabled). Using coordinate liftover fallback for all variants.")
        cdna_to_coords = {}
    else:
        cdna_to_coords = query_genebe_batch(valid_cdnas)

    # 4. Map coordinates back to the original datasets
    def map_coords(df, is_oracle=True, fa=None):
        chroms, poses, refs, alts = [], [], [], []
        fallback_count = 0
        for idx, row in df.iterrows():
            cdna = row.get('std_cdna')
            c, p, r, a = None, None, None, None
            
            # Try cDNA-based GeneBe mapping first
            if pd.notna(cdna) and cdna in cdna_to_coords:
                c, p, r, a = cdna_to_coords[cdna]
                
            # Fallback to hg19 coordinate liftover if GeneBe mapping failed
            if c is None or p is None:
                if is_oracle:
                    old_chrom = row.get('old_chr')
                    old_pos = row.get('old_pos')
                    old_ref = row.get('old_ref')
                    old_alt = row.get('old_alt')
                    abs_desc = row.get('abs_desc')
                else:
                    old_chrom = row.get('old_chromosome')
                    old_pos = row.get('old_position')
                    old_ref = row.get('old_reference')
                    old_alt = row.get('old_alternative')
                    abs_desc = None
                
                is_indel = False
                if pd.isna(old_ref) or str(old_ref).strip() == '-' or pd.isna(old_alt) or str(old_alt).strip() == '-':
                    is_indel = True
                    
                if is_indel and pd.notna(old_chrom) and pd.notna(old_pos):
                    target_chrom, target_pos, target_ref, target_alt = normalize_indel_hg38(
                        old_chrom, old_pos, old_ref, old_alt, abs_desc, lo, fa
                    )
                    if target_chrom and target_pos:
                        c, p, r, a = target_chrom, target_pos, target_ref, target_alt
                        fallback_count += 1
                elif pd.notna(old_chrom) and pd.notna(old_pos):
                    target_chrom, target_pos = liftover_coordinate_hg19_to_hg38(old_chrom, old_pos)
                    if target_chrom and target_pos:
                        c = target_chrom
                        p = target_pos
                        r = old_ref if pd.notna(old_ref) else None
                        a = old_alt if pd.notna(old_alt) else None
                        fallback_count += 1
            
            chroms.append(c)
            poses.append(p)
            refs.append(r)
            alts.append(a)
            
        df['CHROM'] = chroms
        df['POS'] = poses
        df['REF'] = refs
        df['ALT'] = alts
        # Coerce POS to Int64 to handle nulls gracefully
        df['POS'] = df['POS'].astype('Int64')
        return fallback_count

    oracle_fallback = map_coords(df_oracle, is_oracle=True, fa=fa)
    ndg_fallback = map_coords(df_ndg, is_oracle=False, fa=fa)

    if fa:
        fa.close()

    # Count successful conversions
    oracle_mapped = df_oracle['CHROM'].notna().sum()
    ndg_mapped = df_ndg['CHROM'].notna().sum()
    print(f"Mapped coordinates - Oracle: {oracle_mapped}/{len(df_oracle)} (including {oracle_fallback} via hg19 fallback) | NDG: {ndg_mapped}/{len(df_ndg)} (including {ndg_fallback} via hg19 fallback)")

    # Log variants that failed conversion
    oracle_failed = df_oracle[df_oracle['CHROM'].isnull() & df_oracle['std_cdna'].notna()]
    if not oracle_failed.empty:
        print(f"Warning: {len(oracle_failed)} Oracle variants failed cDNA conversion.")
        if len(oracle_failed) <= 5:
            print("Failed: ", oracle_failed['std_cdna'].tolist())

    # De-duplicate datasets prior to genomic join on GRCh38 coordinates
    oracle_before_dedup = len(df_oracle)
    ndg_before_dedup = len(df_ndg)

    df_oracle = deduplicate_dataset(df_oracle, is_oracle=True)
    df_ndg = deduplicate_dataset(df_ndg, is_oracle=False)

    print(f"Row count comparison (Pre -> Post de-duplication):")
    print(f"  Oracle: {oracle_before_dedup} -> {len(df_oracle)} rows")
    print(f"  NDG   : {ndg_before_dedup} -> {len(df_ndg)} rows")

    # Add source indicators before merge
    df_oracle['in_oracle'] = True
    df_ndg['in_ndg'] = True

    # 5. Perform genomic outer join on GRCh38 coordinates
    # For variants with null coordinates, they cannot join and should be kept as unique rows.
    # To handle this properly:
    # Split both datasets into coordinate-valid and coordinate-invalid subsets.
    valid_oracle = df_oracle[df_oracle['CHROM'].notna()]
    invalid_oracle = df_oracle[df_oracle['CHROM'].isnull()]

    valid_ndg = df_ndg[df_ndg['CHROM'].notna()]
    invalid_ndg = df_ndg[df_ndg['CHROM'].isnull()]

    # Join the valid subsets
    df_joined_valid = pd.merge(
        valid_oracle,
        valid_ndg,
        on=['CHROM', 'POS', 'REF', 'ALT'],
        how='outer',
        suffixes=('_oracle', '_ndg')
    )

    # Concatenate the coordinate-invalid subsets (which remain unique to their sources)
    df_combined = pd.concat([df_joined_valid, invalid_oracle, invalid_ndg], ignore_index=True)
    print(f"Combined dataset final size: {len(df_combined)} rows")

    # Resolve common columns
    df_combined['gene'] = df_combined['gene_oracle'].combine_first(df_combined['gene_ndg'])
    df_combined.drop(columns=['gene_oracle', 'gene_ndg'], inplace=True, errors='ignore')
    
    # Fill in_oracle / in_ndg nulls with False
    df_combined['in_oracle'] = df_combined['in_oracle'].fillna(False)
    df_combined['in_ndg'] = df_combined['in_ndg'].fillna(False)

    # 6. Annotate SOURCE_ORIGIN column based on join operation result
    def get_source_origin(row):
        in_o = row['in_oracle']
        in_n = row['in_ndg']
        if in_o and in_n:
            return "HIC_CASES|NDG_CONTROLS"
        elif in_o:
            return "HIC_CASES"
        elif in_n:
            return "NDG_CONTROLS"
        return "UNKNOWN"

    df_combined['SOURCE_ORIGIN'] = df_combined.apply(get_source_origin, axis=1)

    # Clean up indicator columns
    df_combined.drop(columns=['in_oracle', 'in_ndg'], inplace=True)

    # 7. Collect and log unmapped variants (unable to convert to hg38)
    unmapped_records = []
    if not invalid_oracle.empty:
        for _, row in invalid_oracle.iterrows():
            cdna = row.get('std_cdna')
            nm = cdna.split(':')[0] if isinstance(cdna, str) and ':' in cdna else None
            unmapped_records.append({
                'gene': gene,
                'std_cdna': cdna,
                'refseq_nm': nm,
                'source': 'HIC_CASES',
                'old_chrom': row.get('old_chr'),
                'old_pos': row.get('old_pos'),
                'old_ref': row.get('old_ref'),
                'old_alt': row.get('old_alt')
            })
            
    if not invalid_ndg.empty:
        for _, row in invalid_ndg.iterrows():
            cdna = row.get('std_cdna')
            nm = row.get('NM') if pd.notna(row.get('NM')) else (cdna.split(':')[0] if isinstance(cdna, str) and ':' in cdna else None)
            unmapped_records.append({
                'gene': gene,
                'std_cdna': cdna,
                'refseq_nm': nm,
                'source': 'NDG_CONTROLS',
                'old_chrom': row.get('old_chromosome'),
                'old_pos': row.get('old_position'),
                'old_ref': row.get('old_reference'),
                'old_alt': row.get('old_alternative')
            })

    log_dir = os.path.join(out_dir, "log")
    unmapped_csv = os.path.join(log_dir, f"{gene}_unmapped.csv")
    if unmapped_records:
        df_unmapped = pd.DataFrame(unmapped_records)
        os.makedirs(log_dir, exist_ok=True)
        df_unmapped.to_csv(unmapped_csv, index=False)
        print(f"Saved {len(df_unmapped)} unmapped variants to {unmapped_csv}")
    else:
        # Move any stale unmapped file from previous runs to the trash folder
        if os.path.exists(unmapped_csv):
            trash_dir = os.path.join(out_dir, "trash")
            os.makedirs(trash_dir, exist_ok=True)
            import shutil
            shutil.move(unmapped_csv, os.path.join(trash_dir, f"{gene}_unmapped.csv"))
            print(f"Moved stale unmapped file to {os.path.join(trash_dir, f'{gene}_unmapped.csv')}")

    # 8. Serialize to parquet (.parquet only)
    os.makedirs(out_dir, exist_ok=True)
    out_parquet = os.path.join(out_dir, f"{gene}.parquet")

    save_with_retry(df_combined, out_parquet)

    print(f"Serialized output: {len(df_combined)} rows saved to:")
    print(f"  {out_parquet}")
    print(f"  SOURCE_ORIGIN counts: {df_combined['SOURCE_ORIGIN'].value_counts().to_dict()}")

import argparse

def main():
    parser = argparse.ArgumentParser(description="Standardize, liftover, and combine Valencia cases and NDG controls.")
    parser.add_argument(
        "--test",
        action="store_true",
        help="Run in test mode (only process the small gene TNNC1 and output to a combined_test directory)"
    )
    parser.add_argument(
        "--gene",
        type=str,
        default=None,
        help="Process only a specific gene (e.g. TNNC1, BAG3)"
    )
    parser.add_argument(
        "--oracle-dir",
        type=str,
        default="_RAW/data_Valencia_plus_hic/oracle",
        help="Directory containing Oracle raw files"
    )
    parser.add_argument(
        "--ndg-dir",
        type=str,
        default="_RAW/data_Valencia_plus_hic/ndg",
        help="Directory containing NDG raw files"
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default=None,
        help="Output directory (defaults to combined or combined_test if in test mode)"
    )
    parser.add_argument(
        "--overwrite-all",
        action="store_true",
        help="Overwrite existing combined files (by default, existing files are skipped)"
    )
    parser.add_argument(
        "--only-liftover",
        action="store_true",
        help="Bypass GeneBe API queries entirely and perform local coordinate liftover for all variants"
    )
    args = parser.parse_args()

    oracle_dir = args.oracle_dir
    ndg_dir = args.ndg_dir

    if args.out_dir:
        out_dir = args.out_dir
    elif args.test:
        out_dir = "_RAW/data_Valencia_plus_hic/combined_test"
    else:
        out_dir = "_RAW/data_Valencia_plus_hic/combined"

    if not os.path.exists(oracle_dir) or not os.path.exists(ndg_dir):
        print(f"Error: Oracle directory ({oracle_dir}) or NDG directory ({ndg_dir}) does not exist.")
        sys.exit(1)

    # Determine files to process from both folders
    oracle_files = glob.glob(os.path.join(oracle_dir, "*.parquet"))
    ndg_files = glob.glob(os.path.join(ndg_dir, "*.parquet"))

    oracle_genes = {os.path.splitext(os.path.basename(f))[0] for f in oracle_files}
    ndg_genes = {os.path.splitext(os.path.basename(f))[0] for f in ndg_files}

    all_genes = sorted(list(oracle_genes.union(ndg_genes)))
    if not all_genes:
        print("Error: No parquet files found in either Oracle or NDG directory.")
        sys.exit(1)

    # Filter by gene selection (either --gene, or TNNC1 as default for --test)
    selected_gene = args.gene
    if not selected_gene and args.test:
        selected_gene = "TNNC1"
        print(f"Test mode enabled: defaulting to process single small gene '{selected_gene}' in '{out_dir}'")

    processed_count = 0
    for gene in all_genes:
        # If filtering by gene, skip others
        if selected_gene and gene != selected_gene:
            continue

        o_path = os.path.join(oracle_dir, f"{gene}.parquet")
        n_path = os.path.join(ndg_dir, f"{gene}.parquet")

        o_exists = os.path.exists(o_path)
        n_exists = os.path.exists(n_path)

        if not o_exists:
            print(f"Warning: Gene {gene} is missing from Oracle cases directory. Skipping combination.")
            continue
        if not n_exists:
            print(f"Warning: Gene {gene} is missing from NDG controls directory. Skipping combination.")
            continue

        if not args.overwrite_all:
            out_parquet = os.path.join(out_dir, f"{gene}.parquet")
            if os.path.exists(out_parquet) and check_manifest_match(gene, o_path, n_path, out_dir):
                print(f"Output file for gene {gene} already exists and inputs are unchanged. Skipping combination.")
                continue

        process_gene(gene, o_path, n_path, out_dir, only_liftover=args.only_liftover)
        update_manifest(gene, o_path, n_path, out_dir)
        processed_count += 1

    if processed_count == 0:
        print("Warning: No genes were processed. Please check your gene selection or test settings.")
    else:
        print("\n==================================================")
        print("Integration and liftover execution finished successfully!")
        print(f"Output saved to: {out_dir}")
        print("==================================================")

if __name__ == "__main__":
    main()
