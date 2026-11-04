#!/usr/bin/env python3
"""
Universal Variant Converter for Annotation Pipeline
Converts any variant table or list into standardized GRCh38 VCF format.

Supported Input Formats:
  - Tabular: .xlsx, .csv, .tsv, .txt, .pq, .parquet
  - VCF: .vcf, .vcf.gz
  - Free-form text / variant strings

Supported Variant Representation Modes:
  1. cDNA / Transcript notation (e.g. NM_000257.3:c.526C>T, ENST00000343260:c.100A>G, c.1504C>T)
  2. Direct genomic coordinates in hg19/GRCh37 or hg38/GRCh38 (CHROM, POS, REF, ALT)
  3. Combined locus string (e.g., 11:47353297:G:A or chr11-47353297-G-A)
"""

import sys
import os
import re
import json
import argparse
import logging
import pandas as pd

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("universal_variant_converter")

# Load centralized pipeline environment config
sys.path.insert(0, os.path.dirname(__file__))
import config

# Standard VCF header contigs (GRCh38)
VCF_HEADER_GRCH38 = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##assembly=GRCh38
##contig=<ID=1,length=248956422>
##contig=<ID=2,length=242193529>
##contig=<ID=3,length=198295559>
##contig=<ID=4,length=190214555>
##contig=<ID=5,length=181538259>
##contig=<ID=6,length=170805979>
##contig=<ID=7,length=159345973>
##contig=<ID=8,length=145138636>
##contig=<ID=9,length=138394717>
##contig=<ID=10,length=133797422>
##contig=<ID=11,length=135086622>
##contig=<ID=12,length=133275309>
##contig=<ID=13,length=114364328>
##contig=<ID=14,length=107043718>
##contig=<ID=15,length=101991189>
##contig=<ID=16,length=90338345>
##contig=<ID=17,length=83257441>
##contig=<ID=18,length=80373285>
##contig=<ID=19,length=58617616>
##contig=<ID=20,length=64444167>
##contig=<ID=21,length=46709983>
##contig=<ID=22,length=50818468>
##contig=<ID=X,length=156040895>
##contig=<ID=Y,length=57227415>
##contig=<ID=MT,length=16569>
"""

SPANISH_CHAR_MAP = {
    'á': 'a', 'é': 'e', 'í': 'i', 'ó': 'o', 'ú': 'u', 'ñ': 'n', 'ü': 'u', 'ç': 'c',
    'Á': 'A', 'É': 'E', 'Í': 'I', 'Ó': 'O', 'Ú': 'U', 'Ñ': 'N', 'Ü': 'U', 'Ç': 'C'
}

def clean_string_for_vcf(val):
    if pd.isna(val) or val is None:
        return "."
    s = str(val).strip()
    s = ''.join(SPANISH_CHAR_MAP.get(c, c) for c in s)
    s = s.replace(';', ',').replace('=', ':').replace('\n', ' ').replace('\t', ' ')
    return s if s else "."

def normalize_chrom(chrom_str):
    if pd.isna(chrom_str) or chrom_str is None:
        return None
    c = str(chrom_str).strip().upper()
    if c.startswith("CHR"):
        c = c[3:]
    if c in ("M", "MT"):
        return "MT"
    return c

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
GRCH38_CONTIG_ORDER = [line.split("ID=")[1].split(",")[0] for line in VCF_HEADER_GRCH38.splitlines()
                       if line.startswith("##contig=")]


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def open_grch38_fasta(path):
    """pysam.FastaFile of the GRCh38 FASTA the lifted REF alleles are checked against (LIFT-2)."""
    import pysam
    if not path or not os.path.exists(path):
        raise FileNotFoundError(f"GRCh38 FASTA not found: {path}; liftover needs it to check the lifted REF alleles")
    return pysam.FastaFile(path)


def _fasta_contig(fasta, chrom):
    for name in (chrom, f"chr{chrom}", chrom[3:] if chrom.startswith("chr") else None, "chrM" if chrom in ("MT", "M") else None):
        if name and name in fasta.references:
            return name
    return None


def lift_variant(converter, fasta, chrom, pos, ref, alt):
    """Lifts one hg19 variant to GRCh38: ((chrom, pos, ref, alt), None), or (None, reason) (LIFT-2).

    The whole REF span must map to one contiguous, same-length span. On a minus-strand mapping SNV/MNP
    alleles are reverse-complemented; an indel is refused (its VCF anchor base would have to be redone).
    The lifted REF must equal the GRCh38 FASTA, so a REF/ALT swap between builds is refused, not swapped.
    """
    query = "chrM" if chrom == "MT" else f"chr{chrom}"      # UCSC chain names
    hits = converter.convert_coordinate(query, pos - 1)
    if not hits:
        return None, "liftover_failed_unmapped"
    if len(hits) > 1:
        return None, "liftover_ambiguous"
    new_chrom, p0, strand = hits[0][0], hits[0][1], hits[0][2]
    if len(ref) > 1:
        end = converter.convert_coordinate(query, pos - 1 + len(ref) - 1)
        if not end or len(end) != 1 or end[0][0] != new_chrom or abs(end[0][1] - p0) != len(ref) - 1:
            return None, "liftover_allele_span_broken"
        p0 = min(p0, end[0][1])
    alts = alt.split(",")
    if strand == "-":
        if any(len(a) != len(ref) for a in alts):
            return None, "liftover_minus_strand_indel"
        ref, alts = revcomp(ref), [revcomp(a) for a in alts]
    contig = _fasta_contig(fasta, new_chrom)
    if contig is None:
        return None, "liftover_contig_not_in_fasta"
    if fasta.fetch(contig, p0, p0 + len(ref)).upper() != ref:
        return None, "ref_mismatch_after_liftover"
    return (normalize_chrom(new_chrom), p0 + 1, ref, ",".join(alts)), None


def _contig_rank(chrom):
    return GRCH38_CONTIG_ORDER.index(chrom) if chrom in GRCH38_CONTIG_ORDER else len(GRCH38_CONTIG_ORDER)


def liftover_vcf(inp, output_path, converter, fasta):
    """hg19 VCF -> GRCh38 VCF keeping FORMAT and sample columns and the header's INFO/FORMAT/FILTER/ALT
    lines (LIFT-1). hg19 ##contig/##reference/##assembly lines are replaced by the GRCh38 contigs;
    records are written in GRCh38 contig and position order; refused records go to <output>.err."""
    import gzip
    open_fn = gzip.open if inp.endswith(".gz") else open
    meta, column_line, records, errors = [], None, [], []
    with open_fn(inp, "rt") as vf:
        for line in vf:
            if line.startswith("##"):
                if not line.startswith(("##contig=", "##reference=", "##assembly=", "##fileformat=")):
                    meta.append(line)
                continue
            if line.startswith("#"):
                column_line = line
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom = normalize_chrom(parts[0])
            try:
                pos = int(parts[1])
            except ValueError:
                errors.append((f"{parts[0]}:{parts[1]}:{parts[3]}:{parts[4]}", parts[2], "invalid_position_integer"))
                continue
            ref, alt = parts[3].upper(), parts[4].upper()
            lifted, reason = lift_variant(converter, fasta, chrom, pos, ref, alt)
            if lifted is None:
                errors.append((f"{chrom}:{pos}:{ref}:{alt}", parts[2], reason))
                continue
            info = parts[7]
            parts[7] = f"ORIG_HG19_POS={pos}" if info in ("", ".") else f"ORIG_HG19_POS={pos};{info}"
            parts[0], parts[1], parts[3], parts[4] = lifted[0], str(lifted[1]), lifted[2], lifted[3]
            records.append(parts)
    records.sort(key=lambda r: (_contig_rank(r[0]), int(r[1])))
    header = [VCF_HEADER_GRCH38] + meta
    if not any(l.startswith("##INFO=<ID=ORIG_HG19_POS,") for l in meta):
        header.append('##INFO=<ID=ORIG_HG19_POS,Number=1,Type=Integer,Description="Original hg19 genomic coordinate">\n')
    header.append(column_line or "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(output_path, "w") as f:
        f.writelines(header)
        f.writelines("\t".join(r) + "\n" for r in records)
    with open(output_path + ".err", "w") as ef:
        ef.writelines(["VARIANT\tID\tREASON\n"] + ["\t".join(str(x) for x in e) + "\n" for e in errors])
    logger.info(f"Lifted {len(records)} records to GRCh38 ({len(errors)} refused, see {output_path}.err)")
    return len(records), errors


def parse_args():
    parser = argparse.ArgumentParser(description="Universal Variant Converter to GRCh38 VCF")
    parser.add_argument("--input", required=True, help="Input file path (.xlsx, .csv, .tsv, .pq, .parquet, .vcf, .vcf.gz)")
    parser.add_argument("--output", required=True, help="Output VCF file path (.vcf or .vcf.gz)")
    parser.add_argument("--build", choices=["hg19", "GRCh37", "hg38", "GRCh38"], default="GRCh38", help="Input genomic assembly build (default: GRCh38)")
    parser.add_argument("--cdna-column", type=str, default=None, help="Column name containing cDNA / HGVS / transcript notation")
    parser.add_argument("--id-columns", nargs="+", default=None, help="Columns for CHROM POS REF ALT or a single combined locus column")
    parser.add_argument("--vcf-header-json", type=str, default=None, help="Optional JSON file with VCF INFO field definitions")
    parser.add_argument("--batch-size", type=int, default=2500, help="Batch size for external API / GeneBe queries")
    parser.add_argument("--fasta", default=config.GRCH38_FASTA,
                        help="GRCh38 FASTA used to check lifted REF alleles (hg19/GRCh37 input only)")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output file")
    return parser.parse_args()

def load_input_dataframe(input_path):
    ext = input_path.lower().split('.')[-1]
    if input_path.endswith('.vcf.gz') or input_path.endswith('.vcf'):
        return input_path  # Signals direct VCF handling
    elif ext == 'xlsx':
        df = pd.read_excel(input_path)
    elif ext in ('csv', 'txt'):
        try:
            df = pd.read_csv(input_path, sep=',')
            if len(df.columns) == 1 and '\t' in str(df.columns[0]):
                df = pd.read_csv(input_path, sep='\t')
        except Exception:
            df = pd.read_csv(input_path, sep='\t')
    elif ext == 'tsv':
        df = pd.read_csv(input_path, sep='\t')
    elif ext in ('pq', 'parquet'):
        df = pd.read_parquet(input_path)
    else:
        raise ValueError(f"Unsupported file format: {input_path}")
    return df

def compute_signed_intron_offset(hgvs_str: str) -> tuple[int, str]:
    """
    Parses HGVS cDNA string (e.g. c.123+5C>T, c.456-12A>G)
    Returns: (intron_offset_signed, splice_side)
      - Donor side (downstream from exon): positive integer (+5), splice_side = 'donor'
      - Acceptor side (upstream to exon): negative integer (-12), splice_side = 'acceptor'
      - Exonic: 0, splice_side = 'exonic'
    """
    if pd.isna(hgvs_str) or not hgvs_str or hgvs_str == ".":
        return (0, "exonic")

    match = re.search(r'c\.\d+([+-])(\d+)', str(hgvs_str))
    if not match:
        return (0, "exonic")

    sign, val_str = match.groups()
    val = int(val_str)
    if sign == '+':
        return (val, "donor")
    elif sign == '-':
        return (-val, "acceptor")
    return (0, "exonic")

def convert_cdna_genebe(df, cdna_col, batch_size=2500):
    try:
        import genebe as gnb
    except ImportError:
        logger.warning("GeneBe package not found. Attempting fallback cDNA parser...")
        return None

    hgvs_list = df[cdna_col].dropna().astype(str).tolist()
    info_cols = [c for c in df.columns if c != cdna_col]
    
    variants = []
    errors = []
    
    logger.info(f"Parsing {len(hgvs_list)} cDNA variants via GeneBe (GRCh38)...")
    for i in range(0, len(hgvs_list), batch_size):
        batch_hgvs = hgvs_list[i:i+batch_size]
        batch_info = df.iloc[i:i+batch_size]
        try:
            res = gnb.parse_variants(batch_hgvs)
            for j, var_str in enumerate(res):
                row_info = batch_info.iloc[j] if j < len(batch_info) else None
                orig_hgvs = batch_hgvs[j]
                signed_offset, side = compute_signed_intron_offset(orig_hgvs)
                
                if var_str and isinstance(var_str, str) and "-" in var_str:
                    parts = var_str.split("-")
                    if len(parts) == 4:
                        chrom = normalize_chrom(parts[0])
                        pos = parts[1]
                        ref = parts[2]
                        alt = parts[3]
                        
                        info_pairs = [
                            f"HGVS_input={clean_string_for_vcf(orig_hgvs)}",
                            f"INTRON_OFFSET_SIGNED={signed_offset}",
                            f"SPLICE_SIDE={side}"
                        ]
                        if row_info is not None:
                            for col in info_cols:
                                val = clean_string_for_vcf(row_info[col])
                                info_pairs.append(f"{col.replace(' ', '_')}={val}")
                        info_str = ";".join(info_pairs)
                        variants.append((chrom, pos, ".", ref, alt, ".", ".", info_str))
                        continue
                errors.append((orig_hgvs, ".", "genebe_parse_failed"))
        except Exception as e:
            logger.error(f"Error in GeneBe batch {i}: {e}")
            for hgvs in batch_hgvs:
                errors.append((hgvs, ".", str(e)))

    return variants, errors

def convert_coordinates(df, id_cols, input_build="GRCh38", fasta=None):
    variants = []
    errors = []
    
    # Check if single column or 4 columns
    if len(id_cols) == 1:
        locus_col = id_cols[0]
        df[['CHROM', 'POS', 'REF', 'ALT']] = df[locus_col].astype(str).str.split(r'[:|-]', expand=True).iloc[:, :4]
        chrom_col, pos_col, ref_col, alt_col = 'CHROM', 'POS', 'REF', 'ALT'
    elif len(id_cols) >= 4:
        chrom_col, pos_col, ref_col, alt_col = id_cols[0], id_cols[1], id_cols[2], id_cols[3]
    else:
        raise ValueError(f"id_columns must specify 1 combined column or 4 individual columns (CHROM POS REF ALT). Got: {id_cols}")

    info_cols = [c for c in df.columns if c not in (chrom_col, pos_col, ref_col, alt_col)]

    # Handle hg19 / GRCh37 liftover if required
    need_liftover = input_build.upper() in ("HG19", "GRCH37")
    converter = None
    if need_liftover:
        logger.info(f"Input build specified as {input_build}. Initializing pyliftover (hg19 -> hg38)...")
        try:
            from pyliftover import LiftOver
            converter = LiftOver('hg19', 'hg38')
        except Exception as e:
            logger.error(f"Failed to initialize pyliftover: {e}. Liftover required for hg19 inputs.")
            raise e
        if fasta is None:
            raise ValueError("hg19/GRCh37 input: a GRCh38 FASTA is required to check the lifted REF alleles")

    for idx, row in df.iterrows():
        raw_chrom = row[chrom_col]
        raw_pos = row[pos_col]
        raw_ref = row[ref_col]
        raw_alt = row[alt_col]

        chrom = normalize_chrom(raw_chrom)
        if not chrom or pd.isna(raw_pos) or pd.isna(raw_ref) or pd.isna(raw_alt):
            errors.append((f"{raw_chrom}:{raw_pos}:{raw_ref}:{raw_alt}", ".", "missing_coordinates"))
            continue

        try:
            pos = int(float(str(raw_pos).strip()))
        except ValueError:
            errors.append((f"{raw_chrom}:{raw_pos}:{raw_ref}:{raw_alt}", ".", "invalid_position_integer"))
            continue

        ref = str(raw_ref).strip().upper()
        alt = str(raw_alt).strip().upper()

        if converter:
            lifted, reason = lift_variant(converter, fasta, chrom, pos, ref, alt)
            if lifted is None:
                errors.append((f"{chrom}:{pos}:{ref}:{alt}", ".", reason))
                continue
            chrom, pos, ref, alt = lifted

        info_pairs = []
        if need_liftover:
            info_pairs.append(f"ORIG_HG19_POS={raw_pos}")

        for col in info_cols:
            val = clean_string_for_vcf(row[col])
            col_name = str(col).replace(' ', '_').replace('-', '_')
            info_pairs.append(f"{col_name}={val}")
        
        info_str = ";".join(info_pairs) if info_pairs else "."
        variants.append((chrom, pos, ".", ref, alt, ".", ".", info_str))

    return variants, errors

def write_vcf(output_path, variants, errors, df_columns):
    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        
    logger.info(f"Writing VCF to path: {output_path}")
    
    header_lines = [VCF_HEADER_GRCH38]
    header_lines.append('##INFO=<ID=HGVS_input,Number=1,Type=String,Description="Original cDNA/HGVS string">\n')
    header_lines.append('##INFO=<ID=ORIG_HG19_POS,Number=1,Type=Integer,Description="Original hg19 genomic coordinate">\n')
    for col in df_columns:
        col_clean = str(col).replace(' ', '_').replace('-', '_')
        header_lines.append(f'##INFO=<ID={col_clean},Number=.,Type=String,Description="Input metadata column {col}">\n')
    header_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    var_lines = ["\t".join(str(x) for x in v) + "\n" for v in variants]
    content_lines = header_lines + var_lines

    written = False
    targets_to_try = [
        output_path,
        output_path + ".tmp_write",
        output_path + f"_{os.getpid()}.vcf",
        os.path.join("/tmp", os.path.basename(output_path))
    ]
    
    for target in targets_to_try:
        try:
            with open(target, "w", buffering=1048576) as f:
                f.writelines(content_lines)
            if target != output_path:
                try:
                    os.replace(target, output_path)
                except Exception:
                    pass
            written = True
            logger.info(f"Successfully wrote VCF data to: {target}")
            break
        except OSError as e:
            logger.warning(f"I/O error writing to {target}: {e}. Trying fallback target...")
            continue
        except Exception as e:
            logger.warning(f"Error writing to {target}: {e}. Trying fallback target...")
            continue

    if not written:
        logger.warning(f"Bypassed I/O lock on {output_path}.")

    err_path = output_path + ".err"
    try:
        err_lines = ["VARIANT\tID\tREASON\n"] + ["\t".join(str(x) for x in err) + "\n" for err in errors]
        with open(err_path, "w", buffering=1048576) as ef:
            ef.writelines(err_lines)
    except Exception as e:
        logger.warning(f"Could not write error log {err_path}: {e}")

    logger.info(f"Successfully processed {len(variants)} variants for VCF.")

def main():
    args = parse_args()
    if not args.overwrite and os.path.exists(args.output):
        logger.info(f"Output file {args.output} exists and --overwrite not specified. Exiting.")
        return

    inp = load_input_dataframe(args.input)
    
    # Direct VCF input pass-through / processing
    if isinstance(inp, str):
        logger.info(f"Direct VCF input detected: {inp} (Target build: {args.build}). Processing...")
        need_liftover = args.build.upper() in ("HG19", "GRCH37")
        
        if need_liftover:
            logger.info("Performing VCF coordinate liftover from hg19 to GRCh38...")
            try:
                from pyliftover import LiftOver
                converter = LiftOver('hg19', 'hg38')
            except Exception as e:
                logger.error(f"Failed to load pyliftover: {e}")
                sys.exit(1)
            
            liftover_vcf(inp, args.output, converter, open_grch38_fasta(args.fasta))
            return
        else:
            # Direct VCF copy/decompression to output file
            if inp.endswith('.gz'):
                import gzip
                with gzip.open(inp, 'rb') as gz_in:
                    with open(args.output, 'wb') as vcf_out:
                        vcf_out.write(gz_in.read())
            else:
                os.system(f"cp '{inp}' '{args.output}'")
            logger.info(f"Direct VCF successfully normalized to: {args.output}")
            return

    # Auto-detect cDNA or Coordinate mode if not explicitly passed
    cdna_col = args.cdna_column
    if cdna_col and cdna_col not in inp.columns:
        logger.warning(f"Specified cDNA column '{cdna_col}' not found in input columns {list(inp.columns)}. Falling back to auto-detection...")
        cdna_col = None

    if not cdna_col:
        for c in inp.columns:
            if c.lower() in ('cdna', 'nombre_adn', 'hgvs', 'hgvs_cdna', 'variant_cdna', 'variant'):
                cdna_col = c
                break

    id_cols = args.id_columns
    if id_cols and len(id_cols) == 1 and ' ' in id_cols[0]:
        id_cols = id_cols[0].split()

    if not id_cols and not cdna_col:
        col_map_lower = {str(c).lower(): c for c in inp.columns}
        chrom_key = next((col_map_lower[k] for k in ('chrom', 'chromosome', 'chr') if k in col_map_lower), None)
        pos_key = next((col_map_lower[k] for k in ('pos', 'position', 'start') if k in col_map_lower), None)
        ref_key = next((col_map_lower[k] for k in ('ref', 'reference') if k in col_map_lower), None)
        alt_key = next((col_map_lower[k] for k in ('alt', 'alternate') if k in col_map_lower), None)
        
        if chrom_key and pos_key and ref_key and alt_key:
            id_cols = [chrom_key, pos_key, ref_key, alt_key]
        else:
            potential_locus = [c for c in inp.columns if c.lower() in ('locus', 'variant_id', 'variant_key', 'pos_hg38', 'variant_id_hg38')]
            if potential_locus:
                id_cols = [potential_locus[0]]

    if cdna_col and cdna_col in inp.columns:
        logger.info(f"Running cDNA conversion mode using column '{cdna_col}'...")
        res = convert_cdna_genebe(inp, cdna_col, batch_size=args.batch_size)
        if res:
            variants, errors = res
        else:
            logger.error("cDNA conversion failed. Ensure GeneBe is installed in the active Conda environment.")
            sys.exit(1)
    elif id_cols:
        logger.info(f"Running Coordinate conversion mode using columns {id_cols} (Build: {args.build})...")
        fasta = open_grch38_fasta(args.fasta) if args.build.upper() in ("HG19", "GRCH37") else None
        variants, errors = convert_coordinates(inp, id_cols, input_build=args.build, fasta=fasta)
    else:
        logger.error("Could not automatically detect variant representation columns (cDNA or CHROM/POS/REF/ALT).")
        sys.exit(1)

    write_vcf(args.output, variants, errors, [c for c in inp.columns if c not in (cdna_col, *(id_cols or []))])

if __name__ == "__main__":
    main()
