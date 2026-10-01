"""
Pedigree Parsing, Trio Compound Het TRANS Phasing, and 5-Model Single-Allele Inheritance Engine.
"""
import json
import numpy as np
import pandas as pd
from . import logger


def _normalize_gt(val):
    s = str(val).upper().strip()
    if s in ["HET", "0/1", "1/0", "0|1", "1|0"]: return "HET"
    if s in ["HOMALT", "1/1", "1|1"]: return "HOMALT"
    if s in ["HOMREF", "0/0", "0|0"]: return "HOMREF"
    return "MISSING"


def read_pedigree_file(genotypes: str, group: str = None) -> dict[str, str]|pd.DataFrame:
    """Read a pedigree file (.ped, .json, or dict string) and return sample information."""
    if genotypes.endswith(".json"):
        if not group:
            logger.error(f"Group not provided for {genotypes}")
            raise ValueError(f"Group not provided for {genotypes}")
        with open(genotypes, "r") as f:
            gt_file = json.load(f)
            try:
                samples = {
                    s: [info["genotype"], info["relationship"]] for s, info in gt_file[group]["samples"].items()
                }
            except ValueError as e:
                logger.error(f"Group {group} not found in {genotypes}")
                raise e
    elif genotypes.endswith(".ped"):
        logger.info(f"Reading pedigree file {genotypes}")
        samples = pd.read_csv(genotypes, sep=r"\s+|\t", engine="python")
        if "FamilyID" not in samples.columns and "IndividualID" not in samples.columns:
            if len(samples.columns) >= 6:
                samples.columns = ["FamilyID", "IndividualID", "PaternalID", "MaternalID", "Sex", "Phenotype"] + list(samples.columns[6:])
        mask = (
            ~samples["FamilyID"].astype(str).str.lower().isin(["familyid", "family_id", "#familyid"])
            & ~samples["IndividualID"].astype(str).str.lower().isin(["individualid", "individual_id", "#individualid"])
            & ~samples["FamilyID"].astype(str).str.startswith("#")
        )
        samples = samples[mask].copy()
        logger.debug(f"Loaded clean pedigree samples: {len(samples)}")
    else:
        try:
            samples = dict(item.split(":") for item in genotypes.split(","))
        except Exception:
            logger.error(
                f"Invalid format for --genotypes {genotypes}, should be a dictionary-like string"
            )
            raise ValueError(
                "Invalid format for genotypes, should be a dictionary-like string"
            )
            
    return samples


def filter_phenotype(data: pd.DataFrame, pedigree: pd.DataFrame):
    """Filter dataframe for rows that have positive phenotype in pedigree file."""
    logger.info("Filtering dataframe by phenotype positives defined in pedigree file")
    if "Phenotype" not in pedigree.columns:
        logger.warning("Column phenotype not found in data, skipping filtering")
        return data
    positives = pedigree[pedigree["Phenotype"] == 1]["IndividualID"].tolist()
    if not positives:
        logger.warning("No individuals with phenotype 1 found in pedigree, skipping filtering")
        return data
    positives = ['GT_' + str(p) for p in positives]
    data = data[data[positives].apply(lambda x: x.isin(["HET", "HOMALT"]).all(), axis=1)]
    logger.info(f"Dataframe shape after phenotype filtering: {data.shape}")
    return data


def filter_by_dominant(data: pd.DataFrame, samples: list[str]):
    return data


def _mark_pairwise_trans(data: pd.DataFrame, proband_id: str, paternal_id: str, maternal_id: str):
    """Pairwise TRANS matching for a specific trio (Proband, Father, Mother)."""
    prob_col = f"GT_{proband_id}" if not str(proband_id).startswith("GT_") else str(proband_id)
    pat_col = f"GT_{paternal_id}" if not str(paternal_id).startswith("GT_") else str(paternal_id)
    mat_col = f"GT_{maternal_id}" if not str(maternal_id).startswith("GT_") else str(maternal_id)

    if not all(col in data.columns for col in [prob_col, pat_col, mat_col]):
        logger.warning(f"Missing required genotype columns ({prob_col}, {pat_col}, {mat_col}) for compound het evaluation.")
        return data

    for symbol, group in data.groupby("SYMBOL"):
        if len(group) < 2:
            continue
        
        het_rows = group[group[prob_col].apply(_normalize_gt) == "HET"]
        if len(het_rows) < 2:
            continue

        paternal_only_indices = []
        maternal_only_indices = []

        for idx, row in het_rows.iterrows():
            pat_gt = _normalize_gt(row[pat_col])
            mat_gt = _normalize_gt(row[mat_col])

            if pat_gt in ["HET", "HOMALT"] and mat_gt == "HOMREF":
                paternal_only_indices.append(idx)
            elif pat_gt == "HOMREF" and mat_gt in ["HET", "HOMALT"]:
                maternal_only_indices.append(idx)

        if paternal_only_indices and maternal_only_indices:
            marked_indices = set(paternal_only_indices + maternal_only_indices)
            data.loc[list(marked_indices), "compound_het"] = True
            logger.info(f"Marked {len(marked_indices)} TRANS compound heterozygous variants in gene {symbol} for proband {proband_id}.")

            for p_idx in paternal_only_indices:
                p_loc = data.loc[p_idx, "Locus"] if "Locus" in data.columns else f"var_{p_idx}"
                mat_partners = []
                for m_idx in maternal_only_indices:
                    m_loc = data.loc[m_idx, "Locus"] if "Locus" in data.columns else f"var_{m_idx}"
                    m_hgvsc = str(data.loc[m_idx, "HGVSc"]) if "HGVSc" in data.columns and pd.notna(data.loc[m_idx, "HGVSc"]) and str(data.loc[m_idx, "HGVSc"]) != "nan" else ""
                    m_hgvsp = str(data.loc[m_idx, "HGVSp"]) if "HGVSp" in data.columns and pd.notna(data.loc[m_idx, "HGVSp"]) and str(data.loc[m_idx, "HGVSp"]) != "nan" else ""
                    anno = f" ({m_hgvsc} / {m_hgvsp})" if m_hgvsc or m_hgvsp else ""
                    mat_partners.append(f"{m_loc}{anno}")
                pair_str = f"Paternal | Partner (Maternal): {'; '.join(mat_partners)}"
                data.loc[p_idx, "COMPOUND_HET_PAIR"] = pair_str

            for m_idx in maternal_only_indices:
                m_loc = data.loc[m_idx, "Locus"] if "Locus" in data.columns else f"var_{m_idx}"
                pat_partners = []
                for p_idx in paternal_only_indices:
                    p_loc = data.loc[p_idx, "Locus"] if "Locus" in data.columns else f"var_{p_idx}"
                    p_hgvsc = str(data.loc[p_idx, "HGVSc"]) if "HGVSc" in data.columns and pd.notna(data.loc[p_idx, "HGVSc"]) and str(data.loc[p_idx, "HGVSc"]) != "nan" else ""
                    p_hgvsp = str(data.loc[p_idx, "HGVSp"]) if "HGVSp" in data.columns and pd.notna(data.loc[p_idx, "HGVSp"]) and str(data.loc[p_idx, "HGVSp"]) != "nan" else ""
                    anno = f" ({p_hgvsc} / {p_hgvsp})" if p_hgvsc or p_hgvsp else ""
                    pat_partners.append(f"{p_loc}{anno}")
                pair_str = f"Maternal | Partner (Paternal): {'; '.join(pat_partners)}"
                data.loc[m_idx, "COMPOUND_HET_PAIR"] = pair_str

    return data


def _mark_unphased_candidates(data: pd.DataFrame):
    gt_cols = [c for c in data.columns if c.startswith("GT_")]
    if not gt_cols:
        return data
    
    for s_col in gt_cols:
        s_name = s_col.replace("GT_", "")
        het_mask = data[s_col].astype(str).str.upper().isin(["HET", "0/1", "1/0", "0|1", "1|0"])
        dup_mask = data[het_mask].duplicated(subset=["SYMBOL"], keep=False)
        if dup_mask.any():
            cand_indices = data[het_mask & dup_mask].index
            data.loc[cand_indices, "compound_het_candidate"] = True
            for symbol, group in data.loc[cand_indices].groupby("SYMBOL"):
                for idx in group.index:
                    partners = [f"{row['Locus']}" if "Locus" in data.columns else f"var_{i}" for i, row in group.iterrows() if i != idx]
                    data.loc[idx, "COMPOUND_HET_PAIR"] = f"Unphased Candidate | Partners ({s_name}): {'; '.join(partners)}"
    return data


def _mark_with_pandas(data: pd.DataFrame, pedigree: pd.DataFrame):
    pheno_col = "Phenotype" if "Phenotype" in pedigree.columns else ("phenotype" if "phenotype" in pedigree.columns else None)
    ind_col = "IndividualID" if "IndividualID" in pedigree.columns else ("sample" if "sample" in pedigree.columns else pedigree.columns[1])
    pat_col = "PaternalID" if "PaternalID" in pedigree.columns else None
    mat_col = "MaternalID" if "MaternalID" in pedigree.columns else None

    if not (pheno_col and pat_col and mat_col):
        _mark_unphased_candidates(data)
        return data

    relevant_pedigree = pedigree[(pedigree[pheno_col].isin([1, 2, "1", "2", "affected"])) & (pedigree[pat_col] != 0) & (pedigree[mat_col] != 0)]
    
    if relevant_pedigree.empty:
        _mark_unphased_candidates(data)
        return data

    for _, row in relevant_pedigree.iterrows():
        prob_id = row[ind_col]
        pat_id = row[pat_col]
        mat_id = row[mat_col]
        data = _mark_pairwise_trans(data, prob_id, pat_id, mat_id)

    return data


def _mark_with_dict(data: pd.DataFrame, pedigree: dict):
    if "proband" in pedigree or "proband_id" in pedigree:
        prob_id = pedigree.get("proband") or pedigree.get("proband_id")
        pat_id = pedigree.get("father") or pedigree.get("paternal_id")
        mat_id = pedigree.get("mother") or pedigree.get("maternal_id")
        if prob_id and pat_id and mat_id:
            data = _mark_pairwise_trans(data, prob_id, pat_id, mat_id)
        else:
            _mark_unphased_candidates(data)
    else:
        for sample, info in pedigree.items():
            if isinstance(info, dict) and info.get("Phenotype") in [1, 2, "1", "2", "affected"]:
                pat_id = info.get("PaternalID") or info.get("father")
                mat_id = info.get("MaternalID") or info.get("mother")
                if pat_id and mat_id and pat_id != 0 and mat_id != 0:
                    data = _mark_pairwise_trans(data, sample, pat_id, mat_id)
                else:
                    _mark_unphased_candidates(data)

    return data


def mark_compound_het(data: pd.DataFrame, pedigree: dict | pd.DataFrame | None = None) -> pd.DataFrame:
    """Identifies compound heterozygous variants and populates COMPOUND_HET_STATUS and COMPOUND_HET_PAIR."""
    if "compound_het" not in data.columns:
        data["compound_het"] = False
    if "compound_het_candidate" not in data.columns:
        data["compound_het_candidate"] = False
    if "COMPOUND_HET_PAIR" not in data.columns:
        data["COMPOUND_HET_PAIR"] = ""

    gt_cols = [c for c in data.columns if c.startswith("GT_")]
    if not gt_cols:
        return data

    if isinstance(pedigree, pd.DataFrame):
        data = _mark_with_pandas(data, pedigree)
    elif isinstance(pedigree, dict):
        data = _mark_with_dict(data, pedigree)
    else:
        _mark_unphased_candidates(data)

    status_list = []
    for idx in range(len(data)):
        if data["compound_het"].iloc[idx]:
            status_list.append("Phased TRANS Pair")
        elif data["compound_het_candidate"].iloc[idx]:
            status_list.append("Unphased Candidate Pair")
        else:
            status_list.append("No Pair")
    data["COMPOUND_HET_STATUS"] = status_list

    return data


def analyze_pedigree_inheritance(
    data: pd.DataFrame,
    pedigree: dict | pd.DataFrame | None = None
) -> pd.DataFrame:
    """Evaluates 5-model single-allele genetic inheritance across multi-sample variant datasets."""
    gt_cols = [c for c in data.columns if c.startswith("GT_")]
    if len(gt_cols) <= 1 and pedigree is None:
        inh_list = []
        for _, r in data.iterrows():
            g = _normalize_gt(r[gt_cols[0]]) if len(gt_cols) == 1 else "MISSING"
            if g == "HET":
                inh_list.append("Autosomal Dominant / Heterozygous")
            elif g == "HOMALT":
                inh_list.append("Homozygous Variant")
            else:
                inh_list.append("Unclassified")
        data["INHERITANCE_MODEL"] = inh_list
        data["SAMPLE_GENOTYPES_SUMMARY"] = ""
        return data

    if "compound_het" not in data.columns:
        data["compound_het"] = False
    if pedigree is not None:
        try:
            data = mark_compound_het(data, pedigree)
        except Exception as e:
            logger.warning(f"Note on mark_compound_het execution: {e}")

    proband_id = None
    paternal_id = None
    maternal_id = None
    affected_ids = []
    unaffected_ids = []

    if isinstance(pedigree, pd.DataFrame):
        pheno_col = "Phenotype" if "Phenotype" in pedigree.columns else ("phenotype" if "phenotype" in pedigree.columns else None)
        ind_col = "IndividualID" if "IndividualID" in pedigree.columns else ("sample" if "sample" in pedigree.columns else pedigree.columns[1])
        pat_col = "PaternalID" if "PaternalID" in pedigree.columns else None
        mat_col = "MaternalID" if "MaternalID" in pedigree.columns else None

        if pheno_col:
            aff_df = pedigree[pedigree[pheno_col].isin([2, 1, "2", "1", "affected", "Affected"])]
            unaff_df = pedigree[pedigree[pheno_col].isin([1, 0, "1", "0", "unaffected", "Unaffected"])]
            affected_ids = [f"GT_{x}" if not str(x).startswith("GT_") else str(x) for x in aff_df[ind_col].dropna().tolist()]
            unaffected_ids = [f"GT_{x}" if not str(x).startswith("GT_") else str(x) for x in unaff_df[ind_col].dropna().tolist()]

        trio_rows = pedigree[(pedigree[pat_col] != 0) & (pedigree[mat_col] != 0)] if pat_col and mat_col else pd.DataFrame()
        if not trio_rows.empty:
            prob_raw = trio_rows[ind_col].values[0]
            pat_raw = trio_rows[pat_col].values[0]
            mat_raw = trio_rows[mat_col].values[0]
            proband_id = f"GT_{prob_raw}" if not str(prob_raw).startswith("GT_") else str(prob_raw)
            paternal_id = f"GT_{pat_raw}" if not str(pat_raw).startswith("GT_") else str(pat_raw)
            maternal_id = f"GT_{mat_raw}" if not str(mat_raw).startswith("GT_") else str(mat_raw)

    elif isinstance(pedigree, dict):
        proband_id = pedigree.get("proband") or pedigree.get("proband_id")
        paternal_id = pedigree.get("father") or pedigree.get("paternal_id")
        maternal_id = pedigree.get("mother") or pedigree.get("maternal_id")
        affected_ids = pedigree.get("affected", [])
        unaffected_ids = pedigree.get("unaffected", [])

        if proband_id and not proband_id.startswith("GT_"): proband_id = f"GT_{proband_id}"
        if paternal_id and not paternal_id.startswith("GT_"): paternal_id = f"GT_{paternal_id}"
        if maternal_id and not maternal_id.startswith("GT_"): maternal_id = f"GT_{maternal_id}"
        affected_ids = [f"GT_{x}" if not str(x).startswith("GT_") else str(x) for x in affected_ids]
        unaffected_ids = [f"GT_{x}" if not str(x).startswith("GT_") else str(x) for x in unaffected_ids]

    if not proband_id and gt_cols:
        proband_id = gt_cols[0]
    if not affected_ids and proband_id:
        affected_ids = [proband_id]

    inheritance_list = []
    summary_list = []
    chrom_col = "CHROM" if "CHROM" in data.columns else None

    for idx, row in data.iterrows():
        gt_summary_parts = []
        for col in gt_cols:
            s_name = col.replace("GT_", "")
            norm_g = _normalize_gt(row[col]) if col in data.columns else "MISSING"
            gt_summary_parts.append(f"{s_name}:{norm_g}")
        summary_str = " | ".join(gt_summary_parts)
        summary_list.append(summary_str)

        prob_gt = _normalize_gt(row[proband_id]) if proband_id and proband_id in data.columns else "MISSING"
        pat_gt = _normalize_gt(row[paternal_id]) if paternal_id and paternal_id in data.columns else "MISSING"
        mat_gt = _normalize_gt(row[maternal_id]) if maternal_id and maternal_id in data.columns else "MISSING"

        model = "Unclassified"
        has_trio = bool(paternal_id and maternal_id and paternal_id in data.columns and maternal_id in data.columns)

        if has_trio and prob_gt in ["HET", "HOMALT"] and pat_gt == "HOMREF" and mat_gt == "HOMREF":
            model = "De Novo"
        elif has_trio and prob_gt == "HOMALT" and pat_gt == "HET" and mat_gt == "HET":
            model = "Autosomal Recessive (Hom)"
        elif has_trio and chrom_col and str(row[chrom_col]).upper().replace("CHR", "") == "X" and prob_gt in ["HET", "HOMALT"] and mat_gt == "HET" and pat_gt == "HOMREF":
            model = "X-Linked"
        elif prob_gt == "HET" or (affected_ids and any(_normalize_gt(row[a]) == "HET" for a in affected_ids if a in data.columns)):
            model = "Autosomal Dominant / Heterozygous"
        elif prob_gt == "HOMALT" or (affected_ids and any(_normalize_gt(row[a]) == "HOMALT" for a in affected_ids if a in data.columns)):
            model = "Homozygous Variant"

        inheritance_list.append(model)

    data["INHERITANCE_MODEL"] = inheritance_list
    data["SAMPLE_GENOTYPES_SUMMARY"] = summary_list
    return data
