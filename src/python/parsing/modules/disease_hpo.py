"""
Disease Gene Curation, HPO Phenotype Matching, and ACMG SF v3.2 Actionability.
"""
import os
import functools
import numpy as np
import pandas as pd
from . import logger

# ACMG SF v3.2 Actionable Cardiac Genes List (PMID: 37347242)
ACMG_SF_V3_2_CARDIAC_GENES = {
    "MYBPC3", "MYH7", "TNNT2", "TNNI3", "TPM1", "MYL2", "MYL3", "ACTC1",
    "LMNA", "FLNC", "PKP2", "DSP", "DSC2", "DSG2", "TMEM43",
    "KCNQ1", "KCNH2", "SCN5A", "RYR2",
    "FBN1", "TGFBR1", "TGFBR2", "SMAD3", "ACTA2", "MYH11", "COL3A1"
}

HPO_GENE_PANELS = {
    "HP:0001638": ["MYBPC3", "MYH7", "TNNT2", "TNNI3", "TPM1", "MYL2", "MYL3", "ACTC1", "PRKAG2", "GLA", "LAMP2", "TTN", "LMNA", "DSP", "FLNC", "PKP2", "DSG2", "DSC2", "BAG3", "RBM20"], # Cardiomyopathy
    "HP:0001639": ["MYBPC3", "MYH7", "TNNT2", "TNNI3", "TPM1", "MYL2", "MYL3", "ACTC1", "PRKAG2", "GLA", "LAMP2", "CSRP3", "TNNC1"], # Hypertrophic cardiomyopathy
    "HP:0001644": ["TTN", "LMNA", "MYH7", "TNNT2", "BAG3", "RBM20", "FLNC", "DSP", "VCL", "DES", "TNNC1", "DMD", "NKX2-5"], # Dilated cardiomyopathy
    "HP:0001678": ["LMNA", "SCN5A", "NKX2-5", "TRPM4", "PRKAG2", "DES", "TNNT2", "GJB2"], # Atrioventricular block
    "HP:0001657": ["PKP2", "DSP", "DSG2", "DSC2", "JUP", "TMEM43", "FLNC", "DES"], # Arrhythmogenic right ventricular cardiomyopathy
    "HP:0001627": ["FBN1", "TGFBR1", "TGFBR2", "SMAD3", "TGFB2", "ACTA2", "MYH11", "COL3A1"], # Abnormal aorta morphology / Aortic aneurysm
    "HP:0003552": ["DMD", "CAPN3", "DYSF", "SGCG", "SGCA", "SGCB", "SGCD", "FKRP", "ANO5"], # Muscle weakness / Muscular dystrophy
    "HP:0001663": ["KCNQ1", "KCNH2", "SCN5A", "CALM1", "CALM2", "CALM3", "TRDN", "KCNJ2", "CACNA1C"], # Long QT syndrome
    "HP:0001685": ["SCN5A"], # Brugada syndrome
    "HP:0001658": ["RYR2", "SCN5A", "KCNH2", "KCNQ1", "LMNA", "FLNC", "RBM20", "DSP", "PKP2", "CASQ2"], # Ventricular tachycardia
    "HP:0001695": ["LMNA", "FLNC", "DSP", "PLN", "RYR2", "SCN5A", "KCNQ1", "KCNH2", "FBN1", "MYBPC3", "MYH7"], # Sudden cardiac arrest / death
}


@functools.lru_cache(maxsize=1)
def load_cardiac_disease_gene_curation() -> dict[str, dict]:
    """Loads and merges curated disease gene evidence databases."""
    v4_path = "/home/mruizp/data_lab_PGP/shared/utils/data/Complete_gene_list_V4.csv"
    v2_path = "/home/mruizp/data_lab_PGP/shared/utils/data/full_cardiac_gene_list_V2.csv"
    
    curation = {}
    if os.path.exists(v4_path):
        df_v4 = pd.read_csv(v4_path)
        for _, r in df_v4.iterrows():
            g = str(r.get("GENE", "")).strip().upper()
            if not g or g == "NAN": continue
            src = str(r.get("SOURCE", "")).strip()
            ev = str(r.get("EVIDENCE", "")).strip()
            curation[g] = {
                "symbol": g,
                "name": str(r.get("NAME", "")).strip(),
                "source": src,
                "v4_evidence": ev,
                "v2_priority": ""
            }

    if os.path.exists(v2_path):
        df_v2 = pd.read_csv(v2_path)
        for _, r in df_v2.iterrows():
            g = str(r.get("Gen", "")).strip().upper()
            if not g or g == "NAN": continue
            prio = str(r.get("Prioridad", "")).strip().rstrip(",")
            if g in curation:
                curation[g]["v2_priority"] = prio
            else:
                curation[g] = {
                    "symbol": g,
                    "name": str(r.get("Nombre de proteína", "")).strip(),
                    "source": "Cardiomyopathies",
                    "v4_evidence": "Candidate",
                    "v2_priority": prio
                }

    DCM_CORE = {"BAG3", "DES", "DMD", "FLNC", "LMNA", "MYH7", "NKX2-5", "RBM20", "SCN5A", "TNNC1", "TNNT2", "TTN", "DSP", "ACTC1", "ACTN2", "JPH2", "NEXN", "PLN", "TPM1", "VCL", "MYBPC3", "PKP2", "DSG2", "DSC2", "PRDM16", "TBX20"}
    HCM_CORE = {"MYBPC3", "MYH7", "TNNT2", "TNNI3", "TPM1", "MYL2", "MYL3", "ACTC1", "CSRP3", "TNNC1", "JPH2", "PLN", "FHOD3", "GLA", "LAMP2", "PRKAG2", "TTR"}
    ARVC_CORE = {"PKP2", "DSP", "DSG2", "DSC2", "JUP", "TMEM43", "FLNC", "DES", "LMNA", "PLN", "CTNNA3", "RYR2"}
    CHANNELOPATHY_CORE = {"KCNQ1", "KCNH2", "SCN5A", "RYR2", "CASQ2", "CALM1", "CALM2", "CALM3", "TRDN", "KCNJ2", "CACNA1C", "TRPM4", "AKAP9", "ANK2", "KCNE1", "KCNE2"}
    AORTIC_CORE = {"FBN1", "TGFBR1", "TGFBR2", "SMAD3", "TGFB2", "ACTA2", "MYH11", "COL3A1", "MYLK", "LOX", "PRKG1", "EFEMP2", "ELN", "FLNA", "NOTCH1", "SLC2A10"}
    STORAGE_INFILTRATIVE_CORE = {"GLA", "LAMP2", "TTR", "PRKAG2", "GAA", "ALPK3", "AGGL"}

    for g, info in curation.items():
        v4_ev = info["v4_evidence"]
        v2_p = info["v2_priority"]
        
        if v4_ev == "Prioritary" or v2_p == "Prioritario" or g in DCM_CORE or g in HCM_CORE or g in ARVC_CORE or g in CHANNELOPATHY_CORE or g in AORTIC_CORE or g in STORAGE_INFILTRATIVE_CORE:
            cat = "Definitive / Primary Evidence"
        elif v2_p == "sarcomere_genes_GO":
            cat = "Strong / Sarcomeric Evidence"
        elif v4_ev == "Candidate" or v2_p == "Candidato":
            cat = "Moderate / Candidate Evidence"
        elif v2_p == "rna_heart_tissue":
            cat = "Heart Expressed (RNA-seq)"
        elif v2_p == "HPO_heart_related_genes" or v4_ev == "Secondary":
            cat = "HPO Cardiac Related"
        else:
            cat = "Moderate / Candidate Evidence"
        info["evidence_tier"] = cat

    return curation


def apply_disease_phenotype_curation(
    data: pd.DataFrame,
    disease_phenotype: str = "🫀 Cardiomyopathies (General Panel - DCM, HCM, ARVC, RCM)"
) -> pd.DataFrame:
    """Annotates dataset with Disease Curation Evidence levels."""
    curation_db = load_cardiac_disease_gene_curation()
    
    evidence_tiers = []
    curation_sources = []
    matches = []
    bonus_scores = []
    high_risk_flags = []
    actionable_flags = []

    phenotype_clean = disease_phenotype.lower()
    symbol_col = "SYMBOL" if "SYMBOL" in data.columns else ("Gene" if "Gene" in data.columns else None)

    DCM_CORE = {"BAG3", "DES", "DMD", "FLNC", "LMNA", "MYH7", "NKX2-5", "RBM20", "SCN5A", "TNNC1", "TNNT2", "TTN", "DSP", "ACTC1", "ACTN2", "JPH2", "NEXN", "PLN", "TPM1", "VCL", "MYBPC3", "PKP2", "DSG2", "DSC2", "PRDM16", "TBX20"}
    HCM_CORE = {"MYBPC3", "MYH7", "TNNT2", "TNNI3", "TPM1", "MYL2", "MYL3", "ACTC1", "CSRP3", "TNNC1", "JPH2", "PLN", "FHOD3", "GLA", "LAMP2", "PRKAG2", "TTR"}
    ARVC_CORE = {"PKP2", "DSP", "DSG2", "DSC2", "JUP", "TMEM43", "FLNC", "DES", "LMNA", "PLN", "CTNNA3", "RYR2"}
    CHANNELOPATHY_CORE = {"KCNQ1", "KCNH2", "SCN5A", "RYR2", "CASQ2", "CALM1", "CALM2", "CALM3", "TRDN", "KCNJ2", "CACNA1C", "TRPM4", "AKAP9", "ANK2", "KCNE1", "KCNE2"}
    AORTIC_CORE = {"FBN1", "TGFBR1", "TGFBR2", "SMAD3", "TGFB2", "ACTA2", "MYH11", "COL3A1", "MYLK", "LOX", "PRKG1", "EFEMP2", "ELN", "FLNA", "NOTCH1", "SLC2A10"}
    STORAGE_INFILTRATIVE_CORE = {"GLA", "LAMP2", "TTR", "PRKAG2", "GAA", "ALPK3", "AGGL"}

    if "dcm" in phenotype_clean:
        primary_panel = DCM_CORE
        high_risk_panel = {"LMNA", "FLNC", "RBM20", "DES", "DSP", "PLN"}
    elif "hcm" in phenotype_clean:
        primary_panel = HCM_CORE
        high_risk_panel = {"MYBPC3", "MYH7", "TNNT2", "PLN", "GLA", "LAMP2", "PRKAG2", "TTR"}
    elif "arvc" in phenotype_clean or "acm" in phenotype_clean:
        primary_panel = ARVC_CORE
        high_risk_panel = {"PKP2", "DSP", "DSG2", "DSC2", "TMEM43", "FLNC", "DES", "LMNA", "PLN"}
    elif "channelopathy" in phenotype_clean or "arrhythmia" in phenotype_clean:
        primary_panel = CHANNELOPATHY_CORE
        high_risk_panel = {"RYR2", "SCN5A", "KCNQ1", "KCNH2", "CALM1", "CALM2", "CALM3", "TRDN"}
    elif "aortic" in phenotype_clean or "vascular" in phenotype_clean:
        primary_panel = AORTIC_CORE
        high_risk_panel = {"FBN1", "TGFBR1", "TGFBR2", "SMAD3", "ACTA2", "COL3A1"}
    elif "infiltrative" in phenotype_clean or "storage" in phenotype_clean:
        primary_panel = STORAGE_INFILTRATIVE_CORE
        high_risk_panel = {"GLA", "LAMP2", "TTR", "PRKAG2"}
    else:
        primary_panel = DCM_CORE | HCM_CORE | ARVC_CORE | CHANNELOPATHY_CORE | AORTIC_CORE | STORAGE_INFILTRATIVE_CORE
        high_risk_panel = {"LMNA", "FLNC", "DSP", "PLN", "RYR2", "SCN5A", "KCNQ1", "KCNH2", "GLA", "TTR", "LAMP2", "FBN1"}

    for idx in range(len(data)):
        gene = str(data[symbol_col].iloc[idx]).strip().upper() if symbol_col else ""
        info = curation_db.get(gene, None)
        
        is_actionable = gene in ACMG_SF_V3_2_CARDIAC_GENES
        actionable_flags.append(is_actionable)
        
        is_high_risk = gene in high_risk_panel
        high_risk_flags.append(is_high_risk)

        if info:
            base_ev = info["evidence_tier"]
            source_val = info["source"] if info["source"] else "Cardiovascular Panel"

            is_match = (gene in primary_panel)
            if not is_match:
                if "cardiomyopathies" in phenotype_clean and ("cardiomyopathies" in source_val.lower()):
                    is_match = True
                elif "aortic" in phenotype_clean and ("aortic" in source_val.lower()):
                    is_match = True
                elif "conduction" in phenotype_clean and ("peripheral" in source_val.lower()):
                    is_match = True

            if is_match:
                ev_tier = "Definitive / Primary Evidence" if gene in primary_panel else base_ev
            else:
                ev_tier = "Moderate / Candidate Evidence" if base_ev == "Definitive / Primary Evidence" else base_ev

            bonus = 0.0
            if is_match:
                if ev_tier == "Definitive / Primary Evidence": bonus = 25.0
                elif ev_tier == "Strong / Sarcomeric Evidence": bonus = 18.0
                elif ev_tier == "Moderate / Candidate Evidence": bonus = 12.0
                elif ev_tier == "Heart Expressed (RNA-seq)": bonus = 6.0
                elif ev_tier == "HPO Cardiac Related": bonus = 3.0

            evidence_tiers.append(ev_tier)
            curation_sources.append(source_val)
            matches.append(is_match)
            bonus_scores.append(bonus)
        else:
            evidence_tiers.append("Uncurated / Other Genomic")
            curation_sources.append("Genome-Wide Background")
            matches.append(False)
            bonus_scores.append(0.0)

    data["DISEASE_GENE_EVIDENCE"] = pd.Categorical(
        evidence_tiers,
        categories=[
            "Definitive / Primary Evidence",
            "Strong / Sarcomeric Evidence",
            "Moderate / Candidate Evidence",
            "Heart Expressed (RNA-seq)",
            "HPO Cardiac Related",
            "Uncurated / Other Genomic"
        ],
        ordered=True
    )
    data["DISEASE_CURATION_SOURCE"] = curation_sources
    data["DISEASE_PHENOTYPE_MATCH"] = matches
    data["DISEASE_EVIDENCE_SCORE"] = bonus_scores
    data["HIGH_RISK_DISEASE_GENE"] = high_risk_flags
    data["ACMG_SF_V3_2"] = actionable_flags

    return data


def apply_hpo_phenotype_weights(data: pd.DataFrame, hpo_terms: list[str] | None = None) -> pd.DataFrame:
    """Matches patient HPO terms against curated gene panel databases."""
    if "HPO_MATCH" not in data.columns:
        data["HPO_MATCH"] = False

    if not hpo_terms or len(hpo_terms) == 0:
        return data

    target_genes = set()
    for term in hpo_terms:
        term_clean = term.split(" - ")[0].strip().upper()
        if not term_clean.startswith("HP:"):
            term_clean = f"HP:{term_clean}"
        if term_clean in HPO_GENE_PANELS:
            target_genes.update(HPO_GENE_PANELS[term_clean])

    if "SYMBOL" in data.columns and target_genes:
        hpo_mask = data["SYMBOL"].isin(target_genes)
        data["HPO_MATCH"] = hpo_mask
        if "VARIANT_PRIORITY_SCORE" in data.columns:
            data.loc[hpo_mask, "VARIANT_PRIORITY_SCORE"] = np.clip(
                data.loc[hpo_mask, "VARIANT_PRIORITY_SCORE"] + 15.0, 0.0, 100.0
            ).round(1)
            logger.info(f"Applied HPO bonus (+15.0 pts) to {hpo_mask.sum()} variants across {len(target_genes)} HPO-matched genes.")

    return data


def get_dynamic_priority(row, context: str):
    logger.debug(f"Evaluating gene priority for row with SYMBOL {row['SYMBOL']} and context {context}")
    if pd.isna(row.get('EVIDENCE')) or str(row.get('EVIDENCE')).strip() in ('', '.', 'nan', 'NaN', 'None'):
        return 'others'
    source = str(row['SOURCE'])
    evidence = str(row['EVIDENCE'])
    
    if context == 'Peripheral_and_conduction':
        if 'Peripheral_and_conduction' in source:
            return 'prioritary'
        return evidence.lower()

    elif context == 'Cardiomyopathies':
        if 'Cardiomyopathies' in source:
            return evidence.lower()
        if 'Peripheral_and_conduction' in source:
            return 'candidate'

    elif context == 'Aortic_and_vascular':
        if 'Aortic_and_vascular' in source:
            return evidence.lower()
        if 'Cardiomyopathies' in source:
            return 'candidate'
        return evidence.lower()
    
    elif context == 'all':
        if 'prioritary' in evidence.lower():
            return 'prioritary'
        return evidence.lower()

    logger.warning(f"Unrecognized context '{context}' for SYMBOL {row['SYMBOL']}; using evidence fallback")
    return evidence.lower()


def prioritize_genes(
    data: pd.DataFrame, gene_file: str, context: str = "all"
) -> pd.DataFrame:
    genes_df = pd.read_csv(gene_file)
    genes_df['GENE'] = genes_df['GENE'].str.strip()

    validate_gene_categories(
        context,
        genes_df['SOURCE'].str.split(";", expand=True).stack().unique()
    )
    
    data = data.merge(genes_df[['GENE', 'SOURCE', 'EVIDENCE']], 
                     left_on="SYMBOL", right_on="GENE", how="left", 
                     suffixes=('', '_genes'))
    
    suffix_cols = {'SOURCE_genes': 'SOURCE', 'EVIDENCE_genes': 'EVIDENCE'}
    for old_col, new_col in suffix_cols.items():
        if old_col in data.columns:
            data[new_col] = data[old_col]
            data = data.drop(columns=[old_col])
    
    data["gene_priority"] = data.apply(lambda row: get_dynamic_priority(row, context), axis=1)
    data["gene_priority"] = data["gene_priority"].fillna("others")
    
    gene_categories = ["prioritary", "secondary", "candidate", "others"]
    data["gene_priority"] = pd.Categorical(
        data["gene_priority"],
        categories=gene_categories,
        ordered=True
    )
    
    return data.drop(columns=["GENE"])


def validate_gene_categories(context, source_categories):
    if context != "all":
        categories = context.split(",")
        choices = source_categories.tolist()
        if any(x not in choices for x in categories):
            logger.error(
                f"Invalid gene category {categories}, should be one of {choices}"
            )
            raise ValueError(
                f"Invalid gene category {categories}, should be one of {choices}"
            )


def filter_by_gene_priority(data: pd.DataFrame, priorities: list[str]):
    if not priorities:
        logger.warning("No priorities provided for filtering, skipping")
        return data
    data = data[data["gene_priority"].isin(priorities)]
    logger.info(f"Dataframe filtered by gene priority {priorities}")
    return data
