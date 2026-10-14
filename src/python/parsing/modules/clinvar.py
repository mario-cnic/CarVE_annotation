"""Select ClinVar classifications and their review metadata from one source per row."""
import re
import numpy as np
import pandas as pd

SOURCES = (
    ("ClinVarVCF_CLNSIG", "ClinVarVCF_CLNREVSTAT"),
    ("clinvar_clnsig", "clinvar_review"),
    ("CLIN_SIG", None),
    ("CLNSIG", "CLNREVSTAT"),
    ("ClinicalSignificance", None),
    ("ClinVar", None),
    ("clinvar", None),
    ("CLINVAR_DISPLAY", None),
)


def _text(series):
    values = series.fillna("").astype(str).str.strip()
    missing = values.str.lower().isin(["", ".", "-", "nan", "none", "<na>", "na", "not_provided", "not provided"])
    return values.mask(missing, "")


def _review_stars(value):
    text = re.sub(r"\s+", " ", value.lower().replace("_", " ").replace("&", ",")).strip()
    text = re.sub(r",\s*", ", ", text)
    levels = {
        "practice guideline": 4,
        "reviewed by expert panel": 3,
        "criteria provided, multiple submitters, no conflicts": 2,
        "criteria provided, single submitter": 1,
        "criteria provided, conflicting classifications": 1,
        "criteria provided, conflicting interpretations": 1,
        "no assertion criteria provided": 0,
        "no assertion provided": 0,
        "no classification provided": 0,
    }
    return levels.get(text, float("nan"))


def resolve_clinvar(data):
    """Prefer custom VCF, dbNSFP, cache, then legacy aliases; fall back only on empty classification.

    VUS, conflicting and other nonempty classifications retain their source. Review
    metadata never falls back independently to another database release. With a `Locus`
    column, all rows of a variant take the source of its best-ranked row.
    """
    classification = pd.Series("", index=data.index, dtype=object)
    source = classification.copy()
    review = classification.copy()
    stars = pd.Series(float("nan"), index=data.index, dtype=float)
    for column, review_column in SOURCES:
        if column not in data.columns:
            continue
        values = _text(data[column])
        use = classification.eq("") & values.ne("")
        classification.loc[use] = values.loc[use]
        source.loc[use] = column
        if review_column and review_column in data.columns:
            reviews = _text(data[review_column])
            review.loc[use] = reviews.loc[use]
            stars.loc[use] = reviews.loc[use].map(_review_stars).astype(float)
        elif column not in ["ClinVarVCF_CLNSIG", "clinvar_clnsig", "CLIN_SIG"] and "ClinVar_Stars" in data.columns:
            numeric = pd.to_numeric(data["ClinVar_Stars"], errors="coerce")
            valid = numeric.isin([0, 1, 2, 3, 4])
            stars.loc[use & valid] = numeric.loc[use & valid]
    _share_best_source_per_variant(data, classification, source, review, stars)
    return classification, source, review, stars


def _share_best_source_per_variant(data, classification, source, review, stars):
    """Give every row of a variant the classification of its highest-priority source (in place).

    Rows without a `Locus` column keep their per-row result. Only variants with at least one
    classified row are looked at.
    """
    if "Locus" not in data.columns:
        return
    classified = classification.ne("").to_numpy()
    if not classified.any():
        return
    order = {column: i for i, (column, _) in enumerate(SOURCES)}
    rank = source.map(order).fillna(len(SOURCES)).to_numpy()
    key = data["Locus"]
    best = (pd.DataFrame({"key": key.to_numpy()[classified], "rank": rank[classified], "pos": np.flatnonzero(classified)})
            .sort_values("rank", kind="stable").drop_duplicates("key").set_index("key"))
    better = (key.map(best["rank"]).to_numpy() < rank)
    if not better.any():
        return
    donor = key.map(best["pos"]).to_numpy()[better].astype(int)
    for series in (classification, source, review, stars):
        series.iloc[np.flatnonzero(better)] = series.to_numpy()[donor]


def annotate_clinvar(data):
    resolved = resolve_clinvar(data)
    for column, values in zip(["CLINVAR_SELECTED_CLNSIG", "CLINVAR_SOURCE", "CLINVAR_SELECTED_REVIEW", "ClinVar_Stars"], resolved):
        data[column] = values
    return resolved


def clinvar_masks(classification):
    normalized = classification.str.lower().str.replace("_", " ", regex=False)
    tokens = normalized.str.split(r"[&/|;,]")
    pathogenic = tokens.map(lambda items: any(t.strip() in ["pathogenic", "likely pathogenic"] for t in items))
    benign = tokens.map(lambda items: any(t.strip() in ["benign", "likely benign"] for t in items))
    conflict = normalized.str.contains("conflicting|uncertain significance", regex=True) | (pathogenic & benign)
    return pathogenic & ~conflict, benign & ~conflict
