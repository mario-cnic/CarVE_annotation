"""Resolve explicit population frequencies without replacing unknown values with zero."""
import pandas as pd

GNOMAD_GRPMAX = "gnomADv4_AF_grpmax_joint"
GNOMAD_JOINT = "gnomADv4_AF_joint"
GNOMAD_AN_JOINT = "gnomADv4_AN_joint"
ABSENT = "absent_from_gnomad"
AF_COLUMNS = ("POPULATION_AF", "POPULATION_AF_SOURCE", "POPULATION_AF_STATUS")


def _missing(raw: pd.Series) -> pd.Series:
    return raw.isna() | raw.astype(str).str.strip().isin(["", ".", "-", "nan", "None", "<NA>", "NA"])


def resolve_population_af(data: pd.DataFrame, freq_col=GNOMAD_GRPMAX):
    """Return AF, source and availability; default grpmax falls back to joint per row.

    Finite frequencies must lie in [0, 1]. Status `absent_from_gnomad`: the gnomAD v4
    joint annotation is in the table but the allele has no record (empty AF and AN);
    AF stays null and site coverage is not checked. A record with AN 0 stays missing.
    Explicit custom columns have no fallback.
    """
    columns = [GNOMAD_GRPMAX, GNOMAD_JOINT] if freq_col == GNOMAD_GRPMAX else [freq_col]
    af = pd.Series(float("nan"), index=data.index, dtype=float)
    source = pd.Series("", index=data.index, dtype=object)
    invalid = pd.Series(False, index=data.index)
    for column in columns:
        if column not in data.columns:
            continue
        raw = data[column]
        missing = _missing(raw)
        numeric = pd.to_numeric(raw, errors="coerce")
        valid = numeric.between(0, 1) & ~missing
        invalid |= ~missing & ~valid
        use = af.isna() & valid
        af.loc[use] = numeric.loc[use]
        source.loc[use] = column
    status = pd.Series("missing_frequency", index=data.index, dtype=object)
    status.loc[invalid & af.isna()] = "invalid_frequency"
    status.loc[af.notna()] = "available_frequency"
    if freq_col in (GNOMAD_GRPMAX, GNOMAD_JOINT) and GNOMAD_JOINT in data.columns:
        no_record = _missing(data[GNOMAD_JOINT])
        if GNOMAD_AN_JOINT in data.columns:
            no_record &= _missing(data[GNOMAD_AN_JOINT])
        status.loc[no_record & ~invalid & af.isna()] = ABSENT
    return af, source, status


def rarity_af(af: pd.Series, status: pd.Series) -> pd.Series:
    """Frequency for rarity thresholds: an allele absent from gnomAD counts as zero observed."""
    return af.mask(status.eq(ABSENT), 0.0)


def annotate_population_af(data: pd.DataFrame, freq_col=GNOMAD_GRPMAX):
    """Attach the frequency used for interpretation and its source/availability."""
    resolved = resolve_population_af(data, freq_col)
    for column, values in zip(AF_COLUMNS, resolved):
        data[column] = values
    return resolved[0]
