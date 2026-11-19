"""Resolve explicit population frequencies without replacing unknown values with zero."""
import math

import pandas as pd

GNOMAD_GRPMAX = "gnomADv4_AF_grpmax_joint"
GNOMAD_JOINT = "gnomADv4_AF_joint"
GNOMAD_AN_JOINT = "gnomADv4_AN_joint"
ABSENT = "absent_from_gnomad"
ABSENT_LOW_AN = "absent_low_an"
AF_COLUMNS = ("POPULATION_AF", "POPULATION_AF_SOURCE", "POPULATION_AF_STATUS")
# gnomAD v4.1 allele number at every callable site (VEP --custom BED): one value per position
# the variant spans, joined by "&"; empty when gnomAD has no callable genotype there.
SITE_AN_SOURCES = ("gnomADv4AN_exomes", "gnomADv4AN_genomes")
SITE_AN = "GNOMAD_SITE_AN"
# PM2 frequency threshold that an absence must rule out: with 0 alleles observed among AN, the
# one-sided 95% Poisson upper bound of the AF is -ln(0.05)/AN, which must not exceed it.
PM2_ABSENCE_AF = 4e-5


def min_absent_an(pm2_af: float | None = None) -> int:
    """Smallest exomes+genomes AN at which absence bounds the AF below the PM2 threshold."""
    return math.ceil(-math.log(0.05) / (PM2_ABSENCE_AF if pm2_af is None else pm2_af))


def site_allele_number(data: pd.DataFrame) -> pd.Series:
    """Exomes + genomes AN at the variant's positions (minimum over positions; none = 0).

    NaN on every row when the table has no site-AN annotation.
    """
    present = [c for c in SITE_AN_SOURCES if c in data.columns]
    if not present:
        return pd.Series(float("nan"), index=data.index, dtype=float)
    total = pd.Series(0.0, index=data.index)
    for column in present:
        raw = data[column].astype("string").fillna("")
        parts = raw.str.split("&", expand=True).apply(pd.to_numeric, errors="coerce")
        total += parts.min(axis=1).fillna(0.0)
    return total


def _missing(raw: pd.Series) -> pd.Series:
    return raw.isna() | raw.astype(str).str.strip().isin(["", ".", "-", "nan", "None", "<NA>", "NA"])


def resolve_population_af(data: pd.DataFrame, freq_col=GNOMAD_GRPMAX):
    """Return AF, source and availability; default grpmax falls back to joint per row.

    Finite frequencies must lie in [0, 1]. Status `absent_from_gnomad`: the gnomAD v4
    joint annotation is in the table but the allele has no record (empty AF and AN);
    AF stays null. When the table has the site AN annotation, an absence whose site AN
    (exomes + genomes) is below `min_absent_an()` becomes `absent_low_an`: gnomAD did not
    genotype enough alleles there for the absence to mean rarity. A record with AN 0 stays
    missing. Explicit custom columns have no fallback.
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
        site_an = site_allele_number(data)
        status.loc[status.eq(ABSENT) & site_an.lt(min_absent_an())] = ABSENT_LOW_AN
    return af, source, status


def rarity_af(af: pd.Series, status: pd.Series) -> pd.Series:
    """Frequency for rarity thresholds: an allele absent from gnomAD counts as zero observed."""
    return af.mask(status.eq(ABSENT), 0.0)


def annotate_population_af(data: pd.DataFrame, freq_col=GNOMAD_GRPMAX):
    """Attach the frequency used for interpretation and its source/availability."""
    resolved = resolve_population_af(data, freq_col)
    for column, values in zip(AF_COLUMNS, resolved):
        data[column] = values
    if any(c in data.columns for c in SITE_AN_SOURCES):
        data[SITE_AN] = site_allele_number(data)
    return resolved[0]
