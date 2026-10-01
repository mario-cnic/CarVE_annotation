"""Gene / transcript identity resolution for predictor entries.

Predictors write one variant-level string listing every gene (SpliceAI, Pangolin) or transcript
(SPiP) a variant touches; the table has one row per variant x VEP transcript. These helpers decide
which entry belongs to a row: by Ensembl gene ID (genes) or RefSeq NM accession (SPiP) -- never by
position in the string, and never by falling back to another gene's entry.

SpliceAI entries are keyed by an old gene *symbol*. Symbols are resolved to Ensembl gene IDs with two
independent routes that must agree:
  1. coordinates: exon overlap between SpliceAI's own annotation and GENCODE 45 (pre-built table,
     src/python/build_gene_identity_tables.py);
  2. HGNC: current / previous / alias symbol -> Ensembl gene ID (HGNC complete set).
Disagreement -> "ambiguous", no route -> "unresolved"; neither is scored.
"""
import csv
import re
from typing import Iterable, Optional

# Per-row match statuses (kept distinct from the predictors' own "not_covered").
MATCHED = "matched"
ALIAS_RESOLVED = "alias_resolved"
NO_ENTRY = "no_entry_for_row_gene"
NO_ROW_GENE = "no_row_gene"
NOT_APPLICABLE_TRANSCRIPT = "not_applicable_transcript"
AMBIGUOUS = "ambiguous"
UNRESOLVED = "unresolved"
NO_PREDICTION = "no_prediction"
PREDICTOR_ERROR = "predictor_error"
MATCHED_STRUCTURE = "matched_by_exon_structure"
MATCH_STATUSES = (MATCHED, MATCHED_STRUCTURE, ALIAS_RESOLVED, NO_ENTRY, NO_ROW_GENE, NOT_APPLICABLE_TRANSCRIPT, PREDICTOR_ERROR,
                  AMBIGUOUS, UNRESOLVED, NO_PREDICTION)

# Entries of a multi-gene value are joined with '&' in the table (',' in the raw VCF INFO).
ENTRY_SEP = re.compile(r"[&,]")
MISSING = (None, "", ".", "nan", "N/A")


def strip_version(identifier) -> str:
    return str(identifier).strip().split(".")[0]


def is_missing(value) -> bool:
    return value is None or (isinstance(value, float) and value != value) or str(value).strip() in MISSING


def alt_from_locus(locus) -> Optional[str]:
    """ALT allele from the table's `Locus` ("chr7:148871403-TA-T"); None if not parseable."""
    if is_missing(locus):
        return None
    parts = str(locus).split("-")
    return parts[-1] if len(parts) >= 3 else None


def split_entries(value) -> list:
    return [e for e in ENTRY_SEP.split(str(value)) if e]


class GeneIdentity:
    """Resolves predictor symbols / accessions to Ensembl gene IDs (versions stripped)."""

    def __init__(self, spliceai_symbol_map: str, hgnc_table: str, enst_spip_map: str):
        for label, path in (("spliceai_symbol_map", spliceai_symbol_map), ("hgnc_table", hgnc_table),
                            ("enst_spip_map", enst_spip_map)):
            if not path:
                raise ValueError(f"GeneIdentity: {label} is required (no silent fallback to unmatched scores)")
        self._coord = {}
        with open(spliceai_symbol_map) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                self._coord[row["symbol"]] = row["ensg"]
        self._hgnc_current, self._hgnc_other = {}, {}
        with open(hgnc_table) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                ens = row.get("ensembl_gene_id", "")
                if not ens:
                    continue
                self._hgnc_current.setdefault(row["symbol"], set()).add(ens)
                for col in ("prev_symbol", "alias_symbol"):
                    for sym in row.get(col, "").split("|"):
                        if sym:
                            self._hgnc_other.setdefault(sym, set()).add(ens)
        # ENST -> {SPiP NM} by identical exon structure (src/python/build_enst_spip_map.py). The curated
        # Gen,NM,ENST table is deliberately NOT used here: for 28 of its genes the NM and the ENST are different
        # transcripts, so pairing them would score the wrong isoform.
        self._spip_by_enst = {}
        with open(enst_spip_map) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                self._spip_by_enst.setdefault(row["enst"], set()).add(row["spip_nm"])
        self._cache = {}

    def resolve_symbol(self, symbol: str):
        """(ensg | None, status) for a SpliceAI entry symbol; status in {alias_resolved, ambiguous, unresolved}."""
        if symbol in self._cache:
            return self._cache[symbol]
        coord = self._coord.get(symbol)
        hgnc = self._hgnc_current.get(symbol) or self._hgnc_other.get(symbol) or set()
        if coord and hgnc:
            out = (coord, ALIAS_RESOLVED) if coord in hgnc else (None, AMBIGUOUS)
        elif coord:
            out = (coord, ALIAS_RESOLVED)
        elif len(hgnc) == 1:
            out = (next(iter(hgnc)), ALIAS_RESOLVED)
        elif hgnc:
            out = (None, AMBIGUOUS)
        else:
            out = (None, UNRESOLVED)
        self._cache[symbol] = out
        return out

    def hgnc_gene_ids(self, symbol: str) -> set:
        """Ensembl gene IDs HGNC lists for a current / previous / alias symbol (empty if unknown)."""
        return self._hgnc_current.get(symbol) or self._hgnc_other.get(symbol) or set()

    def symbol_is_row_gene(self, symbol, row_symbol, row_gene) -> bool:
        """True when a predictor's gene symbol denotes the row's gene (same symbol, or HGNC says so)."""
        if is_missing(symbol) or is_missing(row_gene):
            return False
        if not is_missing(row_symbol) and symbol == row_symbol:
            return True
        return strip_version(row_gene) in self.hgnc_gene_ids(symbol)

    def row_refseq_transcripts(self, feature, mane_select, mane_plus_clinical):
        """(accession_nms, structural_nms) for a row's Ensembl transcript, versions stripped.

        accession_nms: the RefSeq accession(s) VEP reports for the transcript (MANE Select / Plus Clinical).
        structural_nms: SPiP-database transcripts with the same chromosome, strand and intron chain.
        """
        acc = {strip_version(v) for v in (mane_select, mane_plus_clinical) if not is_missing(v)}
        struct = set(self._spip_by_enst.get(strip_version(feature), ())) if not is_missing(feature) else set()
        return acc, struct


def pick_gene_entry(entries: Iterable, row_gene, row_symbol, identity: Optional[GeneIdentity],
                    entry_gene_id=None, entry_symbol=None):
    """Select the entry belonging to the row's gene.

    `entries` are opaque objects; `entry_gene_id(e)` returns an Ensembl ID (Pangolin) or
    `entry_symbol(e)` returns a SpliceAI symbol. Returns (entry | None, status).
    """
    entries = list(entries)
    if is_missing(row_gene):
        return None, NO_ROW_GENE
    row_gene = strip_version(row_gene)
    worst = None
    for e in entries:
        if entry_gene_id is not None:
            if strip_version(entry_gene_id(e)) == row_gene:
                return e, MATCHED
            continue
        sym = entry_symbol(e)
        if not is_missing(row_symbol) and sym == row_symbol:
            return e, MATCHED
        ens, st = identity.resolve_symbol(sym)
        if ens == row_gene:
            return e, ALIAS_RESOLVED
        if st in (AMBIGUOUS, UNRESOLVED):
            worst = st if worst is None or st == AMBIGUOUS else worst
    return None, (worst or NO_ENTRY)
