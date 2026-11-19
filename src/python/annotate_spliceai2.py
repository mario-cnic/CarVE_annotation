#!/usr/bin/env python3
"""SpliceAI2 precomputed-score lookup (GRCh38, informational).

The SpliceAI2 data (Illumina, precomputed_scores_v2.0/GRCh38) holds every SNV and the population-observed
indels inside GENCODE v48 protein-coding genes, one tabix-indexed TSV per variant type and chromosome
({snv,indel}/chrN.tsv.gz). A variant inside genes on both strands has one row per strand; the row with the
higher summary score is kept. Columns written (exact CHROM/POS/REF/ALT match, biallelic records only):
- SpliceAI2_score: spliceai2_summary_score, the maximum of the four delta scores (0-1).
- SpliceAI2_DS_AG/AL/DG/DL: acceptor/donor gain/loss delta score of the strongest change (0-1).
- SpliceAI2_DP_AG/AL/DG/DL: the matching *_dist_0 (position of that splice site relative to the variant).
- SpliceAI2_strand: strand of the gene the kept row was scored on.
A record with no matching row gets no tags; the table step reports it as not_covered, never as 0.
"""

import argparse
import gzip
import os

import pysam

INFO = [
    ("SpliceAI2_score", "Float", "spliceai2_summary_score", "SpliceAI2 summary score (max of the four delta scores, 0-1)"),
    ("SpliceAI2_DS_AG", "Float", "acceptor_gain_delta_score_0", "SpliceAI2 acceptor gain delta score, strongest change"),
    ("SpliceAI2_DS_AL", "Float", "acceptor_loss_delta_score_0", "SpliceAI2 acceptor loss delta score, strongest change"),
    ("SpliceAI2_DS_DG", "Float", "donor_gain_delta_score_0", "SpliceAI2 donor gain delta score, strongest change"),
    ("SpliceAI2_DS_DL", "Float", "donor_loss_delta_score_0", "SpliceAI2 donor loss delta score, strongest change"),
    ("SpliceAI2_DP_AG", "Integer", "acceptor_gain_dist_0", "SpliceAI2 acceptor gain position relative to the variant"),
    ("SpliceAI2_DP_AL", "Integer", "acceptor_loss_dist_0", "SpliceAI2 acceptor loss position relative to the variant"),
    ("SpliceAI2_DP_DG", "Integer", "donor_gain_dist_0", "SpliceAI2 donor gain position relative to the variant"),
    ("SpliceAI2_DP_DL", "Integer", "donor_loss_dist_0", "SpliceAI2 donor loss position relative to the variant"),
    ("SpliceAI2_strand", "String", "strand", "Strand of the gene the SpliceAI2 row was scored on"),
]
SCORE_COLUMN = "spliceai2_summary_score"


class SpliceAI2Lookup:
    """Opens {snv,indel}/<chrom>.tsv.gz under `scores_dir` on first use; a missing file means not covered."""

    def __init__(self, scores_dir: str):
        if not os.path.isdir(os.path.join(scores_dir, "snv")):
            raise FileNotFoundError(f"SpliceAI2 scores directory without snv/: {scores_dir}")
        self.scores_dir = scores_dir
        self._open = {}
        self._columns = None

    def _tabix(self, kind: str, chrom: str):
        name = chrom if chrom.startswith("chr") else f"chr{chrom}"
        if name == "chrMT":
            name = "chrM"
        key = (kind, name)
        if key not in self._open:
            path = os.path.join(self.scores_dir, kind, f"{name}.tsv.gz")
            tbx = pysam.TabixFile(path) if os.path.exists(path) else None
            if tbx is not None and self._columns is None:
                with gzip.open(path, "rt") as fh:   # the header line is skipped by tabix (-S 1), not a '#' line
                    header = fh.readline().rstrip("\n").split("\t")
                self._columns = {c: i for i, c in enumerate(header)}
                missing = [src for _, _, src, _ in INFO if src not in self._columns]
                if missing:
                    raise ValueError(f"{path}: header lacks SpliceAI2 columns {missing}")
            self._open[key] = (tbx, name)
        return self._open[key]

    def lookup(self, chrom: str, pos: int, ref: str, alt: str):
        """Values for the INFO tags (in INFO order) of the best-scoring matching row, or None."""
        kind = "snv" if len(ref) == 1 and len(alt) == 1 else "indel"
        tbx, name = self._tabix(kind, chrom)
        if tbx is None or name not in tbx.contigs:
            return None
        ref, alt = ref.upper(), alt.upper()
        best = None
        for row in tbx.fetch(name, pos - 1, pos):
            f = row.split("\t")
            if int(f[1]) != pos or f[2].upper() != ref or f[3].upper() != alt:
                continue
            score = float(f[self._columns[SCORE_COLUMN]])
            if best is None or score > best[0]:
                best = (score, f)
        if best is None:
            return None
        values = []
        for _, typ, src, _ in INFO:
            raw = best[1][self._columns[src]]
            values.append(float(raw) if typ == "Float" else int(raw) if typ == "Integer" else raw)
        return values

    def close(self):
        for tbx, _ in self._open.values():
            if tbx is not None:
                tbx.close()


def annotate_spliceai2_vcf(in_vcf: str, out_vcf: str, scores_dir: str):
    lookup = SpliceAI2Lookup(scores_dir)
    in_v = pysam.VariantFile(in_vcf)
    header = in_v.header
    for tag, typ, _, desc in INFO:
        if tag not in header.info:
            header.info.add(tag, "1", typ, desc)
    out_v = pysam.VariantFile(out_vcf, "w", header=header)
    total = scored = 0
    for rec in in_v:
        total += 1
        if rec.alts and len(rec.alts) == 1:
            values = lookup.lookup(rec.chrom, rec.pos, rec.ref, rec.alts[0])
            if values is not None:
                for (tag, _, _, _), value in zip(INFO, values):
                    rec.info[tag] = value
                scored += 1
        out_v.write(rec)
    in_v.close()
    out_v.close()
    lookup.close()
    print(f"SpliceAI2 lookup completed. Total: {total}, scored: {scored}")


def main():
    parser = argparse.ArgumentParser(description="SpliceAI2 precomputed-score lookup (GRCh38, informational)")
    parser.add_argument("input_vcf")
    parser.add_argument("output_vcf")
    parser.add_argument("--scores-dir", required=True,
                        help="precomputed_scores_v2.0/GRCh38 folder holding snv/ and indel/ chrN.tsv.gz with tabix indexes")
    args = parser.parse_args()
    annotate_spliceai2_vcf(args.input_vcf, args.output_vcf, args.scores_dir)


if __name__ == "__main__":
    main()
