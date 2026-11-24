#!/usr/bin/env python3
"""Per-transcript exon structure in cDNA coordinates, for the NMD rules in parse_nmd.

Input: a GENCODE GTF (GENCODE 45 = the gene models of the VEP 111 cache). Only transcripts with CDS features are kept.

Output TSV, one row per transcript:
  transcript_id        ENST without version (VEP's Feature)
  transcript_version   GENCODE version (compared with the version in HGVSc)
  chrom, strand
  exon_cdna_ends       cumulative exon ends in transcript order, 1-based cDNA (last value = transcript length)
  cds_start_cdna       cDNA position of the first CDS nucleotide
  cds_end_cdna         cDNA position of the last CDS nucleotide, stop codon included when annotated
  start_incomplete     1 when tagged cds_start_NF (5' end of the CDS unknown)
  end_incomplete       1 when tagged cds_end_NF or mRNA_end_NF (3' exon structure unknown)
A manifest (<out>.manifest.json) records the GTF path and sha256, counts and the output sha256.
Read-only on the input.
"""
import argparse, gzip, hashlib, json, re, sys
from collections import defaultdict
from datetime import datetime, timezone


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_gtf(path):
    """{transcript_id: record} with exons, CDS and stop-codon intervals (1-based inclusive)."""
    tx = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("transcript", "exon", "CDS", "stop_codon"):
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if not m:
                continue
            tid = m.group(1)
            r = tx.setdefault(tid, {"chrom": f[0], "strand": f[6], "exon": [], "CDS": [], "stop_codon": [], "tags": set()})
            if f[2] == "transcript":
                r["tags"].update(re.findall(r'tag "([^"]+)"', f[8]))
            else:
                r[f[2]].append((int(f[3]), int(f[4])))
    return tx


def to_cdna(pos, exons, strand):
    """1-based cDNA position of genomic `pos`; exons in transcript order."""
    done = 0
    for s, e in exons:
        if s <= pos <= e:
            return done + (pos - s + 1 if strand == "+" else e - pos + 1)
        done += e - s + 1
    raise ValueError(f"position {pos} not in exons")


def build_row(tid, r):
    if not r["CDS"] or not r["exon"]:
        return None
    rev = r["strand"] == "-"
    exons = sorted(r["exon"], reverse=rev)
    ends, done = [], 0
    for s, e in exons:
        done += e - s + 1
        ends.append(done)
    coding = r["CDS"] + r["stop_codon"]
    lo, hi = min(s for s, _ in coding), max(e for _, e in coding)
    first, last = (hi, lo) if rev else (lo, hi)
    base, _, version = tid.partition(".")
    return {
        "transcript_id": base,
        "transcript_version": version,
        "chrom": r["chrom"],
        "strand": r["strand"],
        "exon_cdna_ends": ",".join(map(str, ends)),
        "cds_start_cdna": to_cdna(first, exons, r["strand"]),
        "cds_end_cdna": to_cdna(last, exons, r["strand"]),
        "start_incomplete": int("cds_start_NF" in r["tags"]),
        "end_incomplete": int(bool({"cds_end_NF", "mRNA_end_NF"} & r["tags"])),
    }


COLUMNS = ["transcript_id", "transcript_version", "chrom", "strand", "exon_cdna_ends",
           "cds_start_cdna", "cds_end_cdna", "start_incomplete", "end_incomplete"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True, help="GENCODE GTF (gzipped or plain)")
    ap.add_argument("--out", required=True, help="output TSV")
    args = ap.parse_args(argv)

    tx = read_gtf(args.gtf)
    rows = [row for tid, r in tx.items() if (row := build_row(tid, r))]
    rows.sort(key=lambda x: x["transcript_id"])
    ids = [x["transcript_id"] for x in rows]
    if len(ids) != len(set(ids)):
        sys.exit("duplicate unversioned transcript IDs (PAR_Y copies?); refusing to write an ambiguous table")
    with open(args.out, "w") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for x in rows:
            out.write("\t".join(str(x[c]) for c in COLUMNS) + "\n")

    manifest = {
        "built": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "gtf": args.gtf,
        "gtf_sha256": sha256(args.gtf),
        "gtf_transcripts": len(tx),
        "coding_transcripts": len(rows),
        "start_incomplete": sum(x["start_incomplete"] for x in rows),
        "end_incomplete": sum(x["end_incomplete"] for x in rows),
        "single_exon": sum("," not in x["exon_cdna_ends"] for x in rows),
        "out_sha256": sha256(args.out),
        "argv": sys.argv,
    }
    with open(args.out + ".manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=1)
    print(json.dumps({k: v for k, v in manifest.items() if k != "argv"}, indent=1))


if __name__ == "__main__":
    main()
