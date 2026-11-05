#!/usr/bin/env python3
"""
Build the GRCh38 tabix-indexed TSV of LaBranchoR in silico mutagenesis (ISM) scores.

Source: the authors' genome-wide ISM file for GENCODE v19 (hg19) introns (jpaggi/labranchor,
notebooks/supp_data.ipynb; columns chrom, pos, ref, alt, strand, 3pSS, d_best_bp, d_max_bp). Every SNV in
the 70 bases upstream of a 3' splice site was scored by the LaBranchoR model on REF and ALT sequence:
    d_best_bp = P_alt(reference top branchpoint) - max P_ref   (change at the reference top branchpoint)
    d_max_bp  = max P_alt - max P_ref                          (change of the best branchpoint anywhere)
`pos` is 1-based; ref/alt are on the + strand; `3pSS` is the authors' 3' splice-site coordinate (pos lies
70..1 bases from it in transcript orientation).

Each row is lifted to GRCh38 only when it lifts cleanly (counts per reason go to the manifest):
  - `pos` and `3pSS` each map to exactly one position, on the + strand, on the same contig;
  - their distance is unchanged (no chain gap inside the window);
  - the GRCh38 FASTA base at the lifted `pos` equals `ref`.
pyliftover positions are 0-based: the 1-based `pos`/`3pSS` are converted on the way in and out.

Output (1-based `pos`, sorted, bgzipped, tabix -s1 -b2 -e2, header line starting with '#'):
    #chrom pos ref alt strand 3pSS d_best_bp d_max_bp
plus <output>.manifest.json (input / chain / FASTA sha256, counts).

Usage:
    build_labranchor_ism_grch38.py --ism-hg19 lstm.gencode_v19.hg19.ism.v2.tsv.gz \
        --chain hg19ToHg38.over.chain.gz --fasta Homo_sapiens_assembly38.fasta \
        --output labranchor_ism.gencode19.grch38.tsv.gz
"""
import argparse
import gzip
import hashlib
import json
import os
import subprocess
import tempfile
from collections import Counter
from datetime import datetime, timezone

HEADER = "#chrom\tpos\tref\talt\tstrand\t3pSS\td_best_bp\td_max_bp\n"


def lift_one(converter, chrom: str, pos1: int):
    """(contig, 1-based position) or a refusal reason. Unique, + strand mappings only."""
    hits = converter.convert_coordinate(chrom, pos1 - 1)
    if not hits:
        return None, "unmapped"
    if len(hits) > 1:
        return None, "multi_mapped"
    contig, pos0, strand = hits[0][0], hits[0][1], hits[0][2]
    if strand != "+":
        return None, "minus_strand_mapping"
    return (contig, pos0 + 1), None


def _fasta_contig(fasta, chrom: str):
    for name in (chrom, chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"):
        if name in fasta.references:
            return name
    return None


def lift_rows(lines, converter, fasta, counts: Counter):
    """Yields output lines for the rows of the hg19 ISM file that lift cleanly; counts every outcome."""
    pos_cache, pss_cache = {}, {}
    for line in lines:
        if line.startswith(("#", "chrom\t")):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8:
            counts["malformed"] += 1
            continue
        chrom, pos, ref, alt, strand, pss = f[0], int(f[1]), f[2].upper(), f[3].upper(), f[4], int(f[5])
        counts["rows"] += 1
        key = (chrom, pss)
        if key not in pss_cache:
            pss_cache[key] = lift_one(converter, chrom, pss)
        lifted_pss, why = pss_cache[key]
        if lifted_pss is None:
            counts[f"3pSS_{why}"] += 1
            continue
        key = (chrom, pos)
        if key not in pos_cache:
            lifted, why = lift_one(converter, chrom, pos)
            if lifted is not None:
                contig = _fasta_contig(fasta, lifted[0])
                base = fasta.fetch(contig, lifted[1] - 1, lifted[1]).upper() if contig else None
                lifted = (lifted, base)
            pos_cache[key] = (lifted, why)
        lifted, why = pos_cache[key]
        if lifted is None:
            counts[f"pos_{why}"] += 1
            continue
        (contig38, pos38), base38 = lifted
        if contig38 != lifted_pss[0]:
            counts["pos_and_3pSS_on_different_contigs"] += 1
        elif pos38 - lifted_pss[1] != pos - pss:
            counts["distance_to_3pSS_changed"] += 1
        elif base38 is None:
            counts["contig_not_in_fasta"] += 1
        elif base38 != ref:
            counts["ref_mismatch_grch38"] += 1
        else:
            counts["lifted"] += 1
            yield f"{contig38}\t{pos38}\t{ref}\t{alt}\t{strand}\t{lifted_pss[1]}\t{f[6]}\t{f[7]}\n"


def sort_and_index(raw: str, output: str, tmp_dir: str):
    """Sorts the lifted rows by contig and position, adds the header, bgzips and tabix-indexes."""
    import pysam
    plain = output[:-3]
    with open(plain, "w") as out:
        out.write(HEADER)
        out.flush()
        subprocess.run(["sort", "-k1,1", "-k2,2n", "-S", "2G", "-T", tmp_dir, raw],
                       stdout=out, check=True, env=dict(os.environ, LC_ALL="C"))
    pysam.tabix_compress(plain, output, force=True)
    os.remove(plain)
    pysam.tabix_index(output, seq_col=0, start_col=1, end_col=1, meta_char="#", zerobased=False, force=True)


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def build(ism_hg19: str, converter, fasta, output: str, tmp_dir: str | None = None) -> Counter:
    counts = Counter()
    with tempfile.TemporaryDirectory(dir=tmp_dir or os.path.dirname(os.path.abspath(output))) as tmp:
        raw = os.path.join(tmp, "lifted.tsv")
        opener = gzip.open if ism_hg19.endswith(".gz") else open
        with opener(ism_hg19, "rt") as fin, open(raw, "w") as fout:
            fout.writelines(lift_rows(fin, converter, fasta, counts))
        sort_and_index(raw, output, tmp)
    return counts


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ism-hg19", required=True, help="authors' hg19 ISM TSV (lstm.gencode_v19.hg19.ism.v2.tsv.gz)")
    p.add_argument("--chain", required=True, help="UCSC hg19ToHg38.over.chain.gz")
    p.add_argument("--fasta", required=True, help="GRCh38 FASTA the pipeline uses (params.fasta)")
    p.add_argument("--output", required=True, help="output .tsv.gz (a .tbi and a .manifest.json are written beside it)")
    p.add_argument("--tmp-dir", default=None, help="directory for the sort's temporary files")
    a = p.parse_args(argv)
    if not a.output.endswith(".gz"):
        p.error("--output must end in .gz")
    import pysam
    from pyliftover import LiftOver
    counts = build(a.ism_hg19, LiftOver(a.chain), pysam.FastaFile(a.fasta), a.output, a.tmp_dir)
    manifest = {"built": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "inputs": {k: {"path": os.path.abspath(v), "sha256": sha256(v)}
                           for k, v in (("ism_hg19", a.ism_hg19), ("chain", a.chain))},
                "fasta": {"path": os.path.abspath(a.fasta), "size_bytes": os.path.getsize(a.fasta)},
                "source": "LaBranchoR ISM, GENCODE v19 introns (jpaggi/labranchor, notebooks/supp_data.ipynb)",
                "counts": dict(sorted(counts.items())), "output_sha256": sha256(a.output)}
    with open(a.output + ".manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
