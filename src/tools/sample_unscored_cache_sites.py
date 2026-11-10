#!/usr/bin/env python3
"""Seeded random sample of the sites a score-cache version holds without a score, as a fixture VCF.

Reads every segment of <cache>/<predictor>/<version>/segments, keeps records whose first score field
is ".", one record per (CHROM, POS), and writes a sorted, bgzipped and indexed VCF with synthetic
trio genotypes (0/1, 0/0, 0/1). Header lines other than the contigs are copied from --header-from.
"""
import argparse
import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from annotation_cache import PREDICTORS  # noqa: E402

GT = "GT:GQ:DP\t0/1:99:30\t0/0:99:30\t0/1:99:30"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--version-dir", required=True, type=Path, help="<cache>/<predictor>/<version>")
    ap.add_argument("--predictor", default="spliceai", choices=sorted(PREDICTORS))
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--header-from", required=True, help="VCF whose non-contig header lines and #CHROM line are reused")
    ap.add_argument("--out", required=True, type=Path, help="output .vcf.gz")
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--tabix", default="tabix")
    ap.add_argument("--bgzip", default="bgzip")
    a = ap.parse_args()

    field = PREDICTORS[a.predictor]["fields"][0]
    segments = sorted((a.version_dir / "segments").glob("*.vcf.gz"))
    if not segments:
        sys.exit(f"no segments in {a.version_dir / 'segments'}")
    sites, contigs = {}, []
    for seg in segments:
        hdr = subprocess.run([a.bcftools, "view", "-h", str(seg)], check=True, capture_output=True, text=True).stdout
        contigs += [l for l in hdr.splitlines() if l.startswith("##contig") and l not in contigs]
        out = subprocess.run([a.bcftools, "query", "-i", f'INFO/{field}="."', "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(seg)],
                             check=True, capture_output=True, text=True).stdout
        for line in out.splitlines():
            c, p, r, alt = line.split("\t")
            sites.setdefault((c, int(p)), (r, alt))
    keys = sorted(sites)
    print(f"{len(segments)} segments; {len(keys)} unscored positions", file=sys.stderr)
    picked = sorted(random.Random(a.seed).sample(keys, min(a.n, len(keys))))

    hdr = subprocess.run([a.bcftools, "view", "-h", a.header_from], check=True, capture_output=True, text=True).stdout.splitlines()
    plain = a.out.with_suffix("") if a.out.suffix == ".gz" else a.out
    with open(plain, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.writelines(l + "\n" for l in contigs)
        fh.writelines(l + "\n" for l in hdr if l.startswith("##") and not l.startswith(("##contig", "##fileformat", "##bcftools")))
        fh.write(next(l for l in hdr if l.startswith("#CHROM")) + "\n")
        for c, p in picked:
            r, alt = sites[(c, p)]
            fh.write(f"{c}\t{p}\t.\t{r}\t{alt}\t100\tPASS\t.\t{GT}\n")
    subprocess.run([a.bcftools, "sort", "-Oz", "-o", str(a.out), str(plain)], check=True, capture_output=True)
    plain.unlink()
    subprocess.run([a.tabix, "-f", "-p", "vcf", str(a.out)], check=True)
    print(f"wrote {len(picked)} records to {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
