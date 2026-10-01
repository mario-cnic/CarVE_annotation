#!/usr/bin/env python3
"""Classify every dbNSFP column by what it depends on, from the readme AND from raw dbNSFP lines.

Why not the final table: our own parse_missense collapses MetaRNN_score / MetaLR_score / REVEL_score to the
maximum over transcripts, so the table hides that some of them are per-transcript lists in dbNSFP.

For a random sample of variants (loci from a parsed table) the raw lines are read with tabix; for each column
the number of ';'-separated items is compared with the number of items of Ensembl_transcriptid,
Ensembl_proteinid and Uniprot_acc. Classes:
  transcript_aligned  multi-valued rows have as many items as Ensembl_transcriptid in >= --min-fraction
                      (documented as corresponding to Ensembl_transcriptid / Ensembl_proteinid / Uniprot_*)
  unaligned_multi     multi-valued but not positionally aligned with the transcript list
  single              never multi-valued in the sample (variant-level or collapsed by dbNSFP)
Outputs a TSV (all columns) and the aligned-column list (one name per line). Read-only.
"""
import argparse, json, random, re, subprocess, sys
from collections import defaultdict
from datetime import datetime, timezone
import pandas as pd, pyarrow.parquet as pq

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--readme", required=True); ap.add_argument("--dbnsfp", required=True); ap.add_argument("--tabix", required=True)
ap.add_argument("--pq", required=True); ap.add_argument("--row-group", type=int, default=0)
ap.add_argument("--n-variants", type=int, default=400); ap.add_argument("--seed", type=int, default=7)
ap.add_argument("--min-multi-rows", type=int, default=30); ap.add_argument("--min-fraction", type=float, default=0.999)
ap.add_argument("--exclude", default="", help="comma list of columns never treated as aligned (reason is recorded)")
ap.add_argument("--out-tsv", required=True); ap.add_argument("--out-aligned", required=True)
a = ap.parse_args()

# ---- readme: field number, name, documented dependency
fields, cur = {}, None
for line in open(a.readme, errors="replace"):
    m = re.match(r"^(\d+)\t([A-Za-z0-9_.+()\-]+):", line)
    if m:
        cur = m.group(2); fields[cur] = {"n": int(m.group(1)), "text": line}
    elif cur and line.startswith("\t"):
        fields[cur]["text"] += line
    elif cur and line.strip() == "":
        pass
def documented(text):
    t = " ".join(text.split())
    m = re.search(r"(?:correspond(?:s|ing)? to|matching|match(?:es)?)\s+(Ensembl_transcriptid|Ensembl_proteinid|Uniprot_acc|Uniprot_entry)", t)
    if m: return m.group(1)
    if re.search(r"[Mm]ultiple [a-z ]*separated by", t): return "multi_unspecified"
    return "none"

# ---- raw sample
t = pq.ParquetFile(a.pq).read_row_group(a.row_group, columns=["Locus", "Ensembl_transcriptid"]).to_pandas()
loci = sorted(set(t[t.Ensembl_transcriptid.notna()].Locus))
random.Random(a.seed).shuffle(loci); loci = loci[:a.n_variants]
header = subprocess.run(f"zcat {a.dbnsfp} | head -1", shell=True, capture_output=True, text=True).stdout.rstrip("\n").split("\t")
header[0] = header[0].lstrip("#")
ix = {c: i for i, c in enumerate(header)}
stat = defaultdict(lambda: defaultdict(int))
n_lines = 0
for loc in loci:
    chrom, rest = loc.replace("chr", "").split(":"); pos, ref, alt = rest.split("-")
    out = subprocess.run([a.tabix, a.dbnsfp, f"{chrom}:{pos}-{pos}"], capture_output=True, text=True).stdout
    for ln in out.splitlines():
        f = ln.split("\t")
        if f[ix["ref"]] != ref or f[ix["alt"]] != alt: continue
        n_lines += 1
        n = {k: len(f[ix[k]].split(";")) for k in ("Ensembl_transcriptid", "Ensembl_proteinid", "Uniprot_acc")}
        for c, i in ix.items():
            v = f[i]
            if v in (".", ""): continue
            k = len(v.split(";"))
            stat[c]["nonnull"] += 1
            if k > 1:
                stat[c]["multi"] += 1
                stat[c]["eq_T"] += k == n["Ensembl_transcriptid"]
                stat[c]["eq_P"] += k == n["Ensembl_proteinid"]
                stat[c]["eq_U"] += k == n["Uniprot_acc"]
excluded = set(x for x in a.exclude.split(',') if x)
rows = []
for c in header:
    s = stat[c]; multi = s["multi"]
    fT = s["eq_T"] / multi if multi else float("nan")
    cls = "excluded_ambiguous_separator" if c in excluded else "single" if multi == 0 else ("transcript_aligned" if multi >= a.min_multi_rows and fT >= a.min_fraction else
                                       ("unaligned_multi" if multi >= a.min_multi_rows else "multi_too_few_rows"))
    rows.append({"column": c, "documented_dependency": documented(fields.get(c, {"text": ""})["text"]),
                 "nonnull_rows": s["nonnull"], "multi_rows": multi, "frac_len_eq_transcripts": round(fT, 4) if multi else "",
                 "frac_len_eq_proteins": round(s["eq_P"] / multi, 4) if multi else "",
                 "frac_len_eq_uniprot": round(s["eq_U"] / multi, 4) if multi else "", "class": cls})
df = pd.DataFrame(rows)
df.to_csv(a.out_tsv, sep="\t", index=False)
aligned = sorted(set(df[df["class"] == "transcript_aligned"].column) | {"Ensembl_transcriptid"})
open(a.out_aligned, "w").write("\n".join(aligned) + "\n")
json.dump({"derived": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"), "variants_sampled": len(loci),
           "raw_lines_used": n_lines, "seed": a.seed, "excluded": sorted(excluded), "readme": a.readme, "dbnsfp": a.dbnsfp, "argv": sys.argv},
          open(a.out_tsv + ".manifest.json", "w"), indent=1)
print(df["class"].value_counts().to_dict())
print(pd.crosstab(df["documented_dependency"], df["class"]).to_string())
