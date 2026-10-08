"""Inventory VEP fields, code references and Parquet column storage without reading rows."""
import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess

import pyarrow
import pyarrow.parquet as pq
import pysam


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def group(name):
    if name.startswith(("gnomADe_", "gnomADg_")):
        return "cache_gnomad"
    if name in {"MAX_AF", "MAX_AF_POPS"}:
        return "cache_max_af"
    if name in {"AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF"}:
        return "cache_1000genomes"
    if name.startswith("gnomADv4"):
        return "custom_gnomadv4"
    if name.startswith(("gnomAD_exomes_", "gnomAD_genomes_", "1000Gp3_", "ESP6500_", "ExAC_", "TWINSUK_", "ALSPAC_", "UK10K_", "ALFA_")):
        return "dbnsfp_legacy_population"
    if name.startswith(("clinvar_", "ClinVarVCF")) or name in {"CLIN_SIG", "SOMATIC", "PHENO", "PUBMED"}:
        return "clinical_existing_variation"
    if "rankscore" in name:
        return "dbnsfp_rankscore"
    return "other"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--dbnsfp", type=Path, required=True)
    parser.add_argument("--cache-info", type=Path, required=True)
    parser.add_argument("--consumer", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    args.outdir.mkdir(parents=True, exist_ok=False)
    with pysam.VariantFile(str(args.vcf)) as vcf:
        csq = vcf.header.info["CSQ"].description.split("Format: ", 1)[1].split("|")
    with gzip.open(args.dbnsfp, "rt") as fh:
        dbnsfp = set(fh.readline().lstrip("#").split())
    selector_counts = collections.Counter((root / "resources/all_but_old_gnomad_vep_cols.txt").read_text().split())
    selected = set(selector_counts)
    aligned = set((root / "resources/dbnsfp_transcript_aligned_columns.txt").read_text().splitlines())
    references = collections.defaultdict(list)
    code_manifest = []
    for label, repo in [("M2", root), ("M3", args.consumer.resolve())]:
        for rel in git(repo, "ls-files").splitlines():
            path = repo / rel
            if not rel.startswith(("src/", "modules/", "workflows/", "frontend/src/")):
                continue
            if "external/" in rel or path.suffix not in {".py", ".R", ".sh", ".nf", ".ts", ".tsx", ".js", ".jsx"}:
                continue
            content = path.read_text(errors="replace")
            code_manifest.append({"repo": label, "path": rel, "sha256": hashlib.sha256(content.encode()).hexdigest()})
            for number, line in enumerate(content.splitlines(), 1):
                tokens = set(re.findall(r"[A-Za-z_0-9][A-Za-z_0-9.]*", line))
                tokens.update(re.findall(r'''["']([^"'\s]+)["']''', line))
                for token in tokens:
                    references[token].append(f"{label}:{rel}:{number}")
    metadata = pq.ParquetFile(args.table).metadata
    sizes = collections.Counter()
    for r in range(metadata.num_row_groups):
        rg = metadata.row_group(r)
        for c in range(rg.num_columns):
            col = rg.column(c)
            sizes[col.path_in_schema] += col.total_compressed_size
    all_fields = list(dict.fromkeys(csq + list(sizes)))
    rows = []
    for name in all_fields:
        base = re.sub(r"\.\d+$", "", name).lstrip("_")
        refs = references[name] + (references[base] if base != name else [])
        rows.append({"field": name, "group": group(base), "in_csq": name in csq,
                     "csq_occurrences": csq.count(name), "in_dbnsfp_header": base in dbnsfp,
                     "in_parser_selector": base in selected or name in selected,
                     "transcript_alignment_list": base in aligned,
                     "in_table": name in sizes, "compressed_bytes": sizes[name],
                     "m2_references": " ; ".join(sorted(set(x for x in refs if x.startswith("M2:")))),
                     "m3_references": " ; ".join(sorted(set(x for x in refs if x.startswith("M3:"))))})
    with (args.outdir / "fields.tsv").open("w") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    candidates = [r for r in rows if r["in_csq"] and not r["m2_references"] and not r["m3_references"]]
    with (args.outdir / "unused_candidates.tsv").open("w") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(candidates)
    groups = {}
    for row in rows:
        item = groups.setdefault(row["group"], {"csq_fields": 0, "table_fields": 0, "compressed_bytes": 0, "no_code_token_reference": 0})
        item["csq_fields"] += row["in_csq"]
        item["table_fields"] += row["in_table"]
        item["compressed_bytes"] += row["compressed_bytes"]
        item["no_code_token_reference"] += row["in_csq"] and not row["m2_references"] and not row["m3_references"]
    command = (root / "modules/local/vep.nf").read_text()
    custom_fields = re.search(r"short_name=gnomADv4,.*?fields=([^\s]+)", command).group(1).split("%")
    summary = {"csq_slots": len(csq), "csq_unique_fields": len(set(csq)),
               "no_explicit_code_reference_fields": len(candidates),
               "selector_slots": sum(selector_counts.values()), "selector_unique_fields": len(selected),
               "selector_duplicates": {k: v for k, v in selector_counts.items() if v > 1},
               "csq_duplicates": {k: v for k, v in collections.Counter(csq).items() if v > 1},
               "dbnsfp_header_fields": len(dbnsfp), "dbnsfp_fields_in_csq": sum(n.lstrip("_") in dbnsfp for n in set(csq)),
               "table_rows": metadata.num_rows, "table_columns": metadata.num_columns,
               "table_row_groups": metadata.num_row_groups, "total_column_compressed_bytes": sum(sizes.values()),
               "groups": groups, "custom_gnomad_requested_slots": len(custom_fields),
               "custom_gnomad_unique_fields": len(set(custom_fields)),
               "custom_gnomad_duplicate_requests": {k: v for k, v in collections.Counter(custom_fields).items() if v > 1}}
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    manifest = {"timestamp": datetime.datetime.now().astimezone().isoformat(), "command": __import__("sys").argv,
                "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "selector_sha256": hashlib.sha256((root / "resources/all_but_old_gnomad_vep_cols.txt").read_bytes()).hexdigest(),
                "alignment_list_sha256": hashlib.sha256((root / "resources/dbnsfp_transcript_aligned_columns.txt").read_bytes()).hexdigest(),
                "git_commit": git(root, "rev-parse", "HEAD"), "git_status": git(root, "status", "--short"),
                "consumer_commit": git(args.consumer, "rev-parse", "HEAD"), "consumer_status": git(args.consumer, "status", "--short"),
                "versions": {"pyarrow": pyarrow.__version__, "pysam": pysam.__version__},
                "inputs": [{"path": str(p.resolve()), "size": p.stat().st_size, "mtime_ns": p.stat().st_mtime_ns}
                           for p in [args.vcf, args.table, args.dbnsfp, args.cache_info]],
                "cache_info": args.cache_info.read_text(), "code_files": code_manifest,
                "method": "VCF and dbNSFP headers; Parquet footer only; code token references include comments and plumbing, not proof of functional use. No patient rows read or exported."}
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
