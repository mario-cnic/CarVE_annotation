#!/usr/bin/env python3
"""Writes the genotypes of an annotated VCF as a long table: one row per variant x sample.

The VCF is streamed through `bcftools query`; FORMAT fields not declared in the header are skipped.
Every variant x sample pair is written, including homozygous-reference and no-call genotypes, each
with a `gt_status`. A missing FORMAT value is null, never 0. A pair absent from the table means
"no evidence", not homozygous reference.

Columns: Locus (CHROM:POS-REF-ALT, as in the main table), chrom, pos, ref, alt, sample_id, gt_raw,
gt_status, allele_1, allele_2, phased, ploidy, gq, dp, ad, ad_ref, ad_alt and, when declared in the
header, pl, ps, ft.

gt_status: hom_ref, het, hom_alt, no_call, partial_no_call (e.g. ./1), haploid_ref, haploid_alt,
multiallelic_other (e.g. 1/2 on an unsplit record), alt_plus_other_allele (unphased 1/0: see below).

An unphased "1/0" is not a heterozygote with a reference copy. Callers write heterozygotes as 0/1;
`bcftools norm -m -any` writes 1/0 when the sample carries this ALT plus a different ALT. In the S223
chr22 slice every 1/0 genotype had more reads than AD ref + alt (DP > AD sum) and 81% had no
reference reads, while no 0/1 genotype did.

An input with no sample columns, or without a GT FORMAT field, writes no file and exits 0.
"""

import argparse
import re
import subprocess
import sys

import pyarrow as pa
import pyarrow.parquet as pq

FORMAT_FIELDS = ("GT", "GQ", "DP", "AD", "PL", "PS", "FT")
GT_SPLIT = re.compile(r"[/|]")


def parse_gt(raw):
    """Returns (status, allele_1, allele_2, phased, ploidy) for a GT string such as '0/1' or '1|0'."""
    if raw in (None, "", "."):
        return "no_call", None, None, False, 1
    alleles = [None if a == "." else int(a) for a in GT_SPLIT.split(raw)]
    phased = "|" in raw
    ploidy = len(alleles)
    first = alleles[0]
    second = alleles[1] if ploidy > 1 else None
    called = [a for a in alleles if a is not None]
    if not called:
        status = "no_call"
    elif len(called) < ploidy:
        status = "partial_no_call"
    elif ploidy == 1:
        status = "haploid_ref" if first == 0 else "haploid_alt"
    elif all(a == 0 for a in called):
        status = "hom_ref"
    elif ploidy == 2 and not phased and first != 0 and second == 0:
        status = "alt_plus_other_allele"
    else:
        nonzero = {a for a in called if a != 0}
        if len(nonzero) > 1:
            status = "multiallelic_other"
        else:
            status = "hom_alt" if 0 not in called else "het"
    return status, first, second, phased, ploidy


def to_int(value):
    return None if value in (".", "") else int(value)


def to_str(value):
    return None if value in (".", "") else value


def split_ad(value):
    """Returns (ad_raw, ad_ref, ad_alt); ad_ref/ad_alt only for a two-value AD."""
    if value in (".", ""):
        return None, None, None
    parts = value.split(",")
    if all(x == "." for x in parts):
        return None, None, None
    if len(parts) == 2:
        return value, to_int(parts[0]), to_int(parts[1])
    return value, None, None


def header_info(bcftools, vcf):
    """Returns (sample names, declared FORMAT field ids)."""
    samples = subprocess.run([bcftools, "query", "-l", vcf], capture_output=True, text=True, check=True).stdout.split()
    header = subprocess.run([bcftools, "view", "-h", vcf], capture_output=True, text=True, check=True).stdout
    declared = set(re.findall(r"^##FORMAT=<ID=([^,>]+)", header, flags=re.M))
    return samples, declared


def batch_to_table(rows):
    types = {
        "pos": pa.int64(), "allele_1": pa.int32(), "allele_2": pa.int32(), "phased": pa.bool_(),
        "ploidy": pa.int16(), "gq": pa.int32(), "dp": pa.int32(), "ad_ref": pa.int32(), "ad_alt": pa.int32(),
    }
    return pa.table({col: pa.array([r[col] for r in rows], type=types.get(col, pa.string())) for col in rows[0]})


def run(vcf, output, bcftools, batch_rows):
    samples, declared = header_info(bcftools, vcf)
    if not samples:
        print("Genotype mode no_samples: no genotype table written")
        return 0
    if "GT" not in declared:
        print("Genotype mode samples_without_GT: no genotype table written")
        return 0
    fields = [f for f in FORMAT_FIELDS if f in declared]
    fmt = "%CHROM\\t%POS\\t%REF\\t%ALT[" + "".join("\\t%" + f for f in fields) + "]\\n"
    proc = subprocess.Popen([bcftools, "query", "-f", fmt, vcf], stdout=subprocess.PIPE, text=True)
    n_fields = len(fields)
    rows, writer, total = [], None, 0
    try:
        for line in proc.stdout:
            cells = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt = cells[0], int(cells[1]), cells[2], cells[3]
            locus = f"{chrom}:{pos}-{ref}-{alt}"
            for i, sample in enumerate(samples):
                vals = dict(zip(fields, cells[4 + i * n_fields: 4 + (i + 1) * n_fields]))
                gt_raw = vals["GT"]
                status, a1, a2, phased, ploidy = parse_gt(gt_raw)
                ad, ad_ref, ad_alt = split_ad(vals.get("AD", "."))
                row = {
                    "Locus": locus, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "sample_id": sample,
                    "gt_raw": gt_raw, "gt_status": status, "allele_1": a1, "allele_2": a2, "phased": phased,
                    "ploidy": ploidy, "gq": to_int(vals.get("GQ", ".")), "dp": to_int(vals.get("DP", ".")),
                    "ad": ad, "ad_ref": ad_ref, "ad_alt": ad_alt,
                }
                for f in ("PL", "PS", "FT"):
                    if f in vals:
                        row[f.lower()] = to_str(vals[f])
                rows.append(row)
            if len(rows) >= batch_rows:
                writer = write_batch(writer, rows, output, fields)
                total += len(rows)
                rows = []
        if rows:
            writer = write_batch(writer, rows, output, fields)
            total += len(rows)
    finally:
        if writer is not None:
            writer.close()
    if proc.wait() != 0:
        print("ERROR: bcftools query failed", file=sys.stderr)
        return 1
    if writer is None:
        # Samples declared but no records: write an empty table with the schema.
        pq.write_table(batch_to_table([empty_row(fields)]).slice(0, 0), output)
    print(f"Genotype table: {total} rows ({len(samples)} sample(s)) -> {output}")
    return 0


def empty_row(fields):
    row = {c: None for c in ("Locus", "chrom", "pos", "ref", "alt", "sample_id", "gt_raw", "gt_status", "allele_1",
                             "allele_2", "phased", "ploidy", "gq", "dp", "ad", "ad_ref", "ad_alt")}
    row.update({f.lower(): None for f in ("PL", "PS", "FT") if f in fields})
    return row


def write_batch(writer, rows, output, fields):
    table = batch_to_table(rows)
    if writer is None:
        writer = pq.ParquetWriter(output, table.schema)
    writer.write_table(table.cast(writer.schema))
    return writer


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="Annotated VCF (bgzipped)")
    p.add_argument("--output", required=True, help="Output Parquet path (not written when the input has no genotypes)")
    p.add_argument("--bcftools", default="bcftools", help="bcftools executable (default: bcftools on PATH)")
    p.add_argument("--batch-rows", type=int, default=500000, help="Rows per Parquet row group (default 500000)")
    a = p.parse_args()
    sys.exit(run(a.vcf, a.output, a.bcftools, a.batch_rows))


if __name__ == "__main__":
    main()
