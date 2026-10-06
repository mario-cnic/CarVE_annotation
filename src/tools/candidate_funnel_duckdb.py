#!/usr/bin/env python3
"""Candidate funnel over a parsed table, queried in place with DuckDB (no table is loaded into memory).

Counts rows and distinct Locus values that survive successive flags, per-sample genotype statuses and
the candidate loci carried per sample, and times each query. Thresholds are illustrative, not
validated for any disease. Only aggregate counts are printed.
"""
import argparse
import resource
import time

import duckdb  # keep before any pyarrow import (a prior pyarrow import breaks duckdb's libstdc++ here)

AF_COLUMN = "gnomADv4_AF_grpmax_joint"
AF_MAX = 0.001
SPLICE_MIN = 0.2
IMPACTS = ("HIGH", "MODERATE")
MAX_TIER = 2
CARRIER_STATUSES = ("het", "hom_alt", "alt_plus_other_allele")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("table", help="parsed.clean.pq")
    ap.add_argument("--genotypes", help="long genotype table (.pq)")
    ap.add_argument("--memory-limit", default="10GB")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--temp-dir", default="/tmp")
    a = ap.parse_args()

    con = duckdb.connect()
    con.execute(f"set memory_limit='{a.memory_limit}'; set threads={a.threads}; set temp_directory='{a.temp_dir}'")

    def timed(label, sql):
        t0 = time.time()
        rows = con.execute(sql).fetchall()
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
        print(f"[{time.time() - t0:5.1f}s, peak RSS {rss:.1f} GB] {label}")
        for r in rows:
            print("   ", r)

    impacts = ",".join(f"'{i}'" for i in IMPACTS)
    con.execute(f"""create or replace view f as select Locus,
        TRANSCRIPT_PRIORITY_TIER as tier,
        ({AF_COLUMN} is null or {AF_COLUMN} <= {AF_MAX}) as rare,
        (IMPACT in ({impacts})) as impact_hm,
        (coalesce(try_cast(spliceai_custom_MAX as double), 0) >= {SPLICE_MIN}
         or coalesce(try_cast(Pangolin_max_score as double), 0) >= {SPLICE_MIN}) as splice,
        (coalesce(CLIN_SIG, '') ilike '%pathogenic%') as clinvar_p
        from read_parquet('{a.table}')""")
    cand = "(rare and (impact_hm or splice)) or clinvar_p"
    timed("funnel (rows, distinct Locus)", f"""
        select 'all rows' step, count(*) n_rows, count(distinct Locus) n_loci from f
        union all select 'tier <= {MAX_TIER}', count(*), count(distinct Locus) from f where tier <= {MAX_TIER}
        union all select 'rare (AF <= {AF_MAX} or missing)', count(*), count(distinct Locus) from f where rare
        union all select 'rare and IMPACT HIGH or MODERATE', count(*), count(distinct Locus) from f where rare and impact_hm
        union all select 'rare and splice score >= {SPLICE_MIN}', count(*), count(distinct Locus) from f where rare and splice
        union all select 'ClinVar pathogenic term (CLIN_SIG), any AF', count(*), count(distinct Locus) from f where clinvar_p
        union all select 'candidate = rare and (impact or splice), or ClinVar pathogenic', count(*), count(distinct Locus) from f where {cand}
        union all select 'candidate, tier <= {MAX_TIER} rows only', count(*), count(distinct Locus) from f where ({cand}) and tier <= {MAX_TIER}""")
    if a.genotypes:
        con.execute(f"create or replace table cand as select distinct Locus from f where {cand}")
        carriers = ",".join(f"'{s}'" for s in CARRIER_STATUSES)
        timed("genotype statuses (sample order as stored, no sample names)", f"""
            select dense_rank() over (order by sample_id) sample_no, gt_status, count(*)
            from read_parquet('{a.genotypes}') group by sample_id, gt_status order by 1, 2""")
        timed("candidate loci carried, per sample", f"""
            select dense_rank() over (order by sample_id) sample_no, count(distinct Locus)
            from read_parquet('{a.genotypes}') where gt_status in ({carriers}) and Locus in (select Locus from cand)
            group by sample_id order by 1""")
        timed("candidate loci by number of carrier samples", f"""
            with x as (select Locus, count(distinct sample_id) n from read_parquet('{a.genotypes}')
                       where gt_status in ({carriers}) and Locus in (select Locus from cand) group by 1)
            select n carrier_samples, count(*) from x group by 1 order by 1""")


if __name__ == "__main__":
    main()
