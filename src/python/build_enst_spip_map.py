#!/usr/bin/env python3
"""Map Ensembl transcripts (GENCODE 45 = VEP 111) to the RefSeq transcripts SPiP actually scores, by exon structure,
and audit the curated gene_transcript_mapping.txt.

Why: SPiP is transcript-level and keyed by `NM_` accession, but it uses its own (older) RefSeq snapshot
(RefFiles/dataRefSeqhg38.RData). Many current MANE accessions (e.g. NM_001371727.1) are absent from it while an
older accession with the SAME exon structure is present (NM_000813). Accession matching alone therefore leaves the
row's own transcript "without an entry". Two transcripts are treated as the same when they are on the same
chromosome and strand and have an identical intron chain (all splice junctions equal; first/last exon ends may
differ = UTR extent). Single-exon transcripts need identical start and end.

Only pairs verified to be the SAME GENE are kept: the ENST's GENCODE gene must equal the SPiP transcript's symbol, or
resolve to it through HGNC (current/previous/alias symbol). Pairs across different genes (readthrough or bicistronic loci,
e.g. PYURF/PIGY, LRRC51/LRTOMT) and pairs whose SPiP symbol HGNC does not know are dropped and counted in the manifest.
The SPiP database is exported from its RData by this script (--spip-rdata/--rscript) so the input is reproducible.

Outputs
  --out-map    resources/enst_to_spip_nm.grch38.gencode45.tsv : enst, enst_version, gene_id, gene_name, spip_nm,
               spip_symbol, match_type (identical_exons | same_intron_chain | single_exon_same_bounds)
  --out-audit  per-row audit of the curated Gen,NM,ENST table (gene symbol vs GENCODE, ENST present, NM vs MANE,
               NM present in SPiP, structural SPiP matches)
Read-only on all inputs. Coordinates: GTF 1-based inclusive -> converted to 0-based half-open to compare with SPiP's BED12.
"""
import argparse, csv, gzip, hashlib, json, re, subprocess, sys
from collections import Counter, defaultdict
from datetime import datetime, timezone


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def attr(s, key):
    m = re.search(rf'{key} "([^"]+)"', s)
    return m.group(1) if m else None


def read_gtf(path, want_mane=False):
    tx = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] == "transcript":
                tid = attr(f[8], "transcript_id")
                tx[tid] = {"chrom": f[0], "strand": f[6], "gene_id": attr(f[8], "gene_id"), "gene_name": attr(f[8], "gene_name"),
                           "tags": re.findall(r'tag "([^"]+)"', f[8]), "refseq": (attr(f[8], "db_xref") or "").replace("RefSeq:", ""),
                           "exons": []}
            elif f[2] == "exon":
                tid = attr(f[8], "transcript_id")
                tx.setdefault(tid, {"chrom": f[0], "strand": f[6], "gene_id": attr(f[8], "gene_id"),
                                    "gene_name": attr(f[8], "gene_name"), "tags": [], "refseq": "", "exons": []})
                tx[tid]["exons"].append((int(f[3]) - 1, int(f[4])))
    for t in tx.values():
        t["exons"].sort()
    return tx


def key_of(chrom, strand, exons):
    if len(exons) == 1:
        return (chrom, strand, "single", exons[0])
    return (chrom, strand, "chain", tuple((exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--gtf", required=True); ap.add_argument("--mane-gtf", required=True)
    ap.add_argument("--spip-db-tsv", required=True, help="SPiP RefSeq database as TSV (written from --spip-rdata when given)")
    ap.add_argument("--spip-rdata", help="SPiP RefFiles/dataRefSeqhg38.RData; exported to --spip-db-tsv with --rscript")
    ap.add_argument("--rscript"); ap.add_argument("--hgnc", required=True); ap.add_argument("--curated", required=True)
    ap.add_argument("--out-map", required=True); ap.add_argument("--out-audit", required=True)
    a = ap.parse_args()

    if a.spip_rdata:
        if not a.rscript:
            sys.exit("--spip-rdata needs --rscript")
        expr = ('e<-new.env(); load("%s", envir=e); x<-get("dataRefSeq",e); colnames(x)<-c("chrom","start0","end","nm","score",'
                '"strand","thickStart","thickEnd","rgb","nExons","blockSizes","blockStarts","symbol"); '
                'write.table(x,"%s",sep="\t",quote=FALSE,row.names=FALSE)' % (a.spip_rdata, a.spip_db_tsv))
        subprocess.run([a.rscript, "-e", expr], check=True)
    hgnc_cur, hgnc_other = defaultdict(set), defaultdict(set)
    for r in csv.DictReader(open(a.hgnc), delimiter="\t"):
        if not r["ensembl_gene_id"]:
            continue
        hgnc_cur[r["symbol"]].add(r["ensembl_gene_id"])
        for k in ("prev_symbol", "alias_symbol"):
            for sym in r[k].split("|"):
                if sym:
                    hgnc_other[sym].add(r["ensembl_gene_id"])

    def gene_check(t, spip_symbol):
        if t["gene_name"] == spip_symbol:
            return "same_symbol"
        known = hgnc_cur.get(spip_symbol) or hgnc_other.get(spip_symbol)
        if not known:
            return "dropped_symbol_unknown_to_hgnc"
        return "same_gene_hgnc" if t["gene_id"].split(".")[0] in known else "dropped_different_gene"

    gen = read_gtf(a.gtf)
    mane = read_gtf(a.mane_gtf)
    mane_nm = {t.split(".")[0]: v["refseq"].split(".")[0] for t, v in mane.items() if v["refseq"]}
    mane_nm_full = {t.split(".")[0]: v["refseq"] for t, v in mane.items() if v["refseq"]}

    spip = {}
    by_key = defaultdict(list)
    with open(a.spip_db_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sizes = [int(x) for x in r["blockSizes"].strip(",").split(",")]
            starts = [int(x) for x in r["blockStarts"].strip(",").split(",")]
            s0 = int(r["start0"])
            exons = sorted((s0 + st, s0 + st + sz) for st, sz in zip(starts, sizes))
            spip[r["nm"]] = {"symbol": r["symbol"], "exons": exons, "chrom": r["chrom"], "strand": r["strand"]}
            by_key[key_of(r["chrom"], r["strand"], exons)].append(r["nm"])

    n_pairs = 0
    enst_to_spip = defaultdict(list)
    with open(a.out_map, "w") as out:
        out.write("enst\tenst_version\tgene_id\tgene_name\tspip_nm\tspip_symbol\tmatch_type\tgene_check\n")
        gc = Counter()
        for tid, t in sorted(gen.items()):
            if not t["exons"]:
                continue
            for nm in by_key.get(key_of(t["chrom"], t["strand"], t["exons"]), []):
                s = spip[nm]
                chk = gene_check(t, s["symbol"])
                gc[chk] += 1
                if chk.startswith("dropped"):
                    continue
                mt = ("single_exon_same_bounds" if len(t["exons"]) == 1 else
                      "identical_exons" if s["exons"] == t["exons"] else "same_intron_chain")
                out.write(f"{tid.split('.')[0]}\t{tid}\t{t['gene_id']}\t{t['gene_name']}\t{nm}\t{s['symbol']}\t{mt}\t{chk}\n")
                enst_to_spip[tid.split(".")[0]].append((nm, mt))
                n_pairs += 1

    # ---- audit of the curated table
    audit = []
    seen_gene, seen_enst = defaultdict(int), defaultdict(int)
    rows = list(csv.DictReader(open(a.curated)))
    for r in rows:
        seen_gene[r["Gen"].strip().upper()] += 1; seen_enst[r["ENST"].strip().split(".")[0]] += 1
    gen_unv = {k.split(".")[0]: v for k, v in gen.items()}
    for r in rows:
        gene, nm, enst = r["Gen"].strip(), r["NM"].strip(), r["ENST"].strip().split(".")[0]
        g = gen_unv.get(enst)
        mane_for_enst = mane_nm_full.get(enst)
        issues = []
        if g is None: issues.append("ENST_not_in_GENCODE45")
        elif (g["gene_name"] or "").upper() != gene.upper(): issues.append(f"gene_symbol_differs(GENCODE={g['gene_name']})")
        if mane_for_enst is None:
            issues.append("ENST_not_MANE")
        elif mane_for_enst.split(".")[0] != nm.split(".")[0]:
            issues.append(f"NM_differs_from_MANE({mane_for_enst})")
        elif mane_for_enst != nm:
            issues.append(f"NM_version_differs_from_MANE({mane_for_enst})")
        in_spip = nm.split(".")[0] in spip
        struct = enst_to_spip.get(enst, [])
        if not in_spip and not struct: issues.append("no_SPiP_transcript_for_this_ENST")
        if seen_gene[gene.upper()] > 1: issues.append("gene_listed_more_than_once")
        if seen_enst[enst] > 1: issues.append("ENST_listed_more_than_once")
        audit.append({"gene": gene, "NM": nm, "ENST": enst, "gene_in_GENCODE45": g["gene_name"] if g else "",
                      "MANE_NM_for_ENST": mane_for_enst or "", "NM_accession_in_SPiP": in_spip,
                      "SPiP_by_structure": ";".join(f"{n}({m})" for n, m in struct), "issues": ";".join(issues) or "ok"})
    with open(a.out_audit, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(audit[0]), delimiter="\t"); w.writeheader(); w.writerows(audit)

    manifest = {"built": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
                "gtf": a.gtf, "gtf_sha256": sha256(a.gtf), "mane_gtf": a.mane_gtf, "mane_gtf_sha256": sha256(a.mane_gtf),
                "spip_db_tsv": a.spip_db_tsv, "spip_db_tsv_sha256": sha256(a.spip_db_tsv), "curated": a.curated,
                "gencode_transcripts": len(gen), "spip_transcripts": len(spip), "pairs": n_pairs, "gene_check_counts": dict(gc),
                "hgnc": a.hgnc, "hgnc_sha256": sha256(a.hgnc), "spip_rdata": a.spip_rdata, "spip_rdata_sha256": sha256(a.spip_rdata) if a.spip_rdata else None,
                "enst_with_spip_match": len(enst_to_spip), "map_sha256": sha256(a.out_map), "argv": sys.argv}
    json.dump(manifest, open(a.out_map + ".manifest.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in manifest.items() if k not in ("argv",)}, indent=1))
    print("curated audit rows with issues other than version-only:", sum(1 for r in audit if any(not i.startswith("NM_version") and i != "ok" for i in r["issues"].split(";"))))


if __name__ == "__main__":
    main()
