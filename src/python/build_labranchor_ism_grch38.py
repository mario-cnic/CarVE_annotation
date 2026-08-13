#!/usr/bin/env python3
"""
Build GRCh38 tabix-indexed TSV of LaBranchoR In Silico Mutagenesis (ISM) variant delta scores.
"""
import os
import sys
import gzip
import subprocess
from pyliftover import LiftOver

CHAIN_FILE = "/home/mruizp/data_lab_PGP/resources/UCSC_chain/hg19ToHg38.over.chain.gz"
HG19_ISM = "/home/mruizp/data_lab_PGP/resources/annotation/labranchor/lstm.gencode_v19.hg19.ism.v2.tsv.gz"
OUT_DIR = "/home/mruizp/data_lab_PGP/resources/annotation/labranchor"
OUT_RAW_TSV = os.path.join(OUT_DIR, "labranchor_grch38_ism.raw.tsv")
OUT_FINAL_TSV_GZ = os.path.join(OUT_DIR, "labranchor_grch38_ism.tsv.gz")

def build_grch38_ism():
    print(f"Loading LiftOver chain from {CHAIN_FILE}...")
    lo = LiftOver(CHAIN_FILE)
    
    print(f"Streaming and lifting over {HG19_ISM} to GRCh38...")
    total = 0
    lifted = 0
    
    with gzip.open(HG19_ISM, "rt") as f_in, open(OUT_RAW_TSV, "w") as f_out:
        # Read header
        header = f_in.readline()
        # Header for GRCh38: #chrom pos ref alt strand 3pSS d_best_bp d_max_bp
        f_out.write("#chrom\tpos\tref\talt\tstrand\t3pSS\td_best_bp\td_max_bp\n")
        
        for line in f_in:
            total += 1
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, ref, alt, strand, pss, d_best, d_max = parts[0], int(parts[1]), parts[2], parts[3], parts[4], int(parts[5]), float(parts[6]), float(parts[7])
            
            # Note: pos in hg19 is 1-based or 0-based? Let's check: in TSV pos is 1-based
            res_pos = lo.convert_coordinate(chrom, pos)
            res_pss = lo.convert_coordinate(chrom, pss)
            
            if res_pos and res_pss:
                chr38, pos38 = res_pos[0][0], res_pos[0][1]
                pss38 = res_pss[0][1]
                
                f_out.write(f"{chr38}\t{pos38}\t{ref}\t{alt}\t{strand}\t{pss38}\t{d_best:.6f}\t{d_max:.6f}\n")
                lifted += 1
                
            if total % 1000000 == 0:
                print(f"Processed {total:,} variants (lifted: {lifted:,})...")
                
    print(f"LiftOver complete: {lifted}/{total} variants mapped ({lifted/total*100:.1f}%).")
    
    bgzip_bin = "/home/mruizp/conda_envs/spliceai_env/bin/bgzip"
    tabix_bin = "/home/mruizp/conda_envs/spliceai_env/bin/tabix"
    
    print("Sorting, bgzipping, and tabix-indexing GRCh38 ISM table...")
    # Sort skipping header, then bgzip
    cmd_sort = f"(head -n 1 {OUT_RAW_TSV} && tail -n +2 {OUT_RAW_TSV} | sort -k1,1 -k2,2n) | {bgzip_bin} -c > {OUT_FINAL_TSV_GZ}"
    subprocess.run(cmd_sort, shell=True, check=True)
    
    # Index with tabix (1-based or 0-based: -s 1 -b 2 -e 2)
    cmd_tabix = f"{tabix_bin} -s 1 -b 2 -e 2 -f -S 1 {OUT_FINAL_TSV_GZ}"
    subprocess.run(cmd_tabix, shell=True, check=True)
    
    if os.path.exists(OUT_RAW_TSV):
        os.remove(OUT_RAW_TSV)
        
    print(f"Successfully generated {OUT_FINAL_TSV_GZ} and {OUT_FINAL_TSV_GZ}.tbi!")

if __name__ == "__main__":
    build_grch38_ism()
