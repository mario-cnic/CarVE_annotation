#!/usr/bin/env python3
"""
Build GRCh38 tabix-indexed BED file of top LaBranchoR predicted branch points per 3' splice site.
"""
import os
import sys
import gzip
import subprocess
from pyliftover import LiftOver

CHAIN_FILE = "/home/mruizp/data_lab_PGP/resources/UCSC_chain/hg19ToHg38.over.chain.gz"
HG19_BED = "/home/mruizp/data_lab_PGP/resources/annotation/labranchor/lstm.gencode_v19.hg19.top.bed.gz"
OUT_DIR = "/home/mruizp/data_lab_PGP/resources/annotation/labranchor"
OUT_RAW_BED = os.path.join(OUT_DIR, "labranchor_grch38_top.raw.bed")
OUT_FINAL_BED_GZ = os.path.join(OUT_DIR, "labranchor_grch38_top.bed.gz")

def build_grch38_labranchor():
    print(f"Loading LiftOver chain from {CHAIN_FILE}...")
    lo = LiftOver(CHAIN_FILE)
    
    print(f"Reading and lifting over {HG19_BED} to GRCh38...")
    total = 0
    lifted = 0
    
    with gzip.open(HG19_BED, "rt") as f_in, open(OUT_RAW_BED, "w") as f_out:
        for line in f_in:
            total += 1
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            chrom, start, end, acceptor, score, strand = parts[0], int(parts[1]), int(parts[2]), int(parts[3]), float(parts[4]), parts[5]
            
            res_bp = lo.convert_coordinate(chrom, start)
            res_acc = lo.convert_coordinate(chrom, acceptor)
            
            if res_bp and res_acc:
                chr38, start38 = res_bp[0][0], res_bp[0][1]
                end38 = start38 + (end - start)
                acc38 = res_acc[0][1]
                
                # Write: chrom, start, end, acceptor_pos, score, strand
                f_out.write(f"{chr38}\t{start38}\t{end38}\t{acc38}\t{score:.6f}\t{strand}\n")
                lifted += 1
                
    print(f"LiftOver complete: {lifted}/{total} branch points mapped ({lifted/total*100:.1f}%).")
    
    bgzip_bin = "/home/mruizp/conda_envs/spliceai_env/bin/bgzip"
    tabix_bin = "/home/mruizp/conda_envs/spliceai_env/bin/tabix"
    
    print("Sorting, bgzipping, and tabix-indexing GRCh38 BED file...")
    cmd_sort = f"sort -k1,1 -k2,2n {OUT_RAW_BED} | {bgzip_bin} -c > {OUT_FINAL_BED_GZ}"
    subprocess.run(cmd_sort, shell=True, check=True)
    
    cmd_tabix = f"{tabix_bin} -p bed -f {OUT_FINAL_BED_GZ}"
    subprocess.run(cmd_tabix, shell=True, check=True)
    
    if os.path.exists(OUT_RAW_BED):
        os.remove(OUT_RAW_BED)
        
    print(f"Successfully generated {OUT_FINAL_BED_GZ} and {OUT_FINAL_BED_GZ}.tbi!")

if __name__ == "__main__":
    build_grch38_labranchor()
