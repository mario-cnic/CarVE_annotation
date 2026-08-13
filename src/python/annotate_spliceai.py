#!/usr/bin/env python3
import sys
import os
import argparse

# Suppress TensorFlow verbose logging
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

import warnings
warnings.filterwarnings('ignore')

try:
    import tensorflow as tf
    tf.get_logger().setLevel('ERROR')
    try:
        tf.config.threading.set_inter_op_parallelism_threads(2)
        tf.config.threading.set_intra_op_parallelism_threads(2)
    except Exception:
        pass
except Exception:
    pass

from spliceai.utils import Annotator, get_delta_scores
import pysam

def annotate_vcf(in_vcf: str, out_vcf: str, fasta_path: str, distance: int = 10000, mask: int = 0):
    """
    Annotates VCF variants with SpliceAI deep learning predictions at distance -d (default: 10000).
    """
    print(f"Initializing SpliceAI model weights and FASTA reference ({fasta_path})...")
    annotator = Annotator(fasta_path, "grch38")
    
    in_v = pysam.VariantFile(in_vcf)
    header = in_v.header
    if "SpliceAI" not in header.info:
        header.info.add("SpliceAI", ".", "String", "SpliceAI predicted splice scores (SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL)")
    
    out_v = pysam.VariantFile(out_vcf, "w", header=header)
    
    count = 0
    annotated = 0
    for rec in in_v:
        count += 1
        try:
            scores = get_delta_scores(rec, annotator, distance, mask)
            if scores:
                rec.info["SpliceAI"] = scores
                annotated += 1
        except Exception as e:
            pass
        out_v.write(rec)
        
    in_v.close()
    out_v.close()
    print(f"SpliceAI annotation finished. Total variants: {count}, Annotated with SpliceAI: {annotated}")

def main():
    parser = argparse.ArgumentParser(description="SpliceAI Local Annotation Engine at -D 10000")
    parser.add_argument("input_vcf", help="Input VCF file")
    parser.add_argument("output_vcf", help="Output VCF file")
    parser.add_argument("fasta", help="GRCh38 FASTA genome reference file")
    parser.add_argument("-d", "--distance", type=int, default=10000, help="Maximum intronic distance for splice predictions (Default: 10000)")
    args = parser.parse_args()
    
    annotate_vcf(args.input_vcf, args.output_vcf, args.fasta, distance=args.distance)

if __name__ == "__main__":
    main()
