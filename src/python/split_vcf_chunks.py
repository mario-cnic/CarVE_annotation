#!/usr/bin/env python3
"""
Header-Preserving Streaming VCF Chunker
Splits large VCF files into smaller, indexed VCF chunks for parallel compute.
"""

import os
import sys
import argparse
import pysam
import subprocess


def split_vcf_into_chunks(input_vcf: str, output_dir: str, chunk_size: int = 20000) -> list:
    """
    Splits an input VCF into indexed chunks of `chunk_size` records each.
    Returns a list of created chunk file paths.
    """
    if not os.path.exists(input_vcf):
        raise FileNotFoundError(f"Input VCF not found: {input_vcf}")

    os.makedirs(output_dir, exist_ok=True)

    in_vcf = pysam.VariantFile(input_vcf)
    header = in_vcf.header

    chunk_idx = 0
    record_count = 0
    current_chunk_records = 0
    created_chunks = []

    current_out_vcf_path = None
    current_out_vcf = None

    def start_new_chunk(idx: int):
        nonlocal current_out_vcf_path, current_out_vcf, current_chunk_records
        chunk_name = f"chunk_{idx:04d}.vcf"
        current_out_vcf_path = os.path.join(output_dir, chunk_name)
        current_out_vcf = pysam.VariantFile(current_out_vcf_path, "w", header=header)
        current_chunk_records = 0

    def close_and_compress_chunk():
        nonlocal current_out_vcf, current_out_vcf_path
        if current_out_vcf is not None:
            current_out_vcf.close()
            current_out_vcf = None
            
            # Compress and index with pysam tabix
            gz_path = f"{current_out_vcf_path}.gz"
            pysam.tabix_compress(current_out_vcf_path, gz_path, force=True)
            if os.path.exists(current_out_vcf_path):
                os.remove(current_out_vcf_path)
            pysam.tabix_index(gz_path, preset="vcf", force=True)
            created_chunks.append(gz_path)

    for record in in_vcf:
        if current_out_vcf is None:
            start_new_chunk(chunk_idx)

        current_out_vcf.write(record)
        record_count += 1
        current_chunk_records += 1

        if current_chunk_records >= chunk_size:
            close_and_compress_chunk()
            chunk_idx += 1

    # Close any open trailing chunk
    if current_out_vcf is not None:
        close_and_compress_chunk()

    in_vcf.close()

    print(f"[VCF Chunker] Split {record_count} records from {input_vcf} into {len(created_chunks)} chunks (size ~{chunk_size}) in {output_dir}")
    return created_chunks


def main():
    parser = argparse.ArgumentParser(description="Split VCF into indexed chunks for parallel processing.")
    parser.add_argument("--input", "-i", required=True, help="Input VCF or VCF.gz file")
    parser.add_argument("--output-dir", "-o", required=True, help="Output directory to store chunk files")
    parser.add_argument("--chunk-size", "-c", type=int, default=20000, help="Number of variant records per chunk (default: 20000)")
    args = parser.parse_args()

    chunks = split_vcf_into_chunks(args.input, args.output_dir, args.chunk_size)
    for c in chunks:
        print(c)


if __name__ == "__main__":
    main()
