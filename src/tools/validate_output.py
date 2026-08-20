import os
import glob
import argparse
import pandas as pd

def check_file_status(file_path):
    try:
        # Try to read the first few rows to quickly check if it's readable
        if file_path.endswith('.csv'):
            df = pd.read_csv(file_path, nrows=5)
        elif file_path.endswith('.pq') or file_path.endswith('.parquet'):
            df = pd.read_parquet(file_path)
            df = df.head(5)
        else:
            df = pd.read_excel(file_path, nrows=5)
        if df.empty:
            return "EMPTY"
        return "VALID"
    except Exception:
        return "CORRUPTED"

def validate_pipeline(raw_dir):
    raw_folder = os.path.basename(os.path.normpath(raw_dir))
    
    tmp_master_dir = f"_tmp/{raw_folder}"
    annotation_master_dir = f"annotation/{raw_folder}_annotated"
    gnomad_genes_dir = "annotation/gnomAD_subset"
    results_master_dir = f"results/{raw_folder}"
    
    report = []
    report.append(f"Validation Report for {raw_folder}")
    report.append("-" * 40)
    
    # Collect all input files
    input_files = glob.glob(os.path.join(raw_dir, "*.xlsx"))
    
    if not input_files:
        report.append(f"No .xlsx input files found in {raw_dir}")
        print(report[-1])
        return
    
    genes_missing_files = 0
    total_genes = len(input_files)

    for input_file in input_files:
        base_name = os.path.basename(input_file).split('.')[0]
        gene_name = base_name.split('_')[0]
        
        report.append(f"\nGene: {gene_name} (Input: {os.path.basename(input_file)})")
        
        expected_files = {
            "Step 1 (variant_converter)": os.path.join(tmp_master_dir, f"{gene_name}.vcf.gz"),
            "Step 2 (gnomAD_subset)": os.path.join(gnomad_genes_dir, f"gnomAD.v4.1.{gene_name}.vcf.gz"),
            "Step 3 (merge_gnomAD)": os.path.join(annotation_master_dir, f"{gene_name}.merge.gnomAD.vcf.gz"),
            "Step 4.1 (annotate_spip)": os.path.join(annotation_master_dir, f"{gene_name}.annSPiP.vcf.gz"),
            "Step 4.2 (annotate_vep)": os.path.join(annotation_master_dir, f"{gene_name}.annVEP.vcf.gz"),
            "Step 5 (merge_vep_spip)": os.path.join(annotation_master_dir, f"{gene_name}.annotated.vcf.gz"),
        }
        
        all_present = True
        for step, expected_file in expected_files.items():
            if not os.path.exists(expected_file):
                report.append(f"  [MISSING] {step}: {expected_file}")
                all_present = False

        # Step 6 (tsv or parquet)
        step6_name = "Step 6 (vcf2tsv)"
        expected_tsv = os.path.join(results_master_dir, f"{gene_name}.parsed.tsv")
        expected_parquet = os.path.join(results_master_dir,'compressed_tsv_parquet', f"{gene_name}.parsed.parquet")
        if os.path.exists(expected_tsv):
            step6_file = expected_tsv
        elif os.path.exists(expected_parquet):
            step6_file = expected_parquet
        else:
            step6_file = None

        if not step6_file:
            report.append(f"  [MISSING] {step6_name}: Neither .tsv nor .parquet found in {results_master_dir}.")
            all_present = False
        # Step 7 (final output)
        step7_name = "Step 7 (Final clean output)"
        expected_xlsx = os.path.join(results_master_dir, f"{gene_name}.parsed.clean.xlsx")
        expected_csv = os.path.join(results_master_dir, f"{gene_name}.parsed.clean.csv")
        expected_pq = os.path.join(results_master_dir, f"{gene_name}.parsed.clean.pq")
        
        if os.path.exists(expected_xlsx):
            final_file = expected_xlsx
        elif os.path.exists(expected_csv):
            final_file = expected_csv
        elif os.path.exists(expected_pq):
            final_file = expected_pq
        else:
            final_file = None
            
        if not final_file:
            report.append(f"  [MISSING] {step7_name}: Neither .xlsx, .csv nor .pq found in {results_master_dir}.")
            all_present = False
        else:
            # Extra check for the final file to ensure it is not corrupted
            status = check_file_status(final_file)
            if status == "CORRUPTED":
                report.append(f"  [CORRUPTED] {step7_name}: {final_file} is corrupted or not valid.")
                all_present = False
            elif status == "EMPTY":
                report.append(f"  [WARNING] {step7_name}: {final_file} is valid but contains no data rows (only header).")

        if all_present:
            report.append("  -> All steps completed successfully.")
        else:
            genes_missing_files += 1

    report.append(f"\nSummary:")
    report.append(f"Total genes evaluated: {total_genes}")
    report.append(f"Genes with missing or corrupted files: {genes_missing_files}")
    report.append(f"Genes fully processed: {total_genes - genes_missing_files}")
            
    report_text = "\n".join(report)
    
    # Save the report to a log file
    report_filename = f"results/{raw_folder}/validation_report_{raw_folder}.txt"
    os.makedirs(os.path.dirname(report_filename), exist_ok=True)
    with open(report_filename, "w") as f:
        f.write(report_text)
        
    print(f"Validation complete. Findings have been saved to {report_filename}")
    print(f"({genes_missing_files} out of {total_genes} genes had missing/corrupted files)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate output of pipeline steps.")
    parser.add_argument("--raw-dir", required=True, help="Path to the raw directory containing input .xlsx files")
    args = parser.parse_args()
    
    validate_pipeline(args.raw_dir)
