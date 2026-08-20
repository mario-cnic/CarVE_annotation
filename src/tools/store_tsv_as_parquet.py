# script to transform all tsv in the results folder into parquet files, to save space and speed up downstream analyses
import pandas as pd
from pathlib import Path

def batch_convert_tsv_to_parquet(folder_path='./results/data_HIC', output_folder='compressed_tsv_parquet'):
    # Initialize paths
    src_dir = Path(folder_path)
    dest_dir = Path(folder_path) / output_folder
    dest_dir.mkdir(exist_ok=True)

    # Find all .tsv files
    tsv_files = list(src_dir.glob('*.tsv'))
    
    if not tsv_files:
        print("No .tsv files found in the directory.")
        return

    print(f"Found {len(tsv_files)} files. Starting conversion...")

    for tsv_path in tsv_files:
        try:
            # Construct output filename
            parquet_path = dest_dir / tsv_path.with_suffix('.parquet').name
            
            # Read TSV and write to Parquet
            # engine='pyarrow' is generally faster and more memory-efficient
            df = pd.read_csv(tsv_path, sep='\t', low_memory=False)
            df.to_parquet(parquet_path, engine='pyarrow', compression='snappy', index=False)
            
            print(f"Successfully converted: {tsv_path.name}")
        except Exception as e:
            print(f"Failed to convert {tsv_path.name}: {e}")

if __name__ == "__main__":
    batch_convert_tsv_to_parquet()