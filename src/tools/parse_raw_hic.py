# Script to create one file per gene called `GENENAME_enrichment.xlsx`
import os
import pandas as pd
import shutil
import csv
import mygene

def parse_hic_data(src_dir='_RAW/raw_data_HIC/', output_dir='_RAW/data_HIC/'):
    input_xlsx = os.path.join(src_dir, 'Enrichment_45genes.xlsx')
    transcript_file = os.path.join(output_dir, 'data_hic_transcripts.txt')
    mapping_path = 'resources/gene_transcript_mapping.txt'

    if not os.path.exists(input_xlsx):
        print(f"Warning: Input file {input_xlsx} does not exist. Exiting.")
        return

    print(f'Creating output directory: {output_dir}')
    os.makedirs(output_dir, exist_ok=True)

    print(f'Reading Excel file: {input_xlsx}')
    df = pd.read_excel(input_xlsx)
    print('Excel file loaded.')
    if 'GEN' not in df.columns:
        raise ValueError("Column 'GEN' not found in the Excel file.")
    print(f"Found {df['GEN'].nunique()} unique genes in the data.")
    for gene, gene_df in df.groupby('GEN'):
        gene = str(gene)
        gene = gene.replace('/', '_').replace('_','').replace(' ', '-')
        print(f'Processing gene: {gene} (n={len(gene_df)})')
        out_path = os.path.join(output_dir, f'{gene}_enrichment.xlsx')
        if os.path.exists(out_path):
            print(f'Warning: {out_path} already exists, skipping write.')
            continue
        gene_df.to_excel(out_path, index=False)
        print(f'Wrote {out_path}')

    mg = mygene.MyGeneInfo()
    if 'NOMBRE_ADN' not in df.columns:
        raise ValueError("Column 'NOMBRE_ADN' not found in the Excel file.")
    print('Extracting transcript names from NOMBRE_ADN column...')
    transcripts = set(df['NOMBRE_ADN'].dropna().apply(lambda x: str(x).split(':')[0]))
    print(f'Found {len(transcripts)} unique transcript names.')

    print(f'Writing gene, transcript, and count to {transcript_file}')
    df_trans = df[['GEN', 'NOMBRE_ADN']].dropna().copy()
    df_trans['transcript'] = df_trans['NOMBRE_ADN'].apply(lambda x: str(x).split(':')[0])
    grouped = df_trans.groupby(['GEN', 'transcript']).size().reset_index(name='count')
    grouped = grouped.sort_values(['GEN', 'transcript'])
    grouped.to_csv(transcript_file, sep='\t', header=False, index=False)
    print(f'Saved gene, transcript, and count to {transcript_file}')

    print('Querying mygene.info for ENSEMBL transcript mappings...')
    results = mg.querymany(list(transcripts), scopes='refseq', fields='ensembl.transcript', species='human')
    print('Query complete.')

if __name__ == '__main__':
    parse_hic_data()

