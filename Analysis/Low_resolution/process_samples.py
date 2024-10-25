import pandas as pd
import scanpy as sc
import numpy as np
import scipy.sparse as spp
import re
from multiprocessing import Pool
from tqdm import tqdm

print("Packages loaded successfully.")

def add_spatial_coordinates(adata, csv_path):
    print(f"Adding spatial coordinates from {csv_path}")
    spatial_df = pd.read_csv(csv_path, sep='\t', index_col=0, names=['x', 'y'])
    adata_spatial = pd.merge(adata.obs, spatial_df, left_index=True, right_index=True, how='left')
    adata_spatial[['x', 'y']] = adata_spatial[['x', 'y']].where(pd.notna(adata_spatial[['x', 'y']]), np.nan)
    adata.obsm['spatial'] = adata_spatial[['x', 'y']].values
    return adata

def process_sample(sample_info, host_unidentified=None):
    print(f"Processing sample: {sample_info['SampleID']}")
    kraken_df = pd.read_csv(sample_info['KrakenFilePath'], compression='gzip', sep='\t', header=None, names=['col1', 'col2', 'col3', 'col4', 'col5', 'col6', 'col7', 'col8'])
    print(f"Kraken data loaded for sample: {sample_info['SampleID']}")
    
    kraken_df['col7'] = kraken_df['col6'].astype(str) + kraken_df['col8'].astype(str)
    kraken_df['col3'] = kraken_df['col3'].apply(lambda x: 'taxid_' + re.search(r'taxid (\d+)', str(x)).group(1) if pd.notna(x) and re.search(r'taxid (\d+)', str(x)) else x)
    kraken_df = kraken_df.drop_duplicates(subset='col7')
    kraken_df = kraken_df[~kraken_df['col1'].isin(['U'])]
    kraken_df = kraken_df[~kraken_df['col7'].str.contains('-')]
    print(f"Data cleaned for sample: {sample_info['SampleID']}")
    
    cmtx = pd.crosstab(kraken_df['col6'], kraken_df['col3'])
    if host_unidentified:
        cmtx.drop(columns=host_unidentified, inplace=True)
    adata = sc.AnnData(spp.csr_matrix(cmtx))
    adata.var_names = cmtx.columns
    adata.obs_names = cmtx.index
    adata = add_spatial_coordinates(adata, sample_info['SpatialCoordsPath'])
    adata.write(sample_info['SavePath'])
    print(f"Processed microbial adata saved to {sample_info['SavePath']}.")

def main(csv_path, host_unidentified=None):
    print("Starting sample processing.")
    samples_df = pd.read_csv(csv_path)
    print("CSV file loaded successfully.")
    
    tasks = [(row, host_unidentified) for index, row in samples_df.iterrows()]

    with Pool(processes=24) as pool:  # Utilizing 24 cores
        for _ in tqdm(pool.starmap(process_sample, tasks), total=len(tasks)):
            pass

    print("All samples processed successfully.")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python script.py <csv_path> [optional: <host_unidentified_list>]")
        sys.exit(1)

    csv_path = sys.argv[1]
    host_unidentified = sys.argv[2].split(",") if len(sys.argv) > 2 else None
    
    main(csv_path, host_unidentified)

