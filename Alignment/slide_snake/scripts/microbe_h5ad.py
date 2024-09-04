import sys
import gzip
import pandas as pd
import argparse
from scanpy import read_mtx
from numpy import intersect1d

# Example usage:
## python script.py \
#   --mat_in input_matrix.mtx \
#   --feat_in input_features.tsv \
#   --bc_in input_barcodes.txt \
#   --bb_map input_spatial_map.tsv \
#   --ad_out output_anndata.h5ad \
#   --feat_col 1 \
#   --remove_zero_features


def main(kraken_bcds, bb_map, ad_out_microbe, remove_zero_features=False):
    kraken_df = pd.read_csv('kraken_bcds', sep='\t', header=None, names=['col1', 'col2', 'col3', 'col4', 'col5', 'col6', 'col7', 'col8'])
    print(kraken_df.head())
    kraken_df['col7'] = kraken_df['col6'].astype(str) + kraken_df['col8'].astype(str)
    print(kraken_df.head())
    kraken_dmx = kraken_df.drop_duplicates(subset='col7')
    kraken_df.shape
    kraken_dmx.shape
    mask1 = (kraken_dmx['col1'] == 'U')
    kraken_dmx_filtered = kraken_dmx[~mask1]
    mask2 = kraken_dmx_filtered['col7'].str.contains('-')
    kraken_dmx_filtered = kraken_dmx_filtered[~mask2]
    kraken_dmx_filtered.shape
    cmtx = pd.crosstab(kraken_dmx_filtered['col6'],kraken_dmx_filtered['col3'])
    adata = ad.AnnData(spp.csr_matrix(cmtx))
    adata.var_names = cmtx.columns
    adata.obs_names = cmtx.index

    # Add spatial location
    spatial_data = pd.read_csv(bb_map, sep="\t", header=None, names=["barcode", "x", "y"])

    # Set the cell barcode as index
    spatial_data.set_index("barcode", inplace=True)

    # Check if the barcodes in the AnnData object match the ones in the spatial data
    if not all(adata.obs_names.isin(spatial_data.index)):
        print("Warning: Not all barcodes in the AnnData object match the ones in the spatial data.")

    # Add the spatial coordinates to the AnnData object
    common_labels = intersect1d(adata.obs_names, spatial_data.index)
    adata = adata[common_labels, :]
    spatial_coord = spatial_data.loc[common_labels, ['x', 'y']]

    print(f"{len(common_labels)} / {len(adata.obs_names)} found in barcode map...")

    # spatial_coord = spatial_data.loc[adata.obs_names, ['x', 'y']] #w/ out filtering for intersection
    adata.obsm['spatial'] = spatial_coord.to_numpy()

    # Remove observations with zero features detected
    if remove_zero_features:
        adata = adata[adata.X.sum(axis=1) > 0, :]

    # Write output
    adata.write(ad_out_microbe)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process spatial transcriptomics data.')
    parser.add_argument('--kraken_bcds', required=True, help='Input kraken_bcds (unzipped)')
    parser.add_argument('--bb_map', required=True, help='Input spatial map file (tsv format)')
    parser.add_argument('--ad_out_microbe', required=True, help='Output AnnData file (h5ad format)')
    parser.add_argument('--remove_zero_features', action='store_true', help='Remove observations with zero features detected (default: False)')

    args = parser.parse_args()

    main(args.kraken_bcds, args.bb_map, args.ad_out_microbe, args.remove_zero_features)


















os.chdir("/workdir/in68/microSTRS/Stereo-seq_Analysis/")
kraken_df = pd.read_csv('results_out/New_SS2_071723_solo/Kraken_barcodes.txt', sep='\t', header=None, names=['col1', 'col2', 'col3', 'col4', 'col5', 'col6', 'col7', 'col8'])
print(kraken_df.head())
kraken_df['col7'] = kraken_df['col6'].astype(str) + kraken_df['col8'].astype(str)
print(kraken_df.head())
kraken_dmx = kraken_df.drop_duplicates(subset='col7')
kraken_df.shape
kraken_dmx.shape
mask1 = (kraken_dmx['col1'] == 'U')
kraken_dmx_filtered = kraken_dmx[~mask1]
mask2 = kraken_dmx_filtered['col7'].str.contains('-')
kraken_dmx_filtered = kraken_dmx_filtered[~mask2]
kraken_dmx_filtered.shape
cmtx = pd.crosstab(kraken_dmx_filtered['col6'],kraken_dmx_filtered['col3'])
adata = ad.AnnData(spp.csr_matrix(cmtx))
adata.var_names = cmtx.columns
adata.obs_names = cmtx.index


# main(mat_in, feat_in, bc_in, bb_map, ad_out, feat_col=1, remove_zero_features=False):
#spatial_data = spatial_data_B2_uncommon
#adata = ad.AnnData(spp.csr_matrix(cmtx))
#adata.var_names = cmtx.columns
#adata.obs_names = cmtx.index
bb_mapdir = '/workdir/in68/microSTRS/data/A02077F2.txt'
spatial_data = pd.read_csv(bb_mapdir, sep="\t", header=None, names=["barcode", "x", "y"])

    # # Set the cell barcode as index
spatial_data.set_index("barcode", inplace=True)

    # # Check if the barcodes in the AnnData object match the ones in the spatial data
if not all(adata.obs_names.isin(spatial_data.index)):
    print("Warning: Not all barcodes in the AnnData object match the ones in the spatial data.")

    # # Add the spatial coordinates to the AnnData object
common_labels = intersect1d(adata.obs_names, spatial_data.index)
adata = adata[common_labels, :]
spatial_coord = spatial_data.loc[common_labels, ['x', 'y']]

print(f"{len(common_labels)} / {len(adata.obs_names)} found in barcode map...")

# spatial_coord = spatial_data.loc[adata.obs_names, ['x', 'y']] #w/ out filtering for intersection
adata.obsm['spatial'] = spatial_coord.to_numpy()
sto1_microbe_unk = adata 