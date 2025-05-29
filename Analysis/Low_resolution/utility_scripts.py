# Utility scripts 09/02/2024
# This script contains utility functions for working with AnnData objects,
# including adding spatial coordinates, normalizing them, and calculating
# the percentage of different gene biotypes.

import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import issparse
from typing import Optional
import anndata as ad
import json
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.sparse import issparse
from scipy.spatial import distance_matrix

def add_spatial_coordinates(adata, csv_path):
    """
    Adds spatial coordinates from a CSV file to an AnnData object.

    Parameters:
    - adata: AnnData object to which the spatial coordinates will be added.
    - csv_path: Path to the CSV file containing cell names and spatial coordinates.
    
    Returns:
    - AnnData: The original AnnData object with spatial coordinates added to obsm['spatial'].
    """
    # Read the CSV file containing cell names and coordinates
    spatial_df = pd.read_csv(csv_path, sep='\t', index_col=0, names=['x', 'y'])

    # Merge the spatial information with the original AnnData object
    adata_spatial = pd.merge(adata.obs, spatial_df, left_index=True, right_index=True, how='left')

    # Handle NaN values for cells without coordinates
    adata_spatial[['x', 'y']] = adata_spatial[['x', 'y']].where(pd.notna(adata_spatial[['x', 'y']]), np.nan)

    # Print a message about the number of cells with added coordinates
    print(f"Added spatial coordinates for {adata_spatial[['x', 'y']].count().min()} cells.")

    # Add the 'spatial' component to the AnnData object
    adata.obsm['spatial'] = adata_spatial[['x', 'y']].values
    nan_indices = np.isnan(adata.obsm['spatial']).any(axis=1)

    # Drop instances where any element is NaN
    adata = adata[~nan_indices]

    return adata


def normalize_spatial_coordinates(adata, coordinates_file):
    """
    Normalizes the spatial coordinates in an AnnData object to a range of [0, 6500].

    Parameters:
    - adata: AnnData object containing spatial coordinates.
    - coordinates_file: Path to the file with spatial coordinates to add and normalize.

    Returns:
    - AnnData: The AnnData object with normalized spatial coordinates in obsm['spatial'].
    """
    adata = add_spatial_coordinates(adata, coordinates_file)
    adata.obsm['spatial'][:, 1] = (adata.obsm['spatial'][:, 1] - adata.obsm['spatial'][:, 1].min()) / adata.obsm['spatial'][:, 1].max() * 6500
    adata.obsm['spatial'][:, 0] = (adata.obsm['spatial'][:, 0] - adata.obsm['spatial'][:, 0].min()) / adata.obsm['spatial'][:, 0].max() * 6500
    return adata


def high_exp_list(adata, n_top=30, gene_symbols=None, sort_by='mean'):
    """
    Identify the top n genes with the highest mean or sum expression.

    Parameters:
    - adata: AnnData object.
    - n_top: Number of top genes to return (default is 30).
    - gene_symbols: Key for gene symbols if different from var_names.
    - sort_by: Sorting method, either 'mean' or 'sum' (default is 'mean').

    Returns:
    - List of top gene names.
    """
    # Compute the percentage of each gene per cell
    norm_dict = sc.pp.normalize_total(adata, target_sum=100, inplace=False)
    
    # Choose sorting method
    if sort_by not in ['mean', 'sum']:
        raise ValueError("sort_by must be either 'mean' or 'sum'")
    
    if sort_by == 'mean':
        # Identify the genes with the highest mean
        if issparse(norm_dict['X']):
            metric = norm_dict['X'].mean(axis=0).A1
        else:
            metric = norm_dict['X'].mean(axis=0)
    elif sort_by == 'sum':
        # Identify the genes with the highest sum
        if issparse(norm_dict['X']):
            metric = norm_dict['X'].sum(axis=0).A1
        else:
            metric = norm_dict['X'].sum(axis=0)
    
    # Get the indices of the top genes
    top_idx = np.argsort(metric)[::-1][:n_top]
    
    # Get the gene names
    gene_names = (
        adata.var_names[top_idx]
        if gene_symbols is None
        else adata.var[gene_symbols][top_idx]
    )

    return list(gene_names)


def add_biotypes_pct(
    adata: ad.AnnData,
    biomart: Optional[pd.DataFrame] = None,
    gene_colname: str = "GeneSymbol",
    biotype_colname: str = "Biotype",
    prefix: str = "pct.",
    scale: int = 100,
    verbose: bool = True
) -> ad.AnnData:
    """
    Adds gene biotype percentage values to an AnnData object.

    Parameters:
    - adata: AnnData object containing gene expression data.
    - biomart: DataFrame containing gene biotypes.
    - gene_colname: Column name in biomart DataFrame for gene identifiers (default is "GeneSymbol").
    - biotype_colname: Column name in biomart DataFrame for biotype (default is "Biotype").
    - prefix: Prefix for column names added to the AnnData object (default is "pct.").
    - scale: Scaling factor for the percentage (default is 100).
    - verbose: If True, print messages during function execution (default is True).

    Returns:
    - AnnData: The AnnData object with added gene biotype percentage values.
    """
    if biomart is None:
        if verbose: 
            print("Need a list of gene biotypes! Nothing done.")
        return adata

    biotypes = biomart[biotype_colname].unique()

    for biotype in biotypes:
        tmp_mart = biomart[biomart[biotype_colname] == biotype]
        tmp_feat = tmp_mart[gene_colname][tmp_mart[gene_colname].isin(adata.var_names)].unique()

        if len(tmp_feat) == 0:
            if verbose: 
                print(f"No {biotype} genes found...")
            continue

        col_name = prefix + biotype
        gene_pct = adata[:, tmp_feat].X.sum(axis=1) / adata.X.sum(axis=1)

        if scale == 100: 
            gene_pct *= 100
        elif scale != 1:
            if verbose:
                print(f"Given scale value '{scale}' is not valid. Defaulting to scale=100.")
            gene_pct *= 100

        adata.obs[col_name] = np.asarray(gene_pct).flatten()

    return adata

## Add location classification from manual spot annotations json
def add_histology_info_visium(json_files, barcode_file):
    # Read barcode information from text file
    barcode_data = {}
    with open(barcode_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            barcode_data[(int(parts[1]), int(parts[2]))] = parts[0]

    # Initialize dictionary to store barcode presence in each region
    barcode_presence = {barcode: {'A': False, 'B': False, 'C': False} for barcode in barcode_data.values()}

    # Process each JSON file
    for region, json_file in enumerate(json_files, start=ord('A')):
        # Load JSON data from file
        with open(json_file, 'r') as f:
            json_data = json.load(f)
        
        # Filter the list of dictionaries to include only those with tissue: True
        tissue_data = json_data['oligo']
        tissue_data = [d for d in tissue_data if d.get('tissue') == True]

        # Update barcode presence based on tissue data
        for d in tissue_data:
            row_col_pair = (d['col'], d['row'])
            if row_col_pair in barcode_data:
                barcode = barcode_data[row_col_pair]
                barcode_presence[barcode][chr(region)] = True

    # Ensure all barcodes are included in the DataFrame, with False where not present
    all_barcodes = list(barcode_data.values())
    for barcode in all_barcodes:
        if barcode not in barcode_presence:
            barcode_presence[barcode] = {'A': False, 'B': False, 'C': False}

    # Convert barcode presence dictionary to DataFrame
    df = pd.DataFrame.from_dict(barcode_presence, orient='index')
    
    return df

# Functions to extract microbial counts for a given experiment dict 
# Experiment dictionary
# Experiment_dict = {
#     "CTL": ["Visium_heart", "STRS_heart"], 
#     "PS": ["Visium_A", "STRS_A"],
#     "IL": ["Visium_B", "STRS_B"], 
#     "CEC": ["Visium_D", "STRS_D"],
#     "CO": ["Visium_C", "STRS_C"]
# }
def extract_data(mudata_dict, experiment_dict):
    data_list = []

    for experiment, samples in experiment_dict.items():
        for sample in samples:
            mudata = mudata_dict[sample]  # Retrieve MuData object
            sample_type = 'STRS' if 'STRS' in sample else 'Visium'

            # Extract 'Total_counts' from 'family' modality and 'Total_molecules_detected' from 'host' modality
            genus_total_counts = mudata.mod['family'].obs['Total_counts_classified']
            host_total_molecules_detected = mudata.mod['host'].obs['Total_molecules_detected']
            
            # Create a DataFrame with the data
            combined_df = pd.DataFrame({
                'Total_counts': genus_total_counts,
                'Total_molecules_detected': host_total_molecules_detected
            })
            
            # Add experiment and type information
            combined_df['Experiment'] = experiment
            combined_df['Type'] = sample_type
            combined_df['Sample'] = sample

            data_list.append(combined_df)

    combined_df = pd.concat(data_list, axis=0)
    return combined_df

# Function to filter data based on classification classification C is for spots covered with Tissue based on the histology
def filter_data(combined_df, mudata_dict):
    filtered_data_list = []

    for sample in combined_df['Sample'].unique():
        mudata = mudata_dict[sample]  # Retrieve MuData object

        # Filter the values where classification_C is True
        classification_C = mudata.mod['host'].obsm['classification']['C']
        filter_condition = classification_C
        
        # Apply the filter condition to the DataFrame
        sample_data = combined_df[combined_df['Sample'] == sample]
        filtered_data = sample_data[filter_condition.values]
        
        filtered_data_list.append(filtered_data)

    filtered_combined_df = pd.concat(filtered_data_list, axis=0)
    return filtered_combined_df
## Functions to extract and filter and make my life easier for now -- for the boxplots
def extract_feature_data(mudata_dict, experiment_dict, feature):
    data_list = []

    for experiment, samples in experiment_dict.items():
        for sample in samples:
            mudata = mudata_dict[sample]  # Retrieve MuData object
            sample_type = 'STRS' if 'STRS' in sample else 'Visium'

            # Extract the feature from 'host' modality
            feature_data = mudata.mod['host'].obs[feature]
            
            # Create a DataFrame with the feature data
            combined_df = pd.DataFrame({
                feature: feature_data
            })
            
            # Add experiment and type information
            combined_df['Experiment'] = experiment
            combined_df['Type'] = sample_type
            combined_df['Sample'] = sample

            data_list.append(combined_df)

    combined_df = pd.concat(data_list, axis=0)
    return combined_df
## Filter the data to keep spots covered by Tissue  (~classification_B) & (classification_C)
def filter_feature_data(combined_df, mudata_dict, feature):
    filtered_data_list = []

    for sample in combined_df['Sample'].unique():
        mudata = mudata_dict[sample]
        host_obs = mudata.mod['host'].obs

        # Get classification filters (must be same index as host.obs)
        classification_B = mudata.mod['host'].obsm['classification']['B']
        classification_C = mudata.mod['host'].obsm['classification']['C']
        filter_condition = (~classification_B) & (classification_C)

        # Get the subset of combined_df corresponding to this sample
        sample_df = combined_df[combined_df['Sample'] == sample].copy()

        # Align the sample_df to the index of host.obs (which is also index of filter_condition)
        # This assumes combined_df uses the same index as mudata.mod['host'].obs
        aligned_df = sample_df.loc[filter_condition.index.intersection(sample_df.index)]

        # Now apply the filter condition
        filtered_data = aligned_df.loc[filter_condition.loc[aligned_df.index]]
        filtered_data_list.append(filtered_data)

    filtered_combined_df = pd.concat(filtered_data_list, axis=0).reset_index(drop=True)
    return filtered_combined_df

def calculate_statistics(filtered_combined_df, feature):
    # Calculate statistics for each sample
    sample_stats = filtered_combined_df.groupby(['Experiment', 'Type', 'Sample'])[feature].agg(['mean', 'median', 'std']).reset_index()
    
    # Calculate statistics across Type (excluding Control)
    group_stats = filtered_combined_df[filtered_combined_df['Experiment'] != "CTL"].groupby('Type')[feature].agg(['mean', 'median']).reset_index()
    
    return sample_stats, group_stats

def calculate_distances_and_bins(mudata_dict, sample, modality, num_bins):
    """
    Assigns each spatial spot to a bin based on its distance to the lumen (C but not B spots are distance=0).
    """
    adata = mudata_dict[sample][modality]
    classification_B = adata.obsm['classification']['B']
    classification_C = adata.obsm['classification']['C']

    coords = adata.obsm['spatial']
    B_coords = coords[classification_B]
    C_not_B_coords = coords[classification_C & ~classification_B]

    if len(B_coords) > 0 and len(C_not_B_coords) > 0:
        dist_matrix = distance_matrix(B_coords, C_not_B_coords)
        min_distances = dist_matrix.min(axis=1)

        B_barcodes = adata.obs.index[classification_B]
        distance_B_df = pd.DataFrame(min_distances, index=B_barcodes, columns=['min_distance'])

        C_not_B_barcodes = adata.obs.index[classification_C & ~classification_B]
        distance_C_not_B_df = pd.DataFrame(0, index=C_not_B_barcodes, columns=['min_distance'])

        distance_df = pd.concat([distance_B_df, distance_C_not_B_df])
        distance_df['min_distance'] = np.clip(distance_df['min_distance'], a_min=0, a_max=None)

        distance_df['bin'] = pd.cut(distance_df['min_distance'], bins=num_bins-1, labels=range(1, num_bins))
        distance_df['bin'] = distance_df['bin'].astype('category').cat.set_categories(range(num_bins))
        distance_df.loc[C_not_B_barcodes, 'bin'] = 0

        adata.obs['bin'] = distance_df['bin']
        return adata
    else:
        return adata
    
# Calculate relative abundances
def calculate_relative_abundances(adata, genes):
    adata = adata[adata.obsm['classification']['C']]
    X = adata.X  # Expression data matrix
    var_names = adata.var_names  # Gene names
    expression_data = pd.DataFrame(X.toarray() if issparse(X) else X, index=adata.obs.index, columns=var_names)
    expression_data['location'] = adata.obs['location']
    aggregated_expression = expression_data.groupby('location').sum()
    proportions = aggregated_expression.div(aggregated_expression.sum(axis=1), axis=0).fillna(0)
    
    # Filter out taxa below 2% and group them into 'Other'
    if 'unclassified' in proportions.columns:
        proportions['Other'] = proportions.pop('unclassified')
    else:
        proportions['Other'] = 0

    taxa_to_keep = proportions.columns[(proportions.max(axis=0) >= 0.02)]
    proportions['Other'] += 1 - proportions[taxa_to_keep].sum(axis=1)
    proportions = proportions[taxa_to_keep]
    
    return proportions

# Extract unique taxa per spot
def get_unique_taxa_data(mudata_dict, samples):
    unique_taxa_data = {}
    for sample in samples:
        ## Add a filter to exclude lowly occuring taxa before calculating Richness
        adata = mudata_dict[sample]['genus'].copy()
        adata.X = adata.raw.X
        adata = adata[adata.obsm['classification']['C']]
        taxa_counts = adata.X.sum(axis=0) / adata.obs['Total_counts'].sum()
        mask = taxa_counts > 0.0001
        adata = adata[:,mask]
        spots = adata.X > 0 
        adata.obs['unique_taxa_per_spot'] = spots.sum(axis=1)
        
 
        unique_taxa_data[sample] = adata.obs['unique_taxa_per_spot'].values
    return unique_taxa_data


# Example of how to import and use these functions from a .py script
# Save the script as utils.py

# To import and use these functions in another script or interactive session:
# from utils import add_spatial_coordinates, normalize_spatial_coordinates, high_exp_list, add_biotypes_pct

# Example usage:
# adata_with_spatial = add_spatial_coordinates(adata, 'path/to/coordinates_file.csv')
# normalized_adata = normalize_spatial_coordinates(adata, 'path/to/coordinates_file.csv')
# top_genes = high_exp_list(adata, n_top=50, sort_by='sum')
# adata_with_biotypes = add_biotypes_pct(adata, biomart_df)
# Function to extract data from MuData objects
