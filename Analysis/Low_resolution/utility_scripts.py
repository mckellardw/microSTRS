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
