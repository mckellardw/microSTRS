import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.sparse import issparse
from scipy.spatial import distance_matrix
import scanpy as sc
from scipy.stats import pearsonr

def dot_plot_normalized_zscore(adata, genes, sample_name, num_bins=5, min_counts=1, output_file=None):
    """
    Generate a dot plot showing percent and z-scored relative abundance of taxa or genes across spatial bins.
    Includes a 'Percent' gene to visually indicate the dot size scale.
    """
    adata = adata[adata.obsm['classification']['C']].copy()
    distance_df = adata.obs['bin']

    if distance_df.empty:
        print("No valid distances calculated.")
        return

    X = adata.raw.X
    var_names = adata.var.index
    gene_indices = [var_names.get_loc(g) for g in genes if g in var_names]
    if len(gene_indices) != len(genes):
        print("Missing genes:", set(genes) - set(var_names))

    cell_indices = adata.obs.index.isin(distance_df.index)

    if issparse(X):
        expression_data = pd.DataFrame(X[cell_indices, :][:, gene_indices].toarray(), index=adata.obs[cell_indices].index, columns=genes)
    else:
        expression_data = pd.DataFrame(X[cell_indices, :][:, gene_indices], index=adata.obs[cell_indices].index, columns=genes)

    gene_sums = expression_data.sum(axis=0)
    filtered_genes = gene_sums[gene_sums >= min_counts].index.tolist()
    if not filtered_genes:
        print("No genes meet threshold.")
        return

    expression_data = expression_data[filtered_genes]
    distance_df = distance_df.cat.codes.fillna(num_bins - 1).astype(int)
    expression_data['bin'] = distance_df
    aggregated_expression = expression_data.groupby('bin').sum().reindex(range(num_bins), fill_value=0)

    bin_counts = aggregated_expression.sum(axis=1)
    percent_counts = (aggregated_expression.T / bin_counts).T * 100
    z_scores = (percent_counts - percent_counts.mean()) / percent_counts.std()

    # Add a dummy "Percent" feature to visually encode dot size scale
    legend_gene = 'Percent'
    proportions = [1, 5, 20, 50, 100]
    for i in range(num_bins):
        percent_counts.loc[i, legend_gene] = proportions[i]
        z_scores.loc[i, legend_gene] = (proportions[i] - np.mean(proportions)) / np.std(proportions)

    clustered_genes = filtered_genes + [legend_gene]
    data = [{'Gene': g, 'Bin': b, 'Z-score': z_scores.loc[b, g], 'Percent': percent_counts.loc[b, g]} for g in clustered_genes for b in range(num_bins)]
    df = pd.DataFrame(data)

    # Plot
    fig_height = 0.4 * len(clustered_genes)
    fig, ax = plt.subplots(figsize=(2, fig_height))
    fig.patch.set_facecolor('none')
    ax.set_facecolor('none')

    scatter = sns.scatterplot(
        data=df,
        x='Bin',
        y='Gene',
        size='Percent',
        hue='Z-score',
        palette='coolwarm',
        edgecolor='k',
        ax=ax,
        sizes=(1,100),
        legend=False
    )

    # Set x-ticks explicitly to bin labels
    ax.set_xticks(range(num_bins))
    ax.set_xticklabels([str(i) for i in range(num_bins)], fontsize=7)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)
    plt.grid(False)
    plt.xlabel('Bin')
    plt.ylabel('')

    # Add border box
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(0.5)

    # Add colorbar
    norm = plt.Normalize(df['Z-score'].min(), df['Z-score'].max())
    sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.02)
    cbar.ax.tick_params(labelsize=6)

    if output_file:
        plt.savefig(output_file, format='pdf', dpi=300, transparent=True)
    plt.show()

def pseudobulk_comparisons(adata_dict, experiment_dict, modality):
    # Global style
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['figure.dpi'] = 300

    num_experiments = len(experiment_dict)
    fig, axs = plt.subplots(num_experiments, 1, figsize=(1.5, 1.4 * num_experiments))
    if num_experiments == 1:
        axs = [axs]

    for ax, (exp_name, (sample1, sample2)) in zip(axs, experiment_dict.items()):
        # Load and sum counts
        dfs = {}
        for sample in [sample1, sample2]:
            adata = adata_dict[sample][modality].copy()
            adata.X = adata.raw.X
            df = pd.DataFrame(adata.X.toarray().sum(axis=0), index=adata.var_names, columns=[sample])
            dfs[sample] = df

        # Merge, fill, normalize
        merged_df = pd.concat([dfs[sample1], dfs[sample2]], axis=1).fillna(0)
        mols_1 = adata_dict[sample1]['host'].obs['Total_molecules_detected'].sum()
        mols_2 = adata_dict[sample2]['host'].obs['Total_molecules_detected'].sum()
        tpm_df = merged_df.div([mols_1, mols_2], axis=1) * 1e6

        # Log-log scatter with pseudocount
        pseudo = 1
        x = tpm_df[sample1] + pseudo
        y = tpm_df[sample2] + pseudo
        max_val = max(x.max(), y.max())

        # Pearson correlation
        r, _ = pearsonr(np.log10(x), np.log10(y))

        # Plot
        ax.scatter(x, y, c='black', s=1, alpha=1, rasterized=True)
        ax.plot([pseudo, max_val], [pseudo, max_val], 'k--', lw=1)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(pseudo, max_val)
        ax.set_ylim(pseudo, max_val)

        ax.set_title(f'{exp_name} (r = {r:.2f})', fontsize=8)
        ax.set_xlabel('Visium', fontsize=7)
        ax.set_ylabel('Visium+PAP', fontsize=7)
        ax.tick_params(labelsize=7)

        # Hide minor ticks
        ax.xaxis.set_minor_locator(plt.NullLocator())
        ax.yaxis.set_minor_locator(plt.NullLocator())

        ax.legend().set_visible(False)

    plt.tight_layout()
    plt.show()
    
def plot_host_vs_microbe_scatter(mdata, sample_name, log=False):
    host = mdata['host']
    bacteria = mdata['family'] if 'family' in mdata.mod else None

    if bacteria is None:
        print(f"No bacterial modality for {sample_name}")
        return

    # Align indices and fill with zeros
    total_counts = host.obs['Total_counts'].reindex(host.obs_names).fillna(0)
    total_molecules = host.obs['Total_molecules_detected'].reindex(host.obs_names).fillna(0)
    microbial_counts = bacteria.obs['Total_counts_classified'].reindex(host.obs_names).fillna(0)

    # Classification masks
    is_B = host.obsm['classification']['B']
    is_C_not_B = host.obsm['classification']['C'] & ~is_B

    # Pearson r (exclude nan/inf values)
    valid = (~np.isnan(total_counts)) & (~np.isnan(microbial_counts)) & \
            (~np.isinf(total_counts)) & (~np.isinf(microbial_counts))
    r, _ = pearsonr(total_counts[valid], microbial_counts[valid])

    # Axis scaling
    x = total_counts
    y = microbial_counts
    if log:
        x = x + 1
        y = y + 1
        xscale = yscale = 'log'
        max_lim = max(x.max(), y.max())
        xlim = ylim = (1, max_lim)
    else:
        xscale = yscale = 'linear'
        xlim = (0, total_counts.max())
        ylim = (0, microbial_counts.max())

    # Plot
    plt.figure(figsize=(2, 2), dpi=300)
    plt.scatter(x[is_C_not_B], y[is_C_not_B], color='orange', s=2, label='Tissue', rasterized=True)
    plt.scatter(x[is_B], y[is_B], color='blue', s=2, label='Lumen', rasterized=True)

    plt.xlabel("Host Total_counts", fontsize=7)
    plt.ylabel("Microbial Counts", fontsize=7)
    plt.title(f"{sample_name} (r = {r:.2f})", fontsize=8)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.grid(False)
    plt.legend(fontsize=6, loc='upper left', frameon=False)
    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.xlim(*xlim)
    plt.ylim(*ylim)
    plt.tight_layout()
    plt.show()
