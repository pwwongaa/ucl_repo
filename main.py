import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ttest_ind  # t-test

# from statsmodels.stats.multitest import fdrcorrection

#configuration
EXPR_FILE = "data/GSE138852_pseudobulk_astro_counts.csv"
META_FILE = "data/GSE138852_pseudobulk_astro_metadata.csv"
OUTPUT_DEG_CSV = "DEG_AD_vs_CT.csv"
OUTPUT_DEG_AD_CSV = "DEG_AD_vs_CT_ADgene.csv"
AD_RELEVANT_GENES = ["APOE", "CLU", "TREM2", "BIN1", "PICALM", 
                    "MAPT", "PSEN1", "PSEN2", "APP", "CR1"]

#?helper functions
#expression matrix: genes VS samples_id
def load_expression_data(expr_file: str) -> pd.DataFrame: #type hints
    """Load and preprocess expression table"""
    if not os.path.exists(expr_file):
        raise FileNotFoundError(f"file not found: {expr_file}")
    df = pd.read_csv(expr_file, index_col=0)
    df = df.dropna(how='any') #dropna row
    print(f" Expression matrix loaded. Shape: {df.shape}")
    return df

#meta data: samples_id VS labels - AD and Health group
def load_metadata(meta_file: str) -> pd.DataFrame:
    """Load and preprocess metadata table"""
    if not os.path.exists(meta_file):
        raise FileNotFoundError(f"file not found: {meta_file}")
    df = pd.read_csv(meta_file, index_col=0)
    df['batchCond'] = df['batchCond'].str.upper() #consistent uppercase
    df = df.dropna(how='any') #dropna row
    print(f" Metadata loaded. Shape: {df.shape}")
    return df

#Differential expressed genes
def calculate_deg(expr_df: pd.DataFrame, meta_df: pd.DataFrame) -> pd.DataFrame:
    """To perform DE analysis using t-test and log2 fold change"""
    # Identify sample and control groups
    ad_samples = meta_df[meta_df['batchCond'] == 'AD'].index.intersection(expr_df.columns).tolist()
    ct_samples = meta_df[meta_df['batchCond'] == 'CT'].index.intersection(expr_df.columns).tolist()

    print(f"AD samples: {ad_samples}")
    print(f"CT samples: {ct_samples}")

    # Analyse each gene: AD vs CT
    results = []
    for gene in expr_df.index:
        ad_vals = expr_df.loc[gene, ad_samples].astype(float)
        ct_vals = expr_df.loc[gene, ct_samples].astype(float)

        #?compute Log2FC: AD - CT,healthy control
        ##Positive values → gene is upregulated in AD; negative → downregulated.
        log2fc = np.log2(ad_vals.mean() + 1e-6) - np.log2(ct_vals.mean() + 1e-6) #+1e-6 to avoid log0

        #?perform t-test
        stat, pval = ttest_ind(ad_vals, ct_vals, equal_var=False)

        results.append((gene, log2fc, pval))

    #aggregate the result into df:deg_df & output as a CSV file
    ##pval: p-value
    deg_df = pd.DataFrame(results, columns=["gene", "log2FC", "pval"])
    deg_df["adj_pval"] = deg_df["pval"] * len(deg_df)  #* p-value correction, control false positive
    deg_df = deg_df.sort_values(by="log2FC", ascending=False) #sorting
    deg_df.to_csv(OUTPUT_DEG_CSV, index=False)
    print(f" DEG results saved: {deg_df.shape[0]} genes tested.")
    return deg_df

#VOLCANO plot
def plot_volcano(deg_df: pd.DataFrame):
    """Draw volcano plot. - to visualise the UP and DOWN regulated genes"""
    deg_df['adj_pval'] = deg_df['adj_pval'].replace(0, 1e-300) # replce 0
    #label significance
    deg_df['significance'] = 'Not significant'
    deg_df.loc[(deg_df['log2FC'] > 1) & (deg_df['pval'] < 0.05), 'significance'] = 'Up'
    deg_df.loc[(deg_df['log2FC'] < -1) & (deg_df['pval'] < 0.05), 'significance'] = 'Down'

    #plot all genes, with grouping the signficant DEGs in UP and DOWN groups
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=deg_df, x='log2FC', y=-np.log10(deg_df['adj_pval']),
                    hue='significance', palette={'Up': 'red', 'Down': 'blue', 'Not significant': 'grey'},
                    alpha=0.7)
    plt.axvline(1, color='black', linestyle='--') #Add threshold lines
    plt.axvline(-1, color='black', linestyle='--') #Add threshold lines
    plt.axhline(-np.log10(0.05), color='black', linestyle='--')

    plt.title("Volcano Plot: AD vs Control")
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-Log10 Adjusted P-Value")
    plt.legend(title="Significance")
    plt.tight_layout()
    plt.show()

def plot_heatmap(deg_df: pd.DataFrame, expr_df: pd.DataFrame):
    """Plot heatmap for top genes ."""
    # Extract top genes
    deg_df['adj_pval'] = deg_df['adj_pval'].replace(0, 1e-300)
    up_genes = deg_df[(deg_df['log2FC'] > 1) & (deg_df['pval'] < 0.05)].sort_values(by="log2FC", ascending=False).head(10)
    down_genes = deg_df[(deg_df['log2FC'] < -1) & (deg_df['pval'] < 0.05)].sort_values(by="log2FC").head(10)
    top_combined = pd.concat([up_genes, down_genes])

    # Heatmap - top genes
    selected_genes = top_combined["gene"].tolist()
    heatmap_data = expr_df.loc[selected_genes]

    # Z-score normalisation
    heatmap_data_z = (heatmap_data - heatmap_data.mean(axis=1).values[:, None]) / heatmap_data.std(axis=1).values[:, None]

    plt.figure(figsize=(12, 8))
    sns.heatmap(heatmap_data_z, cmap="vlag", xticklabels=True, yticklabels=True, cbar_kws={"label": "Z-score"})
    plt.title("Expression Heatmap of Top DEGs")
    plt.xlabel("Sample")
    plt.ylabel("Gene")
    plt.tight_layout()
    plt.show()

def plot_heatmap_adgene(deg_df: pd.DataFrame, expr_df: pd.DataFrame, focus_genes: list = None):
    """Plot heatmap for selected AD-related genes ."""
    deg_df['adj_pval'] = deg_df['adj_pval'].replace(0, 1e-300)
    
    #AD related genes
    if focus_genes:
        top_combined = deg_df[deg_df["gene"].isin(focus_genes)]
    else:
        up_genes = deg_df[(deg_df['log2FC'] > 1) & (deg_df['pval'] < 0.05)].sort_values(by="log2FC", ascending=False).head(10)
        down_genes = deg_df[(deg_df['log2FC'] < -1) & (deg_df['pval'] < 0.05)].sort_values(by="log2FC").head(10)
        top_combined = pd.concat([up_genes, down_genes])

    selected_genes = top_combined["gene"].tolist()
    # Heatmap
    heatmap_data = expr_df.loc[selected_genes]
    # Z-score normalisation
    heatmap_data_z = (heatmap_data - heatmap_data.mean(axis=1).values[:, None]) / heatmap_data.std(axis=1).values[:, None]

    plt.figure(figsize=(12, 8))
    sns.heatmap(heatmap_data_z, cmap="vlag", xticklabels=True, yticklabels=True, cbar_kws={"label": "Z-score"})
    plt.title("Expression Heatmap of Selected Genes")
    plt.xlabel("Sample")
    plt.ylabel("Gene")
    plt.tight_layout()
    plt.show()


#? Main function
def main():
    # I. Data loading and preprocessing
    #i. load expression profile table: gene vs samples
    expr_df = load_expression_data(EXPR_FILE)
    # expr_df = sort_each_numeric_column_desc(expr_df)
    print(expr_df.head())
    print(expr_df.describe())
    print(f" Shape: {expr_df.shape}")

    #ii. load sample metadata: samples vs condition
    meta_df = load_metadata(META_FILE)
    print(meta_df.head())

    #iii. match the sample id with expression table
    unmatched = set(expr_df.columns) - set(meta_df.index)
    if unmatched != set(): 
        print(f"unmatched sample if any: {unmatched}")

    #iv. DEG
    deg_df = calculate_deg(expr_df, meta_df)
    # print(deg_df)
    #* focus on common AD-related genes: 
    deg_filtered = deg_df[deg_df["gene"].isin(AD_RELEVANT_GENES)].copy()
    deg_filtered.to_csv(OUTPUT_DEG_AD_CSV)
    print("Filtered AD-relevant genes:")
    print(deg_filtered)

    #v. plot all graphs
    plot_volcano(deg_df) #volcano: UP and DOWN regulation
    plot_heatmap(deg_df, expr_df) #heatmap: AD vs control groups
    plot_heatmap_adgene(deg_df, expr_df, focus_genes=AD_RELEVANT_GENES) #heatmap: AD vs control groups


#! Main execution
if __name__ == "__main__":
    main()

