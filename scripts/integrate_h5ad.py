import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
import argparse
import scipy.sparse

# rRNA gene dictionaries with lengths
wMel_rRNA={
    "GQX67_00940": 2772,
    "GQX67_00945": 107,
    "GQX67_05945": 1505
}

Dmel_rRNA={
    "FBgn0085802": 1995,
    "FBgn0267504": 3970,
    "FBgn0267512": 123,
    "FBgn0267500": 30,
    "FBgn0267508": 821,
    "FBgn0267509": 123,
    "FBgn0267510": 30,
    "FBgn0267498": 1995,
    "FBgn0267502": 123,
    "FBgn0267503": 30,
    "FBgn0267505": 6863,
    "FBgn0267506": 4191,
    "FBgn0267511": 1721,
    "FBgn0085753": 2775,
    "FBgn0267507": 8118,
    "FBgn0267501": 1995,
    "FBgn0267496": 30,
    "FBgn0267497": 2715,
    "FBgn0267499": 123
}

Dmel_mt={
    "FBgn0013686": 1324,
    "FBgn0013688": 786
}

# Set plotting settings
sc.settings.set_figure_params(dpi=100, frameon=False)
sc.settings.verbosity = 1

def calculate_wolbachia_titer(adata):
    """
    Calculate Wolbachia titer for each cell in the AnnData object.
    The titer is calculated as the ratio of wMel rRNA counts per length to Dmel rRNA counts per length.
    """
    print("Calculating Wolbachia titer...")
    
    # Check if we need to look in var_names or gene_ids column
    if 'gene_ids' in adata.var.columns:
        # Create a mapping from gene names/indices to gene_ids
        gene_id_map = dict(zip(adata.var.index, adata.var['gene_ids']))
        
        # Find which genes are present in our dictionaries
        wMel_genes_present = []
        for gene_id in wMel_rRNA.keys():
            # Check if this gene ID is in the gene_ids column
            if gene_id in adata.var['gene_ids'].values:
                wMel_genes_present.append(gene_id)
                
        Dmel_genes_present = []
        for gene_id in Dmel_rRNA.keys():
            # Check if this gene ID is in the gene_ids column
            if gene_id in adata.var['gene_ids'].values:
                Dmel_genes_present.append(gene_id)
        
        # Create masks for the genes
        wMel_mask = [gene_id_map.get(idx, "") in wMel_genes_present for idx in adata.var.index]
        Dmel_mask = [gene_id_map.get(idx, "") in Dmel_genes_present for idx in adata.var.index]
    else:
        # Use original approach with var_names
        wMel_genes_present = [gene for gene in wMel_rRNA.keys() if gene in adata.var_names]
        Dmel_genes_present = [gene for gene in Dmel_rRNA.keys() if gene in adata.var_names]
        
        # Create masks based on var_names
        wMel_mask = [gene in wMel_genes_present for gene in adata.var_names]
        Dmel_mask = [gene in Dmel_genes_present for gene in adata.var_names]
    
    if not wMel_genes_present:
        print("Warning: No wMel rRNA genes found in the dataset!")
        adata.obs['wolbachia_titer'] = np.nan
        return adata
        
    if not Dmel_genes_present:
        print("Warning: No Dmel rRNA genes found in the dataset!")
        adata.obs['wolbachia_titer'] = np.nan
        return adata
    
    print(f"Found {len(wMel_genes_present)} wMel rRNA genes and {len(Dmel_genes_present)} Dmel rRNA genes")
    
    # Get gene indices from the masks
    wMel_indices = np.where(wMel_mask)[0]
    Dmel_indices = np.where(Dmel_mask)[0]
    
    # Convert sparse matrix to dense if necessary
    is_sparse = scipy.sparse.issparse(adata.X)
    
    # Calculate counts per length for each gene
    wMel_counts_per_length = np.zeros((adata.n_obs, len(wMel_indices)))
    Dmel_counts_per_length = np.zeros((adata.n_obs, len(Dmel_indices)))
    
    # Calculate counts per length for wMel genes
    for i, idx in enumerate(wMel_indices):
        gene_idx = adata.var.index[idx]
        gene_id = gene_id_map.get(gene_idx, gene_idx) if 'gene_ids' in adata.var.columns else gene_idx
        length = wMel_rRNA[gene_id]
        
        if is_sparse:
            counts = adata.X[:, idx].toarray().flatten()
        else:
            counts = adata.X[:, idx]
            
        wMel_counts_per_length[:, i] = counts / length
    
    # Calculate counts per length for Dmel genes
    for i, idx in enumerate(Dmel_indices):
        gene_idx = adata.var.index[idx]
        gene_id = gene_id_map.get(gene_idx, gene_idx) if 'gene_ids' in adata.var.columns else gene_idx
        length = Dmel_rRNA[gene_id]
        
        if is_sparse:
            counts = adata.X[:, idx].toarray().flatten()
        else:
            counts = adata.X[:, idx]
            
        Dmel_counts_per_length[:, i] = counts / length
    
    # Calculate the mean counts per length for each organism and cell
    wMel_mean_per_cell = np.mean(wMel_counts_per_length, axis=1)
    Dmel_mean_per_cell = np.mean(Dmel_counts_per_length, axis=1)
    
    # Calculate the titer (ratio of wMel to Dmel rRNA counts per length)
    with np.errstate(divide='ignore', invalid='ignore'):
        titer = np.where(Dmel_mean_per_cell > 0, 
                         wMel_mean_per_cell / Dmel_mean_per_cell, 
                         np.nan)
    
    # Add the titer to the AnnData object
    adata.obs['wolbachia_titer'] = titer
    adata.obs['log1p_wolbachia_titer'] = np.log1p(titer)
    
    # Count cells with Wolbachia
    n_infected = np.sum(titer > 0)
    print(f"Detected Wolbachia in {n_infected} out of {adata.n_obs} cells ({n_infected/adata.n_obs*100:.2f}%)")
    
    return adata

def integrate_h5ad_files(directory_path, output_path, batch_key='batch', 
                         min_cells=3, min_genes=200, n_pcs=30, n_neighbors=15,
                         method='bbknn', calculate_titer=True):
    """
    Integrate multiple h5ad files with batch correction
    
    Parameters:
    -----------
    directory_path : str
        Path to directory containing h5ad files
    output_path : str
        Path to save integrated h5ad file
    batch_key : str
        Name of column to use for batch correction
    min_cells : int
        Minimum number of cells expressing a gene
    min_genes : int
        Minimum number of genes expressed in a cell
    n_pcs : int
        Number of principal components to use
    n_neighbors : int
        Number of neighbors for neighborhood graph
    method : str
        Batch correction method: 'bbknn', 'harmony', or 'both'
    calculate_titer : bool
        Whether to calculate Wolbachia titer before integration
    """
    # Set the figure directory to the directory containing the output file
    output_dir = os.path.dirname(output_path)
    figure_dir = os.path.join(output_dir, "figures")
    # Make sure the directory exists
    os.makedirs(figure_dir, exist_ok=True)
    # Set Scanpy figure directory
    sc.settings.figdir = figure_dir
    # Get all h5ad files in the directory
    
    h5ad_files = glob(os.path.join(directory_path, "*.h5ad"))
    
    if len(h5ad_files) == 0:
        print(f"No h5ad files found in {directory_path}")
        return None
    
    print(f"Found {len(h5ad_files)} h5ad files")
    
    # List to store individual datasets
    adatas = []
    
    # Load each dataset and add batch information
    for i, file_path in enumerate(h5ad_files):
        file_name = os.path.basename(file_path)
        batch_id = file_name.split('.')[0]  # Use filename without extension as batch ID
        print(f"Processing file {i+1}/{len(h5ad_files)}: {file_name}")
        
        # Load the dataset
        adata = sc.read_h5ad(file_path)
        adata.obs['Sample'] = os.path.splitext(os.path.basename(file_path))[0]
        
        # Check if var index contains tab characters and fix if needed
        if any('\t' in idx for idx in adata.var_names):
            print("Fixing var index with tab characters...")
            # Extract gene IDs and gene names from the index
            gene_ids = []
            gene_names = []
            
            for idx in adata.var_names:
                if '\t' in idx:
                    parts = idx.split('\t')
                    gene_id = parts[0]
                    gene_name = parts[1]
                    gene_ids.append(gene_id)
                    gene_names.append(gene_name)
                else:
                    # If there's no tab, use the same value for both id and name
                    gene_ids.append(idx)
                    gene_names.append(idx)
            
            # Create a new var DataFrame with gene names as index
            new_var = pd.DataFrame(index=gene_names)
            new_var['gene_ids'] = gene_ids
            new_var['feature_types'] = 'Gene Expression'
            
            # Copy over other columns from original var
            for col in adata.var.columns:
                new_var[col] = adata.var[col].values
            
            # Create a new AnnData object with fixed var
            # We need to get the X matrix and other components
            import anndata as ad
            new_adata = ad.AnnData(
                X=adata.X,
                obs=adata.obs,
                var=new_var,
                uns=adata.uns,
                obsm=adata.obsm if hasattr(adata, 'obsm') else None,
                varm=adata.varm if hasattr(adata, 'varm') else None,
                obsp=adata.obsp if hasattr(adata, 'obsp') else None,
                varp=adata.varp if hasattr(adata, 'varp') else None
            )
            adata = new_adata
        
        # Add batch information based on filename
        batch_id = file_name.split('.')[0]  # Use filename without extension as batch ID
        adata.obs[batch_key] = batch_id
        
        # Calculate Wolbachia titer if requested
        if calculate_titer:
            adata = calculate_wolbachia_titer(adata)
            
        adatas.append(adata)
        
    # Concatenate all datasets
    print("Concatenating datasets...")
    try:
        # Try using anndata.concat as recommended by the FutureWarning
        import anndata as ad
        combined = ad.concat(adatas, join='outer', merge='same', label=batch_key, index_unique='-')
    except:
        # Fall back to concatenate method if concat fails
        combined = adatas[0].concatenate(adatas[1:], join='outer', batch_key=batch_key)
    
    print(f"Combined data shape: {combined.shape}")
    
    # Basic preprocessing
    print("Performing basic preprocessing...")
    sc.pp.filter_cells(combined, min_genes=min_genes)
    sc.pp.filter_genes(combined, min_cells=min_cells)
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(combined, inplace=True)
    
    # Normalize the data
    print("Normalizing data...")
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)
    
    # Find highly variable genes
    print("Finding highly variable genes...")
    sc.pp.highly_variable_genes(combined, batch_key=batch_key)
    combined = combined[:, combined.var.highly_variable]
    
    # Run PCA
    print("Running PCA...")
    sc.pp.pca(combined, n_comps=n_pcs)
    
    # Save a copy of the unintegrated data for comparison
    combined_unintegrated = combined.copy()
    sc.pp.neighbors(combined_unintegrated, n_pcs=n_pcs)
    sc.tl.umap(combined_unintegrated)
    sc.pl.umap(combined_unintegrated, color=batch_key, save='_before_batch_correction.pdf')
    
    # Optional: Also plot Wolbachia titer on unintegrated data
    if 'wolbachia_titer' in combined_unintegrated.obs.columns:
        sc.pl.umap(combined_unintegrated, color='wolbachia_titer', save='_unintegrated_wolbachia_titer.pdf')
    
    # Run batch correction using the specified method
    if method == 'bbknn' or method == 'both':
        print("Performing batch correction with BBKNN...")
        try:
            import bbknn
            bbknn.bbknn(combined, batch_key=batch_key, n_pcs=n_pcs, neighbors_within_batch=5)
            
            # Save BBKNN-corrected object
            if method == 'both':
                combined_bbknn = combined.copy()
        except ImportError:
            print("BBKNN package not available. Please install it with: pip install bbknn")
            if method == 'bbknn':
                print("Falling back to standard PCA without batch correction.")
                sc.pp.neighbors(combined, n_pcs=n_pcs)
    
    if method == 'harmony' or method == 'both':
        print("Performing batch correction with Harmony...")
        # If we're using both methods, we need to start from scratch for harmony
        if method == 'both':
            combined = combined_unintegrated.copy()
        
        try:
            sc.external.pp.harmony_integrate(combined, batch_key)
            # Use harmony embedding for downstream analysis
            combined.obsm['X_pca_harmony'] = combined.obsm['X_pca_harmony']
            sc.pp.neighbors(combined, use_rep='X_pca_harmony')
        except ImportError:
            print("Harmony package not available. Please install it with: pip install harmonypy")
            if method == 'harmony':
                print("Falling back to standard PCA without batch correction.")
                sc.pp.neighbors(combined, n_pcs=n_pcs)
    
    # Run UMAP 
    print("Running UMAP and Leiden clustering...")
    sc.tl.umap(combined)
    sc.tl.leiden(combined, resolution=0.8)
    
    # Save the integrated object
    print(f"Saving integrated object to {output_path}")
    combined.write(output_path)
    
    # Generate diagnostic plots
    print("Generating diagnostic plots...")
    sc.pl.umap(combined, color=batch_key, save='_batch.pdf')
    sc.pl.umap(combined, color='leiden', save='_leiden.pdf')
    
    # If wolbachia_titer exists, plot it too
    if 'wolbachia_titer' in combined.obs.columns:
        sc.pl.umap(combined, color='wolbachia_titer', save='_wolbachia_titer.pdf')
        sc.pl.umap(combined, color='log1p_wolbachia_titer', save='_log1p_wolbachia_titer.pdf')
        
        # Create a violin plot of titer by batch
        sc.pl.violin(combined, 'wolbachia_titer', groupby=batch_key, save='_wolbachia_titer_by_batch.pdf')
        
        # Create a violin plot of titer by cluster
        sc.pl.violin(combined, 'wolbachia_titer', groupby='leiden', save='_wolbachia_titer_by_cluster.pdf')
    
    # If both methods are used, compare them
    if method == 'both':
        # Run UMAP on BBKNN object for comparison
        sc.tl.umap(combined_bbknn)
        sc.tl.leiden(combined_bbknn, resolution=0.8)
        
        # Generate comparison plots
        sc.pl.umap(combined_bbknn, color=batch_key, save='_bbknn_batch.pdf')
        sc.pl.umap(combined_bbknn, color='leiden', save='_bbknn_leiden.pdf')
        
        if 'wolbachia_titer' in combined_bbknn.obs.columns:
            sc.pl.umap(combined_bbknn, color='wolbachia_titer', save='_bbknn_wolbachia_titer.pdf')
        
        # Save the BBKNN object as well
        combined_bbknn.write(output_path.replace('.h5ad', '_bbknn.h5ad'))
        
        print("Saved both Harmony and BBKNN integrated objects for comparison")
    
    print("Integration complete!")
    return combined

def main():
    parser = argparse.ArgumentParser(description='Integrate multiple h5ad files with batch correction')
    
    # Required arguments
    parser.add_argument('--input_dir', type=str, required=True, 
                        help='Directory containing h5ad files to integrate')
    parser.add_argument('--output_file', type=str, required=True, 
                        help='Path to save the integrated h5ad file')
    
    # Optional arguments
    parser.add_argument('--batch_key', type=str, default='batch',
                        help='Name of column to use for batch correction')
    parser.add_argument('--min_cells', type=int, default=3,
                        help='Minimum number of cells expressing a gene')
    parser.add_argument('--min_genes', type=int, default=200,
                        help='Minimum number of genes expressed in a cell')
    parser.add_argument('--n_pcs', type=int, default=30,
                        help='Number of principal components to use')
    parser.add_argument('--n_neighbors', type=int, default=15,
                        help='Number of neighbors for neighborhood graph')
    parser.add_argument('--method', type=str, default='bbknn', choices=['bbknn', 'harmony', 'both'],
                        help='Batch correction method to use')
    parser.add_argument('--calculate_titer', action='store_true', 
                        help='Calculate Wolbachia titer for each cell before integration')
    
    args = parser.parse_args()
    
    # Run the integration
    integrated_adata = integrate_h5ad_files(
        directory_path=args.input_dir,
        output_path=args.output_file,
        batch_key=args.batch_key,
        min_cells=args.min_cells,
        min_genes=args.min_genes,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        method=args.method,
        calculate_titer=args.calculate_titer
    )
    
    # Print summary
    if integrated_adata is not None:
        print("Summary of integrated data:")
        print(f"Number of cells: {integrated_adata.n_obs}")
        print(f"Number of genes: {integrated_adata.n_vars}")
        print(f"Number of batches: {integrated_adata.obs[args.batch_key].nunique()}")
        print(f"Number of clusters: {integrated_adata.obs['leiden'].nunique()}")
        
        if 'wolbachia_titer' in integrated_adata.obs.columns:
            # Calculate percentage of infected cells
            n_infected = np.sum(integrated_adata.obs['wolbachia_titer'] > 0)
            print(f"Number of cells with Wolbachia: {n_infected} ({n_infected/integrated_adata.n_obs*100:.2f}%)")
            
            # Calculate average titer
            mean_titer = np.nanmean(integrated_adata.obs['wolbachia_titer'])
            median_titer = np.nanmedian(integrated_adata.obs['wolbachia_titer'])
            print(f"Average Wolbachia titer: mean={mean_titer:.4f}, median={median_titer:.4f}")
            
            # Calculate titer by batch
            for batch in integrated_adata.obs[args.batch_key].unique():
                batch_cells = integrated_adata[integrated_adata.obs[args.batch_key] == batch]
                n_batch_infected = np.sum(batch_cells.obs['wolbachia_titer'] > 0)
                mean_batch_titer = np.nanmean(batch_cells.obs['wolbachia_titer'])
                print(f"  {batch}: {n_batch_infected}/{batch_cells.n_obs} cells infected ({n_batch_infected/batch_cells.n_obs*100:.2f}%), mean titer={mean_batch_titer:.4f}")

if __name__ == "__main__":
    main()