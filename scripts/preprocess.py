import sys
sys.stdout = sys.stderr = open(snakemake.log[0], "w")

import scanpy as sc

p = snakemake.params

adata = sc.read_10x_mtx(snakemake.input[0], var_names="gene_symbols")
adata.var_names_make_unique()
print(f"Raw: {adata.n_obs} cells x {adata.n_vars} genes")

# Filter low-quality cells and rare genes
sc.pp.filter_cells(adata, min_genes=p.min_genes)
sc.pp.filter_genes(adata, min_cells=p.min_cells)

# Remove dying cells (high % of mitochondrial genes)
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs["pct_counts_mt"] < p.max_pct_mt].copy()
print(f"After filtering: {adata.n_obs} cells x {adata.n_vars} genes")

# Normalize to 10,000 counts per cell and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.write(snakemake.output[0])
