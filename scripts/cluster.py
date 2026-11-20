import sys
sys.stdout = sys.stderr = open(snakemake.log[0], "w")

import scanpy as sc

p = snakemake.params
sc.settings.n_jobs = snakemake.threads

adata = sc.read_h5ad(snakemake.input[0])

# Keep all genes in .raw (for marker genes), analyse only highly variable ones
sc.pp.highly_variable_genes(adata, n_top_genes=p.n_top_genes)
adata.raw = adata
adata = adata[:, adata.var["highly_variable"]].copy()

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=p.n_neighbors, n_pcs=p.n_pcs)
sc.tl.leiden(adata, resolution=p.resolution, flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(adata)
print(f"{adata.obs['leiden'].nunique()} Leiden clusters")

# Annotation: score each cell for every marker set, then give each
# cluster the cell type with the highest mean score
for cell_type, genes in p.markers.items():
    genes = [g for g in genes if g in adata.raw.var_names]
    sc.tl.score_genes(adata, genes, score_name=cell_type, use_raw=True)

scores = adata.obs.groupby("leiden", observed=True)[list(p.markers)].mean()
labels = scores.idxmax(axis=1)
print("Cluster annotation:\n" + labels.to_string())
adata.obs["cell_type"] = adata.obs["leiden"].map(labels).astype("category")

adata.write(snakemake.output[0])
