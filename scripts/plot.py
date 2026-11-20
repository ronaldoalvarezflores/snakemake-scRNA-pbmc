import sys
sys.stdout = sys.stderr = open(snakemake.log[0], "w")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

adata = sc.read_h5ad(snakemake.input[0])
sample = snakemake.wildcards.sample

# Figure 1: UMAP coloured by Leiden cluster and by annotated cell type
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
sc.pl.umap(adata, color="leiden", legend_loc="on data", ax=axes[0], show=False,
           title=f"{sample}: Leiden clusters")
sc.pl.umap(adata, color="cell_type", ax=axes[1], show=False,
           title=f"{sample}: cell types")
plt.savefig(snakemake.output.umap, dpi=150, bbox_inches="tight")
plt.close()

# Figure 2: marker gene expression per cell type
markers = {k: [g for g in v if g in adata.raw.var_names]
           for k, v in snakemake.params.markers.items()}
sc.pl.dotplot(adata, markers, groupby="cell_type", show=False)
plt.savefig(snakemake.output.markers, dpi=150, bbox_inches="tight")
