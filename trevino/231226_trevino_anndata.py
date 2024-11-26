# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "anndata==0.10.9",
#     "numpy==2.0.2",
#     "pandas==2.2.3",
#     "pybiomart==0.2.0",
#     "scipy==1.13.1",
#     "scanpy[leiden]==1.10.3",
#     "seaborn==0.13.2",
#     "scvi-tools",
#     "torch==2.5.1",
# ]
# ///

import marimo

__generated_with = "0.9.20"
app = marimo.App()


@app.cell
def __():
    import os
    import numpy as np
    import pandas as pd
    import scipy
    import anndata
    import scanpy as sc
    import pybiomart
    import scvi
    import torch
    import random
    import seaborn as sns
    return (
        anndata,
        np,
        os,
        pd,
        pybiomart,
        random,
        sc,
        scipy,
        scvi,
        sns,
        torch,
    )


@app.cell
def __(os):
    cwd = os.getcwd()
    cwd
    return (cwd,)


@app.cell
def __(cwd, pd):
    meta = pd.read_csv(cwd+'/data/GSE162170_rna_cell_metadata.txt.gz', compression='gzip', sep='\t')
    meta.index = meta['Cell.ID']
    meta.index.name = None
    meta
    return (meta,)


@app.cell
def __(cwd, pd):
    # magic command not supported in marimo; please file an issue to add support
    # %%time
    counts = pd.read_csv(cwd+'/data/GSE162170_rna_counts.tsv.gz', compression='gzip', sep='\t')
    counts = counts.transpose()
    counts
    return (counts,)


@app.cell
def __(anndata, counts, meta):
    # magic command not supported in marimo; please file an issue to add support
    # %%time
    adata = anndata.AnnData(X=counts,
                            obs=meta,
                            var=counts.columns.to_frame())
    adata
    return (adata,)


@app.cell
def __(a, adata, sc):
    _a = sc.queries.biomart_annotations('hsapiens', ['ensembl_gene_id', 'hgnc_symbol'])
    b = dict(zip(a['ensembl_gene_id'], a['hgnc_symbol']))
    adata.var['hgnc_symbol'] = adata.var[0].map(b)
    adata.var.index = adata.var['hgnc_symbol']
    adata.var.index.name = None
    adata.var.drop(columns=[0, 'hgnc_symbol'], inplace=True)
    adata.var
    return (b,)


@app.cell
def __(adata):
    adata.var_names_make_unique
    return


@app.cell
def __(a, adata):
    _a = ~adata.var.index.isnull()
    adata_1 = adata[:, a].copy()
    return (adata_1,)


@app.cell
def __(adata_1):
    adata_1
    return


@app.cell
def __(adata_1, sc):
    sc.pp.filter_genes(adata_1, min_cells=10)
    sc.pp.filter_cells(adata_1, min_genes=200)
    return


@app.cell
def __():
    #sc.pp.subsample(adata, n_obs=40000, random_state=0, copy=False)
    return


@app.cell
def __(adata_1, scipy):
    adata_1.X = scipy.sparse.csr_matrix(adata_1.X.copy())
    adata_1.layers['counts'] = scipy.sparse.csr_matrix(adata_1.X.copy())
    return


@app.cell
def __(adata_1, random, scvi):
    random.seed(17)
    scvi.model.SCVI.setup_anndata(adata_1, layer='counts', batch_key='Sample.ID')
    scvi_model = scvi.model.SCVI(adata_1, n_layers=2, n_latent=30, n_hidden=128, gene_likelihood='nb')
    scvi_model.train()
    return (scvi_model,)


@app.cell
def __(adata_1, random, scvi_model):
    random.seed(17)
    adata_1.obsm['X_scvi'] = scvi_model.get_latent_representation()
    adata_1.layers['counts_scvi'] = scvi_model.get_normalized_expression(library_size=10000)
    return


@app.cell
def __(adata_1, sc):
    sc.pp.neighbors(adata_1, use_rep='X_scvi', key_added='neighbors_scvi', n_neighbors=20)
    sc.tl.leiden(adata_1, neighbors_key='neighbors_scvi', key_added='leiden_scvi', resolution=3)
    sc.tl.umap(adata_1, neighbors_key='neighbors_scvi')
    sc.pl.umap(adata_1, color=['leiden_scvi'], legend_loc='on data')
    return


@app.cell
def __(adata_1, sc):
    sc.pp.neighbors(adata_1, use_rep='X_scvi', key_added='neighbors_scvi', n_neighbors=20)
    sc.tl.leiden(adata_1, neighbors_key='neighbors_scvi', key_added='leiden_scvi', resolution=3)
    sc.tl.umap(adata_1, neighbors_key='neighbors_scvi')
    sc.pl.umap(adata_1, color=['leiden_scvi'], legend_loc='on data')
    return


@app.cell
def __(adata_1, cwd):
    adata_1.write(cwd + '/outs/231226_trevino_rna_scvi.h5ad')
    return


if __name__ == "__main__":
    app.run()

