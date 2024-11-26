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
def __(sc, sns):
    def PlotUMAP(adata, markers, layer='log2_counts_scvi', size=2, vmin='p0', vmax='p99'):
        for i in range(len(markers)):
            sc.pl.umap(adata,
                      color=markers[i],
                      layer=layer,
                      size=size,
                      cmap=sns.blend_palette(['lightgray', sns.xkcd_rgb['red orange']], as_cmap=True),
                       vmin=vmin, vmax=vmax)
    return (PlotUMAP,)


@app.cell
def __(os, sc):
    cwd = os.getcwd()
    adata = sc.read_h5ad(cwd+'/outs/231226_trevino_rna_scvi.h5ad')
    adata
    return adata, cwd


@app.cell
def __(adata, sc):
    adata.layers['log2_counts'] = sc.pp.log1p(adata.layers['counts'].copy(), base=2)
    adata.layers['log2_counts_scvi'] = sc.pp.log1p(adata.layers['counts_scvi'].copy(), base=2)
    return


@app.cell
def __(adata, sc):
    sc.pl.umap(adata, color=['leiden_scvi'], legend_loc='on data')
    return


@app.cell
def __(adata):
    adata.obs
    return


@app.cell
def __(adata, sc):
    sc.pl.umap(adata, color=['Age'])
    sc.pl.umap(adata, color=['Sample.ID'])
    sc.pl.umap(adata, color=['Tissue.ID'])
    sc.pl.umap(adata, color=['DF_classification'])
    return


@app.cell
def __(adata):
    adata.var_names_make_unique()
    return


@app.cell
def __(PlotUMAP, adata, sc):
    _markers = ['TOP2A', 'MKI67', 'SOX9', 'HES1', 'FBXO32', 'CCN2', 'MOXD1', 'HOPX', 'CRYAB', 'NR4A1', 'FOXJ1']
    sc.pl.dotplot(adata, _markers, groupby='leiden_scvi')
    PlotUMAP(adata, _markers, layer='log2_counts_scvi', size=5)
    return


@app.cell
def __(PlotUMAP, adata, sc):
    _markers = ['ASCL1', 'OLIG2', 'PDGFRA', 'EGFR', 'AQP4', 'APOE', 'SOX10', 'NKX2-2', 'MBP']
    sc.pl.dotplot(adata, _markers, groupby='leiden_scvi')
    PlotUMAP(adata, _markers, layer='log2_counts_scvi', size=5)
    return


@app.cell
def __(PlotUMAP, adata, sc):
    _markers = ['EOMES', 'PPP1R17', 'NEUROG1', 'BCL11B', 'SATB2', 'SLC17A7', 'NR4A2', 'CRYM']
    sc.pl.dotplot(adata, _markers, groupby='leiden_scvi')
    PlotUMAP(adata, _markers, layer='log2_counts_scvi', size=5)
    return


@app.cell
def __(PlotUMAP, adata, sc):
    _markers = ['DLX2', 'GAD2', 'LHX6', 'SST', 'SP8', 'NR2F2', 'MEIS2', 'ETV1']
    sc.pl.dotplot(adata, _markers, groupby='leiden_scvi')
    PlotUMAP(adata, _markers, layer='log2_counts_scvi', size=5)
    return


@app.cell
def __(PlotUMAP, adata, sc):
    _markers = ['AIF1', 'CCL3', 'CLDN5', 'PECAM1', 'FOXC2', 'PDGFRB', 'HOPX', 'COL1A1', 'LUM', 'CRYAB', 'EGFR', 'HEMGN']
    sc.pl.dotplot(adata, _markers, groupby='leiden_scvi')
    PlotUMAP(adata, _markers, layer='log2_counts_scvi', size=5)
    return


@app.cell
def __(adata, sc):
    cluster_labels = {'0':'IN MGE',
                      '1':'IN CGE',
                      '2':'GluN CPN',
                      '3':'GluN SCPN',
                      '4':'GluN SCPN',
                      '5':'GluN SCPN',
                      '6':'GluN SCPN',
                      '7':'GluN SCPN',
                      '8':'GluN SCPN',
                      '9':'GluN CPN/SCPN',
                      '10':'vRG',
                      '11':'GluN SCPN',
                      '12':'IPC',
                      '13':'mGPC',
                      '14':'GluN SCPN',
                      '15':'GluN SCPN',
                      '16':'GluN SCPN',
                      '17':'IN CGE',
                      '18':'oRG',
                      '19':'IPC',
                      '20':'GluN CPN',
                      '21':'GluN CPN',
                      '22':'oRG',
                      '23':'GluN SCPN',
                      '24':'vRG/oRG',
                      '25':'oRG',
                      '26':'OPC/Oligo',
                      '27':'GluN CPN',
                      '28':'SP',
                      '29':'GluN SCPN',
                      '30':'oRG',
                      '31':'GluN CPN',
                      '32':'GluN CPN',
                      '33':'IN CGE/MGE',
                      '34':'GluN SCPN',
                      '35':'tRG',
                      '36':'IN CGE/MGE',
                      '37':'MG',
                      '38':'MG',
                      '39':'Pericyte',
                      '40':'GluN SCPN',
                      '41':'vRG',
                      '42':'GluN CPN',
                      '43':'IN MGE',
                      '44':'EC',
                      '45':'vRG/oRG',
                      '46':'GluN CPN/SCPN',
                      '47':'GluN CPN/SCPN',
                      '48':'IN CGE',
                      '49':'IN MGE',
                      '50':'VLMC',
                      '51':'Unknown',
                      '52':'RBC'}
    adata.obs['CellType'] = adata.obs['leiden_scvi'].map(cluster_labels)
    sc.pl.umap(adata, color=['CellType'], size=5, legend_loc='on data', legend_fontsize=7)
    return (cluster_labels,)


@app.cell
def __(adata, cwd):
    # save output
    adata.write(cwd+'/outs/231227_trevino_rna_scvi.h5ad')
    return


@app.cell
def __(adata, sc):
    _markers = ['FOXG1', 'CRYAB', 'AQP4', 'EGFR', 'OLIG2', 'SOX10', 'HES5', 'HOPX', 'MOXD1', 'TOP2A', 'ASCL1', 'EOMES', 'NEUROG2', 'NEUROD2', 'BCL11B', 'SATB2', 'DLX2', 'LHX6', 'SP8', 'CLDN5', 'PECAM1', 'FOXC2', 'AIF1']
    sc.pl.dotplot(adata, _markers, groupby='CellType')
    return


if __name__ == "__main__":
    app.run()

