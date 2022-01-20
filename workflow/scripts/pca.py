#!/usr/bin/env python
# coding: utf-8

"""
pca 
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

from tools import loadZarrArrays, getCohorts, run_pca, plot_coords, hash_params
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.express as px
import dask.array as da
import seaborn as sns
import matplotlib.pyplot as plt


# Garuds Selection Scans # 
chrom = snakemake.wildcards['chrom']
dataset = snakemake.params['dataset']
genotypePath = snakemake.input['genotypes']
positionsPath = snakemake.input['positions']
siteFilterPath = snakemake.input['siteFilters']

results_dir = snakemake.params['data']

# Read metadata 
metadata = pd.read_csv(snakemake.params['metadata'], sep=",")
metadata['location'] = metadata['location'].str.split(".").str.get(0)

# Load Arrays
snps, pos = loadZarrArrays(genotypePath, positionsPath, siteFilterPath=siteFilterPath, haplotypes=False)

# Determine cohorts
cohorts = getCohorts(metadata, columns=['species_gambiae_coluzzii'])


# choose colours for species
species_palette = px.colors.qualitative.Plotly
species_color_map = {
    'gambiae': species_palette[0],
    'coluzzii': species_palette[1],
    'arabiensis': species_palette[2],
    'intermediate_gambiae_coluzzii': species_palette[3],
    'intermediate_arabiensis_gambiae': species_palette[4],
}


# Run PCA on whole dataset together
data, evr = run_pca(contig=chrom, gt=snps, pos=pos, df_samples=metadata,
    sample_sets=dataset, results_dir=results_dir
)
evr = evr.astype("float").round(4) # round decimals for variance explained % 

plot_coords(data, evr, title=f" PCA | {dataset} | {chrom}", filename=f"results/PCA/{dataset}.{chrom}.html")

fig = plt.figure(figsize=(10, 10))
fig = sns.scatterplot('PC1','PC2', data=data, hue='location')
fig.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title(f"PCA | {dataset} | {chrom}", fontsize=14)
plt.xlabel(f"PC1 ({evr[0]*100} % variance explained)", fontdict={"fontsize":14})
plt.ylabel(f"PC2 ({evr[1]*100} % variance explained)", fontdict={"fontsize":14})
plt.savefig(f"results/PCA/{dataset}.{chrom}.png")



# Loop through each cohort, manipulate genotype arrays and calculate chosen Garuds Statistic
for idx, cohort in cohorts.iterrows():

    # filter to correct loc, year, species individuals
    gt_cohort = snps.take(cohort['indices'], axis=1)
    meta = metadata.take(cohort['indices'])
    
    
    data, evr = run_pca(contig=chrom, gt=gt_cohort, pos=pos, df_samples=meta,
        sample_sets=cohort['cohortNoSpaceText'], results_dir=results_dir
    )
    evr = evr.astype("float").round(4)

    plot_coords(data, evr, title=f" PCA | {cohort['cohortText']} | {chrom}", filename=f"results/PCA/{cohort['cohortNoSpaceText']}.{chrom}.html")

    fig = plt.figure(figsize=(10, 10))
    fig = sns.scatterplot('PC1','PC2', data=data, hue='location')
    fig.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title(f"PCA | {cohort['cohortText']} | {chrom}", fontsize=14)
    plt.xlabel(f"PC1 ({evr[0]*100} % variance explained)", fontdict={"fontsize":14})
    plt.ylabel(f"PC2 ({evr[1]*100} % variance explained)", fontdict={"fontsize":14})
    plt.savefig(f"results/PCA/{cohort['cohortNoSpaceText']}.{chrom}.png")
