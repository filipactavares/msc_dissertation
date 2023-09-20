# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 14:31:11 2023

@author: Utilizador
"""

import os
os.chdir("C:\\Users\\Utilizador\\OneDrive\\Documentos\\tese\\UMaps")
os.getcwd()


#Install, import, configure & load data into a AnnData object#


#!pip install python-louvain
#!pip install scikit-learn
#!pip install scanpy
#!pip install python-igraph
#!pip install louvain
#!pip install pandas
#!pip install umap-learn
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as plt
import umap
import seaborn as sns

# Set verbosity level of scanpy to 3. Level 3 means that scanpy will also show 'hint' messages.
sc.settings.verbosity = 3

# Set some figure parameters for scanpy. Specifically, set the resolution of rendered figures to 150 DPI and the face color of figures to white.
sc.settings.set_figure_params(dpi=150, facecolor="white")

# Print a header with information about scanpy, such as its version number.
sc.logging.print_header()

#import data
df = pd.read_csv('all_just_mutations.csv', index_col=0)
df_all = pd.read_csv('all.csv',index_col=0)
df_all = df_all.dropna()

adata = sc.AnnData(df.astype('float32'))
adata

adata_all = sc.AnnData(df_all.astype('float32'))
adata_all


#UMAP/Louvain using scanpy

#sc.tl.pca(adata, svd_solver='arpack')

#sc.pl.pca_variance_ratio(adata, log=False, n_pcs=50)

#This code is calling the neighbors function from scanpy.pp module. 
#The neighbors function computes a neighborhood graph of observations. 
#In this specific example, it is being called with three arguments: adata, which represents the input data; n_neighbors=15, which specifies that 15 neighboring data points should be used for manifold approximation; and n_pcs=6, which specifies that 6 PCs should be used.
#The choice of n_neighbors depends on your data and analysis goals. 
#A larger value of n_neighbors will result in more global views of the manifold, while smaller values will preserve more local data 1. You can try different values of n_neighbors and see how it affects your downstream analysis.
#sc.pp.neighbors(adata, n_neighbors=500, n_pcs=50)
#Run this on 132 dimensions of `.X`, if you really want this, set `use_rep='X'`.

sc.pp.neighbors(adata, n_neighbors=600,use_rep='X') #n_neighbours=600


# Compute a UMAP (Uniform Manifold Approximation and Projection) embedding of the data stored in 'adata'.
# Set the following parameters for UMAP:
# - 'alpha' controls the level of exponential decay for the graph weights. Set it to 0.3.
# - 'min_dist' controls how tightly UMAP is allowed to pack points together. Set it to 0.6.
# - 'random_state' controls the seed used by the random number generator. Set it to 0.
sc.tl.umap(adata, alpha = 0.3, min_dist = 1, random_state = 0)


sc.tl.louvain(adata, resolution=0.75) #0.6

sc.pl.umap(adata, color=["louvain"], wspace=0.5, ncols=2, size=[50])



subtype=df_all['SUBTIPO']
adata.obs["subtype"] = pd.Categorical(subtype)  # Categoricals are preferred for efficiency

arv=df_all['ARV(2)']
adata.obs["ARV"] = pd.Categorical(arv)

treatment=df_all['TREATMENT']
adata.obs["treatment"] = pd.Categorical(treatment)

viral_load=df_all['RESULTADO_CV_N']
adata.obs["viral load"] = pd.Categorical(viral_load)

cd4 = df_all['RESULTADO_CD4_N']
adata.obs["CD4+"] = pd.Categorical(cd4)

sex=df_all['SEXO']
adata.obs["sex"] = pd.Categorical(sex)

age=df_all['IDADE_N']
adata.obs["age"] = pd.Categorical(age)

sc.pl.umap(adata, color=['subtype'], wspace=0.5, ncols=2, size=[50])
sc.pl.umap(adata, color=['ARV'], wspace=0.5, ncols=2, size=[50])
sc.pl.umap(adata, color=['treatment'], wspace=0.5, ncols=2, size=[50])
sc.pl.umap(adata, color=['sex'], wspace=0.5, ncols=2, size=[50])
sc.pl.umap(adata, color=['viral load'], wspace=0.5, ncols=2, size=[50])
sc.pl.umap(adata, color=['CD4+'], wspace=0.5, ncols=2, size=[50])
sc.pl.umap(adata, color=['age'], wspace=0.5, ncols=2, size=[50])


#rank mutations#

sc.tl.rank_genes_groups(adata, 'louvain', method = 'wilcoxon', n_genes=50)

adata.uns["rank_genes_groups"]

#(without threshold)
#sc.pl.rank_genes_groups_heatmap(adata, show_gene_labels=True)
#sc.pl.rank_genes_groups_dotplot(adata)
#mp=sc.pl.rank_genes_groups_matrixplot(adata,return_fig=True)
#mp.add_totals().style(edge_color='black').show()
#sc.pl.rank_genes_groups_tracksplot(adata)
#sc.pl.heatmap(adata, var_names=adata.var_names, groupby='louvain', show_gene_labels=True, cmap='viridis', dendrogram=True)


# Calculate mean expression value for each variable
mean_expression = adata.X.mean(axis=0)

# Set threshold for filtering low-value variables
threshold = 0.005

# Create list of variable names with mean expression above threshold
high_value_var_names = adata.var_names[mean_expression > threshold]

# Plot heatmap with filtered variable names
sc.tl.dendrogram(adata,n_pcs=5,groupby=['louvain'])
sc.pl.heatmap(adata, var_names=high_value_var_names, groupby='louvain', show_gene_labels=True, cmap='viridis', dendrogram=True)





