# import umap
import numpy as np
import pandas as pd
# import requests
# import os
# import datashader as ds
# import datashader.utils as utils
# import datashader.transfer_functions as tf
# import matplotlib.pyplot as plt
# import seaborn as sns
import scvelo as scv
import anndata as ad
from scipy.sparse import csr_matrix

scv.set_figure_params()

# ndata = scv.datasets.pancreas()
# print(ndata)

frame = ad.read_csv("layer1.csv")
frame2 = ad.read_csv("layer2.csv")
types = ad.read_text("colours.txt")
embedding = ad.read_csv("embedding.csv")

adata = ad.AnnData(frame)
adata.layers["spliced"] = ad.AnnData(frame)
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
adata.layers["unspliced"] = ad.AnnData(frame2)

# adata.layers["spliced"] = ad.AnnData(frame)
adata.obsm["X_umap"] = embedding

adata.obs["clusters"] = types


# print(adata.obs_names[:10])
# print(adata.obs)
print(adata)






# counts = csr_matrix(np.random.poisson(1, size=(100, 2000)), dtype=np.float32)
# adata = ad.AnnData(counts)



# scv.pp.filter_and_normalize(adata, min_shared_counts=20) # , min_shared_cells=20, n_top_genes=2000)
# scv.pp.moments(adata) #, n_pcs=30, n_neighbors=30)

scv.tl.velocity(adata)

scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding_stream(adata, basis='umap')



