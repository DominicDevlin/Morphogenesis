"""
UMAP on the Fashion MNIST Digits dataset using Datashader
---------------------------------------------------------

This is a simple example of using UMAP on the Fashion-MNIST
dataset. The goal of this example is largely to demonstrate
the use of datashader as an effective tool for visualising
UMAP results. In particular datashader allows visualisation
of very large datasets where overplotting can be a serious
problem. It supports coloring by categorical variables
(as shown in this example), or by continuous variables,
or by density (as is common in datashader examples).
"""
import umap
import numpy as np
import pandas as pd
import requests
import os
import datashader as ds
import datashader.utils as utils
import datashader.transfer_functions as tf
import matplotlib.pyplot as plt
import seaborn as sns
import pickle



obj1 = open('embedding', 'rb')
embedding = pickle.load(obj1)
obj1.close()

arrows = []
for i in range(len(embedding)):
  if (i < 20050 and i > 19900):
    new_arrow = []
    x1 = embedding[i-1][0]
    y1=embedding[i-1][1]
    x=embedding[i][0] - embedding[i-1][0]
    y=embedding[i][1] - embedding[i-1][1]
    length = np.sqrt(pow(x,2) + pow(y,2))
    new_arrow.append(x1)
    new_arrow.append(y1)
    new_arrow.append(x)
    new_arrow.append(y)
    new_arrow.append(length)
    arrows.append(new_arrow)
    

fig, ax = plt.subplots(figsize=(12, 10))
plt.scatter(embedding[:, 0], embedding[:, 1], s=0.3)#, cmap="Spectral", s=0.1)


for i in arrows:
  plt.arrow(i[0], i[1], i[2], i[3], head_width=np.log10(i[4])/3, width=(np.log10(i[4])/6), color="black")

plt.show()

# df = pd.DataFrame(embedding, columns=("x", "y"))

# df["class"] = pd.Series([str(x) for x in target], dtype="category")


# cvs = ds.Canvas(plot_width=400, plot_height=400)
# agg = cvs.points(df, "x", "y", ds.count_cat("class"))
# img = tf.shade(agg, min_alpha=255, color_key=color_key, how="eq_hist")

# utils.export_image(img, filename="tester", background="white", fmt=".png")

# image = plt.imread("tester.png")
# fig, ax = plt.subplots(figsize=(6, 6))
# plt.imshow(image)
# plt.setp(ax, xticks=[], yticks=[])
# plt.title(
#     "Cell expression data\n"
#     "into two dimensions by UMAP\n"
#     "visualised with Datashader",
#     fontsize=12,
# )

# plt.show()



# default
# pal = [
#     "#9e0142",
#     "#d8434e",
#     "#f67a49",
#     "#fdbf6f",
#     "#feeda1",
#     "#f1f9a9",
#     "#bfe5a0",
#     "#74c7a5",
#     "#378ebb",
#     "#5e4fa2",
#     "#17becf",
#     "#ff7f0e",
#     "#2ca02c",
#     "#d62728",
#     "#9467bd",
#     "#8c564b",
#     "#F2A2E8",
#     "#F8F6F0",
#     "#C9BE62",
#     "#FDBD01",
#     "#FFE5B4",
#     "#16E2F5",
#     "#A2AD9C",
#     "#00FA9A",
#     "#DAEE01",
#     "#B5EAAA",
#     "#C3FDB8",
#     "#F5FFFA",
#     "#F0FFFF",
#     "#368BC1", 
#     "#2F539B",
#     "#B6B6B4",
#     "#0C090A",
#     "#FF6700",
#     "#FF6700",
#     "#FF6700"
# ]


# this is for pluri56
# pal = [
# "#C2D354",
# "#8F2CFE",
# "#00FEFE",
# "#E4A600",
# "#BA6581"
# ]

# mushroom palette
# pal = [
# "#262626",
# "#E99A28",
# "#485AF2",
# "#8D378D",
# "#78EB4A",
# "#CE72ED",
# "#BDBDBD",
# "#C0C0C0",
# "#74c7a5",
# "#00FEFE",
# "#5e4fa2",
# "#F0D67C",
# "#E14600",
# "#A3ACFF",
# "#d62728",
# "#9467bd",
# "#8c564b",
# "#FE00FE",
# "#70C570",
# "#C9BE62",
# "#FDBD01",
# "#FF6700",
# "#16E2F5",
# "#A2AD9C",
# "#00FA9A",
# "#DAEE01",
# "#B5EAAA",
# "#C3FDB8",
# "#808000",
# "#808000",
# "#368BC1",
# "#2F539B",
# "#B6B6B4",
# "#0C090A"
# ]