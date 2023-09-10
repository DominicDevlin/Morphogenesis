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

sns.set(context="paper", style="white")

# if not os.path.isfile("fashion-mnist.csv"):
#     csv_data = requests.get("https://www.openml.org/data/get_csv/18238735/phpnBqZGZ")
#     with open("fashion-mnist.csv", "w") as f:
#         f.write(csv_data.text)
# source_df = pd.read_csv("fashion-mnist.csv")

source_df = pd.read_csv("cdata.csv")


data = source_df.iloc[:, :27].values.astype(np.float32)
target = source_df["class"].values

pal = [
    "#9e0142",
    "#d8434e",
    "#f67a49",
    "#fdbf6f",
    "#feeda1",
    "#f1f9a9",
    "#bfe5a0",
    "#74c7a5",
    "#378ebb",
    "#5e4fa2",
    "#17becf",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#9e0142",
    "#d8434e",
    "#f67a49",
    "#fdbf6f",
    "#feeda1",
    "#f1f9a9",
    "#bfe5a0",
    "#74c7a5",
    "#378ebb",
    "#5e4fa2",
    "#17becf",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b"
]
color_key = {str(d): c for d, c in enumerate(pal)}

reducer = umap.UMAP(random_state=42, n_neighbors=30, min_dist=0.8)
embedding = reducer.fit_transform(data)


# fig, ax = plt.subplots(figsize=(12, 10))
# # color = mnist.target.astype(int)
# plt.scatter(embedding[:, 0], embedding[:, 1], c=source_df["class"].values, cmap="Spectral", s=0.1)
# plt.show()


df = pd.DataFrame(embedding, columns=("x", "y"))
df["class"] = pd.Series([str(x) for x in target], dtype="category")

cvs = ds.Canvas(plot_width=400, plot_height=400)
agg = cvs.points(df, "x", "y", ds.count_cat("class"))
img = tf.shade(agg, color_key=color_key, min_alpha=255, how="eq_hist")

utils.export_image(img, filename="tester", background="black")

image = plt.imread("tester.png")
fig, ax = plt.subplots(figsize=(6, 6))
plt.imshow(image)
plt.setp(ax, xticks=[], yticks=[])
plt.title(
    "Fashion MNIST data embedded\n"
    "into two dimensions by UMAP\n"
    "visualised with Datashader",
    fontsize=12,
)

plt.show()
