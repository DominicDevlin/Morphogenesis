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
import scvelo as scv

sns.set(context="paper", style="white")

# if not os.path.isfile("fashion-mnist.csv"):
#     csv_data = requests.get("https://www.openml.org/data/get_csv/18238735/phpnBqZGZ")
#     with open("fashion-mnist.csv", "w") as f:
#         f.write(csv_data.text)
# source_df = pd.read_csv("fashion-mnist.csv")

source_df = pd.read_csv("cdata.csv")


data = source_df.iloc[:, :27].values.astype(np.float32)
target = source_df["class"].values


max_val=0
for i in target:
  if i > max_val:
    max_val = i

# color_key = {str(d): c for d, c in enumerate(pal)}


# count = 0
# while max_val+1 > len(color_key):
#   color_key[str(len(color_key))] = pal[count]
#   count += 1
#   if count >= len(pal): 
#     count = 0


color_key = {
'4' : "#0000fe",
'5' : "#ff00ff",
'6' : "#00ffff",
'7' : "#00ff00",
'8' : "#555555",
'9' : "#c67171",
'10' : "#71c671",
'11' : "#8e8e38",
'12' : "#7171c6",
'13' : "#8e388e",
'16' : "#4c719e",
'17' : "#495b5b",
'18' : "#232b2b",
'19' : "#a4adad",
'20' : "#00375b",
'21' : "#7391a5",
'22' : "#ffbf00",
'23' : "#001828",
'24' : "#99913a",
'25' : "#c1c1c1",
'26' : "#272727",
'27' : "#e8e8e8",
'28' : "#707070",
'29' : "#bebebe",
'30' : "#b98e8e",
'31' : "#729c9c",
'32' : "#5680ab",
'33' : "#e2d0d0",
'34' : "#664e4e",
'35' : "#d5d5d5",
'36' : "#999999",
'37' : "#00003f",
'38' : "#d6d6d6",
}

colours = []
for i in target:
  val = str(i)
  colours.append(color_key[val])




reducer = umap.UMAP(random_state=42, n_neighbors=50, min_dist=0.5)
embedding = reducer.fit_transform(data)


arrows = []

for i in range(len(embedding)):
  while (i < 1000 and i > 0):
    new_arrow = []
    x1=embedding[i-1][0]
    y1=embedding[i-1][1]
    x2=embedding[i][0]
    y2=embedding[i][1]
    new_arrow.append(x1)
    new_arrow.append(y1)
    new_arrow.append(x2-x1)
    new_arrow.append(y2-y1)
    arrows.append(new_arrow)
    

print(arrows)

fig, ax = plt.subplots(figsize=(12, 10))
plt.scatter(embedding[:, 0], embedding[:, 1], c=colours, s=0.3)#, cmap="Spectral", s=0.1)

count = 0
for i in range(len(embedding)):
  if (count % 100 == 0):
    embedding[:,]

for i in arrows:
  plt.arrow(i[0], i[1], i[2], i[3], head_width=0.5, width=0.1, color="black")

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