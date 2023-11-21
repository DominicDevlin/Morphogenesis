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
'0'  : "#4c719e",
'16' : "#4c719e",
'24' : "#99913a",
'32' : "#5680ab",
'34' : "#664e4e",
'35' : "#d5d5d5",
'37' : "#87459b",
'38' : "#d6d6d6",
'47' : "#86e3e3",
'59' : "#161616",
'60' : "#2d2d2d",
'62' : "#929292",
'67' : "#b2c0dc",
'73' : "#ff3883",
'76' : "#a8610c",
'78' : "#6dd339",
'79' : "#ab4d00",
'88' : "#c665ff",
'91' : "#2e8c68",
'92' : "#35461e",
'95' : "#a200ff",
'99' : "#b1b115",
'102' : "#9b2eff",
'105' : "#ffa623",
'108' : "#adecc0",
'109' : "#28a7a7",
'110' : "#66c88f",
'112' : "#43ff3f",
'113' : "#646c05",
'114' : "#63a5b3",
'116' : "#191a3d",
'117' : "#756bd0",
'118' : "#d3ffd8",
'119' : "#52dec5",
'120' : "#cccc33",
'121' : "#ff3333",
'122' : "#153dc1",
'124' : "#1b576f",
'125' : "#505050",
'126' : "#5d58cb",
'127' : "#d3aad9",
'128' : "#72e558",
'129' : "#93b5ff",
'130' : "#2c59a2",
'131' : "#58252f",
'132' : "#274373",
'133' : "#da9940",
'135' : "#2ba540",
'136' : "#2c45d8",
'137' : "#b5ff61",
}

colours = []
for i in target:
  val = str(i)
  colours.append(color_key[val])



reducer = umap.UMAP(random_state=42, n_neighbors=50, min_dist=0.5)
embedding = reducer.fit_transform(data)


with open('embedding', 'wb') as f:
  pickle.dump(embedding, f)

# print(arrows)

fig, ax = plt.subplots(figsize=(12, 10))
plt.scatter(embedding[:, 0], embedding[:, 1], c=colours, s=0.3)#, cmap="Spectral", s=0.1)

count = 0
for i in range(len(embedding)):
  if (count % 100 == 0):
    embedding[:,]


# with open('my_plot.pkl', 'wb') as f:
#     pickle.dump(fig, f)

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