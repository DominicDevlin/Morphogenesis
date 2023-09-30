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
source_df = pd.read_csv("cdata.csv")
data = source_df.iloc[:, :27].values.astype(np.float32)
target = source_df["class"].values


color_key = {
'0' : "#00ff00",
'7' : "#00ff00",
'12' : "#7171c6",
'19' : "#a4adad",
'27' : "#e8e8e8",
'32' : "#5680ab",
'33' : "#e2d0d0",
'34' : "#664e4e",
'35' : "#d5d5d5",
'36' : "#999999",
'38' : "#d6d6d6",
'39' : "#adadad",
'41' : "#737373",
'42' : "#4c4c4c",
'44' : "#efefef",
'45' : "#8de38d",
'46' : "#9b9b9b",
'47' : "#86e3e3",
'48' : "#9c9c9c",
'49' : "#bfbfbf",
'50' : "#838383",
'51' : "#353535",
'52' : "#b00000",
'53' : "#212121",
'54' : "#085b08",
'55' : "#bbbbbb",
'59' : "#161616",
'62' : "#929292",
'63' : "#7d7d7d",
'64' : "#0000b0",
'65' : "#005b00",
'66' : "#6920ac",
'278' : "#45b9bc",
'279' : "#c2b23b",
'280' : "#70c570",
'281' : "#c4639c",
'282' : "#ff9800",
'283' : "#727272",
'284' : "#48b5ff",
}

colours = []
for i in target:
  val = str(i)
  colours.append(color_key[val])




obj1 = open('embedding', 'rb')
embedding = pickle.load(obj1)
obj1.close()


df = pd.DataFrame(embedding, columns=("x", "y"))

df["class"] = pd.Series([str(x) for x in target], dtype="category")


cvs = ds.Canvas(plot_width=400, plot_height=400)
agg = cvs.points(df, "x", "y", ds.count_cat("class"))
img = tf.shade(agg, min_alpha=255, color_key=color_key, how="eq_hist")

utils.export_image(img, filename="tester", background="white", fmt=".png")

image = plt.imread("tester.png")
fig, ax = plt.subplots(figsize=(6, 6))
plt.imshow(image)
plt.setp(ax, xticks=[], yticks=[])
plt.title(
    "Cell expression data\n"
    "into two dimensions by UMAP\n"
    "visualised with Datashader",
    fontsize=12,
)

plt.show()



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