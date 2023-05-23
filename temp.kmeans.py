# Kmeans clustering for the data

#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn.metrics import pairwise_distances
from sklearn.metrics import silhouette_score

data = pd.read_csv("yeast_raw_feature.tsv", sep="\t", header=None)

# %%
# Randomly select 100000 data points
new_data = data.sample(n=100000, random_state=1)

#%%
# Kmeans clustering with k = 1024
kmeans = KMeans(n_clusters=1024, random_state=0).fit(new_data)
labels = kmeans.labels_
print(labels)

# %%
# Plot the clustering result
# Magma color map
plt.scatter(new_data.iloc[:, 0], new_data.iloc[:, 1], c=labels, s=50, cmap='magma')
