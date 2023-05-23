# PCA
#%%
import pandas
import numpy
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import scipy.cluster.hierarchy as hcluster

#%%
# Load data
data = pandas.read_csv(
    '/Users/hbk/Projects/Lab/06_FoldMotif/repos/motifsearch/data/temp.vec.data.raw.tsv', sep='\t', header=0)

data_std = StandardScaler().fit_transform(data)

#%%
threshold = 5
clusters = hcluster.fclusterdata(data, threshold, criterion="distance")

#%%

plt.scatter(data_std[:,0], data_std[:,1], c=clusters, cmap='prism')

#%%
# Cluster number
n_clusters = len(set(clusters))
print(n_clusters)


#%%
# PCA
pca = PCA(n_components=2)
pca.fit(data)
data_pca = pca.transform(data)
#%%

# Plot
# plt.scatter(data_pca[:,0], data_pca[:,1], alpha=0.02)

plt.scatter(data.iloc[:,0], data.iloc[:,1], c=clusters, alpha=0.02)
plt.show()


# %%
# Plot PCA with cluster
plt.scatter(data_pca[:,0], data_pca[:,1], c=clusters, cmap='prism', alpha=0.02)
# %%
# Print hierarhical cluster results
# print(clusters)
# Print centroid of each hierarchical cluster
print(hcluster)

#%%
# DBSCAN
from sklearn.cluster import DBSCAN

#%%
# DBSCAN
dbscan = DBSCAN(eps=0.2, min_samples=2)
dbscan.fit(data_std)
data_dbscan = dbscan.fit_predict(data_std)

#%%
# Plot DBSCAN
plt.scatter(data.iloc[:,0], data.iloc[:,1], c=data_dbscan, cmap='prism', alpha=0.02)
#%%
plt.scatter(data_pca[:,0], data_pca[:,1], c=data_dbscan, cmap='prism', alpha=0.02)

# %%
# DBSCAN cluster number
n_dbscan = len(set(data_dbscan))
print(n_dbscan)
# %%
# Apply Gaussian Mixture Model
from sklearn.mixture import GaussianMixture

#%%
# GMM
gmm = GaussianMixture(n_components=1024)
gmm.fit(data_std)
data_gmm = gmm.predict(data_std)

#%%
# Plot GMM
plt.scatter(data.iloc[:,0], data.iloc[:,1], c=data_gmm, cmap='prism', alpha=0.02)

# %%


# %%
