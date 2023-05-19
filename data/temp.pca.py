# PCA
#%%
import pandas
import numpy
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

#%%
# Load data
data = pandas.read_csv(
    '/Users/hbk/Projects/Lab/06_FoldMotif/repos/motifsearch/data/temp.vec.data.raw.tsv', sep='\t', header=0)
#%%
data
#%%

# Standardize data
data_std = StandardScaler().fit_transform(data)
#%%

# PCA
pca = PCA(n_components=2)
pca.fit(data_std)
data_pca = pca.transform(data_std)
#%%

# Plot
# plt.scatter(data_pca[:,0], data_pca[:,1], alpha=0.02)

plt.scatter(data.iloc[:,0], data.iloc[:,1], alpha=0.02)
plt.show()


# %%
