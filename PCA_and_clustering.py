#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import matplotlib
from MDAnalysis.analysis import align
get_ipython().run_line_magic('matplotlib', 'inline')
matplotlib.rcParams['figure.dpi'] = 100

# The line of code sets the default resolution for plots created using matplotlib, specifically for figures. 
# The default value for figure.dpi is 100, meaning that each inch of the figure is 100 pixels wide and 100 pixels tall.


# In[5]:


sys.path.append('/path/to/your/folder')


# In[6]:


from TrajectoryPCA import TrajectoryPCA


# In[7]:


# This is an example with ATG4B system ( 2776 heavy atoms )

model = mda.Universe(os.path.join(path,'trajectory.dcd'))
#residue_indices = [42, 43, 44, 46, 47, 48, 49, 50, 69, 71, 72, 129, 130, 131, 132, 168, 169, 170, 172, 173, 194, 195, 197, 198, 201, 202, 203, 205, 206]
#selection_string = ' or '.join(f'resid {idx}' for idx in residue_indices)
protein = model.select_atoms('backbone') #atom selection to use for superimposition

#x_ray = mda.Universe(os.path.join(path,'GT_min_fixed.pdb'))
#x_ray_selected = x_ray.select_atoms('protein and backbone')


# In[8]:


protein_PCA = TrajectoryPCA(protein)


# In[9]:


newMeanStructure = protein_PCA.superimpose2mean()


# In[10]:


cum_explained_var_ratio, pca_components = protein_PCA.getPCA()


# In[11]:


pcs, pca_2_components_ = protein_PCA.project2PC12()


# In[12]:


plt.scatter(pcs[:, 0], pcs[:, 1], marker='^', color='green')
for i, point in enumerate(pcs):
    plt.text(point[0], point[1], f'{i+1}', fontsize=12, ha='left')

# Add labels and title
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('Scatter plot with point labels')
plt.savefig(os.path.join(path,'frames_pca.png'),dpi=600)


# In[41]:


import matplotlib.pyplot as plt

pca_num = range(1, len(cum_explained_var_ratio)+1)
plt.plot(pca_num, cum_explained_var_ratio)
plt.xlabel('# principal components')
plt.ylabel('cumulative explained variance');

# Annotate a specific point in the figure, here we specify the first two principal components explanation value.

index = 2
point = (index, cum_explained_var_ratio[index-1])
plt.plot(point[0], point[1], 'o', color='red')
plt.annotate(f'Point {index} ({index}, {cum_explained_var_ratio[index-1]:.2f})', 
             xy=point, xytext=(index + 10, cum_explained_var_ratio[index-1] + 0.02),
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.5));

# The second argument, xy, specifies the position of the point on the plot where the annotation will be placed. 
# In this case, the value is point, which is a tuple (index, cum_explained_var_ratio[index]) that represents 
# the x and y values of the point to be annotated.

# The third argument, xytext, specifies the position of the text relative to the point. In this case, the value 
# is (index + 5, cum_explained_var_ratio[index] + 0.05), which places the text to the right and above the point.


# # Implement silhouette_avg score to best choose the number of clusters (clustering method = KMeans)

# In[42]:


import matplotlib.pyplot as plt
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


# In[43]:


# Plot the silhouette_avg score along different number of clusters

import matplotlib
get_ipython().run_line_magic('matplotlib', 'inline')
matplotlib.rcParams['figure.dpi'] = 100

n_clusters = 10
silhouette_avg = []

fig = plt.figure(figsize=(10, 10))

# define the number of rows and columns for the grid
fig_x = 3
fig_y = 3
fig_num = fig_x * fig_y

# create a new figure object
fig = plt.figure(figsize=(8, 8))

# set the spacing between the subplots
fig.subplots_adjust(hspace=1, wspace=1)

# loop over each subplot in the grid
for i in range(2, n_clusters + 1):
    kmeans_fit = KMeans(init = 'k-means++', n_clusters = i, n_init = 20).fit(pcs)
    silhouette_avg.append(silhouette_score(pcs, kmeans_fit.labels_))
    
    # add a subplot to the figure
    ax = fig.add_subplot(fig_x, fig_y, i-1)
    
    # plot things on the subplot
    n_clusters_IDs = []
    for j in range(i):
        n_clusters_IDs.append(np.linalg.norm(pcs-kmeans_fit.cluster_centers_[j], axis=1).argmin())
        ax.scatter(pcs.T[0][kmeans_fit.labels_==j], pcs.T[1][kmeans_fit.labels_==j], marker='.')
        ax.scatter(pcs[n_clusters_IDs[j]][0],pcs[n_clusters_IDs[j]][1],marker='^')
        ax.annotate('pca'+str(j+1), [pcs[n_clusters_IDs[j]][0]+2, pcs[n_clusters_IDs[j]][1]+2])
    
    # set the title for the subplot
    ax.set_title(f"PCA with {i} clusters")
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    
    # print n_cluster_IDs
    print(f"{i}_clusters_IDs is: {n_clusters_IDs}.")
    
# show the figure
plt.show()


    
# n_init is the number of times the KMeans algorithm will be run with different 
# centroid seeds. The final results will be the best output of n_init consecutive 
# runs in terms of inertia, i.e., the sum of squared distances of samples to their 
# closest cluster center.


# In[44]:


plt.plot(range(2, n_clusters + 1), silhouette_avg)

# The xy parameter specifies the coordinates of the annotated point, and the xytext parameter specifies the 
# coordinates of the text label. The arrowprops parameter adds an arrow pointing to the annotated point. 
# Finally, the fontsize parameter sets the font size of the text label.

kmeans_ymax = max(silhouette_avg)
kmeans_xmax = silhouette_avg.index(kmeans_ymax) + 2
plt.annotate('max point: ({:.0f}, {:.2f})'.format(kmeans_xmax, kmeans_ymax), 
             xy=(kmeans_xmax, kmeans_ymax), 
             xytext=(kmeans_xmax, kmeans_ymax+0.05), 
             arrowprops=dict(facecolor='black', shrink=1), 
             fontsize=12)

plt.xticks(range(2, n_clusters + 1), range(2, n_clusters + 1))
plt.xlabel('number of clusters')
plt.ylabel('silhouette_avg')
plt.savefig(os.path.join(path,'Silhoutte_coefficient.png'),dpi=600)
plt.show()


# # For optimal choices of K means

# In[46]:


from sklearn.cluster import KMeans

n_clusters = kmeans_xmax
n_clusters_IDs = []

kmean = KMeans(init = 'k-means++', n_clusters = n_clusters, n_init = 20).fit(pcs)
#kmean = KMeans(init = 'random', n_clusters = n_clusters, n_init = 'auto').fit(pcs)
#print(kmean)

fig = plt.figure() 
ax = fig.add_subplot(111)
ax.set_title('PCA of the GT peptide for 200ns MD simulation') 

for i in range(n_clusters):
    n_clusters_IDs.append(np.linalg.norm(pcs-kmean.cluster_centers_[i], axis=1).argmin())
    ax.scatter(pcs.T[0][kmean.labels_==i], pcs.T[1][kmean.labels_==i], marker='.')
    ax.scatter(pcs[n_clusters_IDs[i]][0], pcs[n_clusters_IDs[i]][1],marker='^')
    ax.annotate('clust'+str(i+1), [pcs[n_clusters_IDs[i]][0] - 1.0, pcs[n_clusters_IDs[i]][1] + 0.7], fontweight='bold')

plt.xlabel('PC1')
plt.ylabel('PC2')
#plt.plot(new_x_ray_x, new_x_ray_y, '^', color='black', label='model5')
#plt.annotate('model5', xy=point, xytext=(new_x_ray_x - 1.5, new_x_ray_y + 0.7), fontweight='bold');

plt.savefig(os.path.join(path,'pca_clustering.png'),dpi=600)
plt.show() 

print("n_clusters_IDs:", n_clusters_IDs)

#print(pcs-kmean.cluster_centers_[1])

