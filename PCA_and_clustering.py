import sys
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import matplotlib
from MDAnalysis.analysis import align

path = sys.argv[1]
top = sys.argv[2]
traj = sys.argv[3]


from TrajectoryPCA import TrajectoryPCA

model = mda.Universe(os.path.join(path,top),os.path.join(path,traj))
protein = model.select_atoms('backbone') #atom selection to use for superimposition

x_ray = mda.Universe(os.path.join(path,'prot_solv_ions.pdb'))
x_ray_selected = x_ray.select_atoms('protein and backbone')

protein_PCA = TrajectoryPCA(protein)

newMeanStructure = protein_PCA.superimpose2mean()

cum_explained_var_ratio, pca_components = protein_PCA.getPCA()

pcs, pca_2_components_ = protein_PCA.project2PC12()

#Experimentally resolved or Predicted structure
# Initiation variables

x_ray_Position = x_ray_selected.positions - protein.center_of_geometry()
R_x_ray, rmsd_x_ray = align.rotation_matrix(x_ray_Position, newMeanStructure)
new_x_ray_Position = np.array(R_x_ray.dot(x_ray_Position.T).T)

x_ray_Position = x_ray_selected.positions - protein.center_of_geometry()
R_x_ray, rmsd_x_ray = align.rotation_matrix(x_ray_Position, newMeanStructure)
new_x_ray_Position = np.array(R_x_ray.dot(x_ray_Position.T).T)

# Reshape (3 x N) to 3N
shape = new_x_ray_Position.shape
flatPos_x_ray = new_x_ray_Position.reshape(shape[0]*shape[1])
flatPos_mean = newMeanStructure.reshape(shape[0]*shape[1])

# Compute the dot products
x_ray_x = flatPos_x_ray.dot(pca_2_components_[0])
x_ray_y = flatPos_x_ray.dot(pca_2_components_[1])
mean_x = flatPos_mean.dot(pca_2_components_[0])
mean_y = flatPos_mean.dot(pca_2_components_[1])

new_x_ray_x = x_ray_x - mean_x
new_x_ray_y = x_ray_y - mean_y

#Plot the experimentally resolved or predicted structure on the PCA space
plt.plot(new_x_ray_x, new_x_ray_y, '^', color='black')

#Plot the MD frames on the PCA space
plt.scatter(pcs[:, 0], pcs[:, 1], marker='^', color='green')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.savefig(os.path.join(path,'frames_pca.png'),dpi=600)


