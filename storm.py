import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from scipy.spatial import ConvexHull

condition_prefix = "135 control"
folder = "/Users/cameron/Downloads/4.22.24"

all_molec_counts = np.zeros((0,3))
all_densities = np.zeros((0,3))
all_intensities = np.zeros((0,3))

for file in os.listdir(folder):
    dir = os.path.join(folder, file)
    if os.path.isfile(dir) and (condition_prefix in file):
        print(file)
        with open(dir, 'rb') as f:
            contents = f.read()

        rows = str(contents,'utf-16').split('\r\n')[1:]
        storm_data = pd.DataFrame(data = [row.split('\t') for row in rows[1:(len(rows)-1)]], columns=rows[0].split('\t'))

        # Extract important parameters from dataframe
        positions = np.array(storm_data.iloc[:,1:4]).astype(float)
        storm_colors = np.array(storm_data.Color).astype(int)
        intensity = np.array(storm_data.Intensity).astype(float)

        ## DBSCAN
        clustering = DBSCAN(eps=0.1, min_samples=30, n_jobs=-1).fit(positions)
        total_clusters = max(clustering.labels_)+1

        ## Calculations
        molecs = np.zeros((total_clusters,3))
        densities = np.zeros((total_clusters,3))
        intensities = np.zeros((total_clusters,3))
        for clust in range(total_clusters):
            for (i, color) in enumerate([255, 65280, 16711935]):
                cluster = positions[(clustering.labels_==clust) & (storm_colors==color)]
                if len(cluster) >= 10: # Need at least 10 molecules in all three colors to count as a synapse for the analysis
                    if np.sum(cluster[:,2]==cluster[0,2])==len(cluster): #Cluster is all in one z-plane (only two dimensional)
                        densities[clust,i] = len(cluster)/ConvexHull(cluster[:,:2]).volume # Only in XY plane
                    else:
                        densities[clust,i] = len(cluster)/ConvexHull(cluster).volume
                    molecs[clust,i] = len(cluster)
                    intensities[clust,i] = np.sum(intensity[(clustering.labels_==clust) & (storm_colors==color)])

        included = (np.sum(molecs==0,axis=1)==0) # Array of indices where all three colors had enough molecules

        all_molec_counts = np.append(all_molec_counts,molecs[included],axis=0)
        all_densities = np.append(all_densities,densities[included],axis=0)
        all_intensities = np.append(all_intensities,intensities[included],axis=0)

print(f"Total synapses: {len(all_molec_counts)}")
print(f"Mean±SEM of molecule counts: {np.mean(all_molec_counts,axis=0)}±{np.std(all_molec_counts,axis=0)/np.sqrt(len(all_molec_counts))}")
print(f"Mean±SEM of cluster densities (#/µm2): {np.mean(all_densities,axis=0)}±{np.std(all_densities,axis=0)/np.sqrt(len(all_densities))}")
print(f"Mean±SEM of cluster intensities: {np.mean(all_intensities,axis=0)}±{np.std(all_intensities,axis=0)/np.sqrt(len(all_intensities))}")
print("(Order is red, green, magenta)")


