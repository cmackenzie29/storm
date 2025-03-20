import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from scipy.spatial import ConvexHull

dir = "/Users/cameron/Downloads/135 control 4.22.24 001"

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

plt.plot()

# Plot colors together
colors = ["r", "g", "m"]
plt.subplot(1,2,1)
for (i, color) in enumerate([255, 65280, 16711935]):
	if color == 255:
		plt.scatter(positions[storm_colors==color][:,0], positions[storm_colors==color][:,1], color=colors[i], s=0.5)

# Plot clusters together
plt.subplot(1,2,2)
for clust in range(total_clusters):
    plt.scatter(positions[(clustering.labels_==clust)][:,0], positions[(clustering.labels_==clust)][:,1], s=0.5)

plt.show()







