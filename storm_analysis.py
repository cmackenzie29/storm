import os
from tqdm import tqdm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from scipy.spatial import ConvexHull
from sklearn.neighbors import NearestNeighbors, KDTree

def load_storm_file(folder, f, colors):
    directory = os.path.join(folder, f)
    with open(directory, 'rb') as f:
        contents = f.read()
    rows = str(contents,'utf-16').split('\r\n')[1:]
    storm_data = pd.DataFrame(data = [row.split('\t') for row in rows[1:(len(rows)-1)]], columns=rows[0].split('\t'))
    
    positions = np.array(storm_data.iloc[:,1:4]).astype(float)
    storm_colors = np.array(storm_data.Color).astype(int)
    intensities = np.array(storm_data.Intensity).astype(float)

    # Filter by user-requested colors
    col_filter = np.column_stack([(storm_colors == x) for x in colors.keys()]).any(axis=1)
    return positions[col_filter, :], storm_colors[col_filter], intensities[col_filter]


def display_clusters(positions, eps, min_samples, storm_colors=None):
    # Pass in no colors for clustering mode "a", where all colors are assumed to exist within a cluster
    if type(storm_colors) == type(None):
        plt.scatter(positions[:,0], positions[:,1], s=0.1, color='k', alpha=0.1)
        clustering = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit(positions)
        for clust in range(max(clustering.labels_)+1):
            plt.scatter(positions[(clustering.labels_==clust)][:,0], positions[(clustering.labels_==clust)][:,1], s=0.5)
        plt.show()
    else:
        plt.scatter(positions[:,0], positions[:,1], s=0.1, color='k', alpha=0.1)
        for c in list(set(storm_colors)):
            clustering = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit(positions[storm_colors==c, :])
            for clust in range(max(clustering.labels_)+1):
                plt.scatter(positions[storm_colors==c, :][(clustering.labels_==clust)][:,0], positions[storm_colors==c, :][(clustering.labels_==clust)][:,1], s=0.5)
        plt.show()


def analyze_basic(folder, files, condition_prefix, colors, output_mode):
    nn_dists = {}
    # Initialize color combos
    for c1 in colors.keys():
        for c2 in colors.keys():
            nn_dists[f"{c1}-{c2}"] = np.zeros(0)

    for file in tqdm(files):
        positions, storm_colors, _ = load_storm_file(folder, file, colors)

        # Nearest neighbors within each color
        for c in colors.keys():
            nbrs = NearestNeighbors(n_neighbors=2, algorithm='kd_tree').fit(positions[storm_colors==c])
            dists, _ = nbrs.kneighbors(positions[storm_colors==c])
            nn_dists[f"{c}-{c}"] = np.append(nn_dists[f"{c}-{c}"], dists[:,1])

        # Nearest neighbor to different colors
        for c2 in colors.keys():
            c2_tree = KDTree(positions[storm_colors==c2, :])
            for c1 in colors.keys():
                if c1 == c2:
                    continue
                dists, _ = c2_tree.query(positions[storm_colors==c1, :], k=1)
                nn_dists[f"{c1}-{c2}"] = np.append(nn_dists[f"{c1}-{c2}"], dists.flatten())

    # Summary statistics
    output_text = ""
    for k in nn_dists.keys():
        c1, c2 = [int(x) for x in k.split("-")]
        output_text += f"Distance to nearest {colors[c2]} molecule from each {colors[c1]} molecule:\n"
        output_text += f"Mean: {nn_dists[k].mean()}\n"
        output_text += f"SEM: {nn_dists[k].std() / np.sqrt(len(nn_dists[k]))}\n"
        output_text += f"N: {len(nn_dists[k])}\n\n"

        if output_mode == "r": # Output csvs of raw data
            np.savetxt(os.path.join(folder, f"results/{condition_prefix}_nndists_basic_from_{colors[c1]}_to_{colors[c2]}.csv"), nn_dists[k], header = f"nn_from_{colors[c1]}_to_{colors[c2]}", delimiter = ",", comments = "")

    with open(os.path.join(folder, f"results/{condition_prefix}_nndists_basic_summary.txt"), "w") as f:
        f.write(output_text)

    print(f"\nAnalysis complete. Results are in \"{folder}/results\"")


def analyze_clustering_colors_together(folder, files, condition_prefix, dim, colors, eps, min_samples, output_mode):
    nn_dists = {}
    # Initialize color combos
    for c1 in colors.keys():
        for c2 in colors.keys():
            nn_dists[f"{c1}-{c2}"] = np.zeros(0)

    molec_counts = {}
    densities = {}
    volumes = {}
    intens = {}
    for c in colors.keys():
        molec_counts[c] = np.zeros(0)
        densities[c] = np.zeros(0)
        volumes[c] = np.zeros(0)
        intens[c] = np.zeros(0)

    for file in tqdm(files):
        positions, storm_colors, intensities = load_storm_file(folder, file, colors)

        # Begin plotting
        plt.clf()
        plt.scatter(positions[:,0], positions[:,1], s=0.1, color='k', alpha=0.1)

        # Clustering
        clustering = np.array(DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit(positions).labels_)
        total_clusters = int(clustering.max())+1

        # Remove noisy molecules
        positions = positions[clustering != -1, :]
        storm_colors = storm_colors[clustering != -1]
        intensities = intensities[clustering != -1]
        clustering = clustering[clustering != -1]

        # Nearest neighbors within each color
        for c in colors.keys():
            nbrs = NearestNeighbors(n_neighbors=2, algorithm='kd_tree').fit(positions[storm_colors==c])
            dists, _ = nbrs.kneighbors(positions[storm_colors==c])
            nn_dists[f"{c}-{c}"] = np.append(nn_dists[f"{c}-{c}"], dists[:,1])

        # Nearest neighbor to different colors
        for c2 in colors.keys():
            c2_tree = KDTree(positions[storm_colors==c2, :])
            for c1 in colors.keys():
                if c1 == c2:
                    continue
                dists, _ = c2_tree.query(positions[storm_colors==c1, :], k=1)
                nn_dists[f"{c1}-{c2}"] = np.append(nn_dists[f"{c1}-{c2}"], dists.flatten())

        # Synapse calculations
        for clust in range(total_clusters):
            pos_clust = positions[clustering == clust, :]
            cols_clust = storm_colors[clustering == clust]
            ints_clust = intensities[clustering == clust]
            
            # Check that cluster contains at least dim+2 points for all colors (min to determine area/volume)
            col_cts = np.array([np.sum(cols_clust == x) for x in colors.keys()])
            if (len(col_cts) == len(colors)) and (col_cts >= dim+2).all(): 
                plt.scatter(pos_clust[:,0], pos_clust[:,1], s=0.5) # Plot cluster
                for c in colors.keys():
                    if dim == 2: # 2D STORM
                        vol = ConvexHull(pos_clust[cols_clust==c, :2]).volume # Only XY plane
                    else: # 3D STORM
                        vol = ConvexHull(pos_clust[cols_clust==c, :]).volume
                    molecs = np.sum(cols_clust==c)

                    molec_counts[c] = np.append(molec_counts[c], molecs)
                    densities[c] = np.append(densities[c], molecs / vol)
                    volumes[c] = np.append(volumes[c], vol)
                    intens[c] = np.append(intens[c], np.sum(ints_clust[cols_clust==c]))

        # Finish plotting for file
        plt.savefig(os.path.join(folder, f"results/{'.'.join(file.split('.')[:-1])}_clustering.png"))

    # Summary statistics for NNs
    output_text = ""
    for k in nn_dists.keys():
        c1, c2 = [int(x) for x in k.split("-")]
        output_text += f"Distance to nearest {colors[c2]} molecule from each {colors[c1]} molecule:\n"
        output_text += f"Mean: {nn_dists[k].mean()}\n"
        output_text += f"SEM: {nn_dists[k].std() / np.sqrt(len(nn_dists[k]))}\n"
        output_text += f"N: {len(nn_dists[k])}\n\n"

        if output_mode == "r": # Output csvs of raw data
            np.savetxt(os.path.join(folder, f"results/{condition_prefix}_nndists_from_{colors[c1]}_to_{colors[c2]}.csv"), nn_dists[k], header = f"nn_from_{colors[c1]}_to_{colors[c2]}", delimiter = ",", comments = "")

    with open(os.path.join(folder, f"results/{condition_prefix}_nndists_summary.txt"), "w") as f:
        f.write(output_text)

    # Summary statistics for synapse calculations
    output_text = ""
    for k in molec_counts.keys():
        output_text += f"{colors[k]} molecule counts per synapse:\n"
        output_text += f"Mean: {molec_counts[k].mean()}\n"
        output_text += f"SEM: {molec_counts[k].std() / np.sqrt(len(molec_counts[k]))}\n"
        output_text += f"N: {len(molec_counts[k])}\n\n"

        output_text += f"{colors[k]} molecule densities per synapse:\n"
        output_text += f"Mean: {densities[k].mean()}\n"
        output_text += f"SEM: {densities[k].std() / np.sqrt(len(densities[k]))}\n"
        output_text += f"N: {len(densities[k])}\n\n"

        output_text += f"{colors[k]} synapse volumes:\n"
        output_text += f"Mean: {volumes[k].mean()}\n"
        output_text += f"SEM: {volumes[k].std() / np.sqrt(len(volumes[k]))}\n"
        output_text += f"N: {len(volumes[k])}\n\n"

        output_text += f"{colors[k]} molecule intensities per synapse:\n"
        output_text += f"Mean: {intens[k].mean()}\n"
        output_text += f"SEM: {intens[k].std() / np.sqrt(len(intens[k]))}\n"
        output_text += f"N: {len(intens[k])}\n\n"

        if output_mode == "r": # Output csvs of raw data
            np.savetxt(os.path.join(folder, f"results/{condition_prefix}_{colors[k]}_molec_counts.csv"), molec_counts[k], header = f"{colors[k]}_molec_counts", delimiter = ",", comments = "")
            np.savetxt(os.path.join(folder, f"results/{condition_prefix}_{colors[k]}_densities.csv"), densities[k], header = f"{colors[k]}_densities", delimiter = ",", comments = "")
            np.savetxt(os.path.join(folder, f"results/{condition_prefix}_{colors[k]}_volumes.csv"), volumes[k], header = f"{colors[k]}_volumes", delimiter = ",", comments = "")
            np.savetxt(os.path.join(folder, f"results/{condition_prefix}_{colors[k]}_intensities.csv"), intens[k], header = f"{colors[k]}_intensities", delimiter = ",", comments = "")

    with open(os.path.join(folder, f"results/{condition_prefix}_synapses_summary.txt"), "w") as f:
        f.write(output_text)

    print(f"\nAnalysis complete. Results are in \"{folder}/results\"")


def analyze_clustering_colors_individual(folder, files, condition_prefix, dim, colors, eps, min_samples, output_mode):
    nn_dists = {}
    # Initialize color combos
    for c1 in colors.keys():
        for c2 in colors.keys():
            nn_dists[f"{c1}-{c2}"] = np.zeros(0)

    molec_counts = {}
    densities = {}
    volumes = {}
    intens = {}
    for c in colors.keys():
        molec_counts[c] = np.zeros(0)
        densities[c] = np.zeros(0)
        volumes[c] = np.zeros(0)
        intens[c] = np.zeros(0)

    for file in tqdm(files):
        positions, storm_colors, intensities = load_storm_file(folder, file, colors)

        # Begin plotting
        plt.clf()
        plt.scatter(positions[:,0], positions[:,1], s=0.1, color='k', alpha=0.1)

        # Take one color at a time for clustering
        clustering = np.zeros(len(positions))
        for c in colors.keys():
            # Clustering
            clustering_c = np.array(DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit(positions[storm_colors==c,:]).labels_)
            if (clustering==-2).all(): # First time through
                clustering[storm_colors==c] = clustering_c
            else: # -1's should stay that way. New color should be incremented above current maximum value
                clustering_c[clustering_c != -1] = clustering_c[clustering_c != -1] + 1 + clustering.max()
                clustering[storm_colors==c] = clustering_c
        
        total_clusters = int(clustering.max())+1

        # Remove noisy molecules
        positions = positions[clustering != -1, :]
        storm_colors = storm_colors[clustering != -1]
        intensities = intensities[clustering != -1]
        clustering = clustering[clustering != -1]

        # Nearest neighbors within each color
        for c in colors.keys():
            nbrs = NearestNeighbors(n_neighbors=2, algorithm='kd_tree').fit(positions[storm_colors==c])
            dists, _ = nbrs.kneighbors(positions[storm_colors==c])
            nn_dists[f"{c}-{c}"] = np.append(nn_dists[f"{c}-{c}"], dists[:,1])

        # Nearest neighbor to different colors
        for c2 in colors.keys():
            c2_tree = KDTree(positions[storm_colors==c2, :])
            for c1 in colors.keys():
                if c1 == c2:
                    continue
                dists, _ = c2_tree.query(positions[storm_colors==c1, :], k=1)
                nn_dists[f"{c1}-{c2}"] = np.append(nn_dists[f"{c1}-{c2}"], dists.flatten())

        # Synapse calculations
        for clust in range(total_clusters):
            pos_clust = positions[clustering == clust, :]
            c = storm_colors[clustering == clust] # Single color
            ints_clust = intensities[clustering == clust]

            # Check that cluster contains at least dim+2 points (min to determine area/volume)
            if len(pos_clust) >= dim+2: 
                c = c[0]
                plt.scatter(pos_clust[:,0], pos_clust[:,1], s=0.5) # Plot cluster
                
                if dim == 2: # 2D STORM
                    vol = ConvexHull(pos_clust[:, :2]).volume # Only XY plane
                else: # 3D STORM
                    vol = ConvexHull(pos_clust).volume
                molecs = len(pos_clust)

                molec_counts[c] = np.append(molec_counts[c], molecs)
                densities[c] = np.append(densities[c], molecs / vol)
                volumes[c] = np.append(volumes[c], vol)
                intens[c] = np.append(intens[c], np.sum(ints_clust))

        # Finish plotting for file
        plt.savefig(os.path.join(folder, f"results/{'.'.join(file.split('.')[:-1])}_clustering_single_color.png"))

    # Summary statistics for NNs
    output_text = ""
    for k in nn_dists.keys():
        c1, c2 = [int(x) for x in k.split("-")]
        output_text += f"Distance to nearest {colors[c2]} molecule from each {colors[c1]} molecule:\n"
        output_text += f"Mean: {nn_dists[k].mean()}\n"
        output_text += f"SEM: {nn_dists[k].std() / np.sqrt(len(nn_dists[k]))}\n"
        output_text += f"N: {len(nn_dists[k])}\n\n"

        if output_mode == "r": # Output csvs of raw data
            np.savetxt(os.path.join(folder, f"results/{condition_prefix}_nndists_single_color_from_{colors[c1]}_to_{colors[c2]}.csv"), nn_dists[k], header = f"nn_from_{colors[c1]}_to_{colors[c2]}", delimiter = ",", comments = "")

    with open(os.path.join(folder, f"results/{condition_prefix}_nndists_single_color_summary.txt"), "w") as f:
        f.write(output_text)

    # Summary statistics for synapse calculations
    output_text = ""
    for k in molec_counts.keys():
        output_text += f"{colors[k]} molecule counts per synapse:\n"
        output_text += f"Mean: {molec_counts[k].mean()}\n"
        output_text += f"SEM: {molec_counts[k].std() / np.sqrt(len(molec_counts[k]))}\n"
        output_text += f"N: {len(molec_counts[k])}\n\n"

        output_text += f"{colors[k]} molecule densities per synapse:\n"
        output_text += f"Mean: {densities[k].mean()}\n"
        output_text += f"SEM: {densities[k].std() / np.sqrt(len(densities[k]))}\n"
        output_text += f"N: {len(densities[k])}\n\n"

        output_text += f"{colors[k]} synapse volumes:\n"
        output_text += f"Mean: {volumes[k].mean()}\n"
        output_text += f"SEM: {volumes[k].std() / np.sqrt(len(volumes[k]))}\n"
        output_text += f"N: {len(volumes[k])}\n\n"

        output_text += f"{colors[k]} molecule intensities per synapse:\n"
        output_text += f"Mean: {intens[k].mean()}\n"
        output_text += f"SEM: {intens[k].std() / np.sqrt(len(intens[k]))}\n"
        output_text += f"N: {len(intens[k])}\n\n"

        if output_mode == "r": # Output csvs of raw data
            np.savetxt(os.path.join(folder, f"results/{condition_prefix}_{colors[k]}_molec_counts_single_color_clusters.csv"), molec_counts[k], header = f"{colors[k]}_molec_counts", delimiter = ",", comments = "")
            np.savetxt(os.path.join(folder, f"results/{condition_prefix}_{colors[k]}_densities_single_color_clusters.csv"), densities[k], header = f"{colors[k]}_densities", delimiter = ",", comments = "")
            np.savetxt(os.path.join(folder, f"results/{condition_prefix}_{colors[k]}_volumes_single_color_clusters.csv"), volumes[k], header = f"{colors[k]}_volumes", delimiter = ",", comments = "")
            np.savetxt(os.path.join(folder, f"results/{condition_prefix}_{colors[k]}_intensities_single_color_clusters.csv"), intens[k], header = f"{colors[k]}_intensities", delimiter = ",", comments = "")

    with open(os.path.join(folder, f"results/{condition_prefix}_synapses_single_color_summary.txt"), "w") as f:
        f.write(output_text)

    print(f"\nAnalysis complete. Results are in \"{folder}/results\"")


