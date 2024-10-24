#!/usr/bin/env python

#SBATCH -p batch
#SBATCH -o myMPI.o%j    # Output and error file name (%j expands to jobID)
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH -t 13-1:00:00
#SBATCH --mem 60G
#SBATCH --output=deformations_refinement

from quippy.descriptors import Descriptor
import ase.io
import numpy as np
import glob
import matplotlib.pyplot as plt
from ase.io import read
from sklearn.cluster import KMeans
import random
import os
import shutil


desc = Descriptor("soap l_max=6 n_max=12 cutoff =5 atom_sigma=.5")
# groups = ["amorphousSilica/3000K","amorphousSilica/4000K","amorphousSilica/5000K","amorphousSilica/6000K","amorphousSilica/4800K"]

# subs = ["COE","CRI","HAQ","HSB","MOG","QUA","STI","TRI","CHA","IRR","MTF","MVY","SOD"]
subs = ["8017444","8025591","8034621","8059807","8067891","8069993","8072794","8076949","8129288","8129504"]
for j in subs:
#     shutil.rmtree(f"/project/palmer/hadi/metagga/kmeans_one_eighth_dataset/deformations/{j}")
    os.makedirs(f"/project/palmer/hadi/metagga/kmeans_one_eighth_dataset/refinement/deform{j}", exist_ok=True)
    soap_vectors = []
    filenames = []
    total = int(len(glob.glob(f"/project/palmer/hadi/nnp/fulldft/calculations/hypotheticalZeolites/deformations/VASP/data*{j}*/POSCAR"))/8)
    for i in glob.glob(f"/project/palmer/hadi/nnp/fulldft/calculations/hypotheticalZeolites/deformations/VASP/data*{j}*/POSCAR"): 
        print(i)
        zeolite = ase.io.read(filename=i,format='vasp')
        D2=desc.calc(zeolite)['data']
        soap_vectors.append(D2.flatten())
        filenames.append(i)

    soap_vectors = np.array(soap_vectors)  
    kmeans = KMeans(n_clusters=total, random_state=0).fit(soap_vectors)
    labels = kmeans.labels_


    cluster_dict = {i: [] for i in range(total)}
    for filename, label in zip(filenames, labels):
        cluster_dict[label].append(filename)

    empty_clusters = [i for i in range(total) if len(cluster_dict[i]) == 0]
    if empty_clusters:
        print(f"Warning: The following clusters are empty: {empty_clusters}")


    selected_structures = []

    results = list(zip(filenames, labels))
    for filename, label in results:
        print(f"Filename: {filename}, Cluster Label: {label}")

    for label in cluster_dict:
        if cluster_dict[label]:  # If the cluster is not empty
            selected_structure = random.choice(cluster_dict[label])
            selected_structures.append(selected_structure)

    while len(selected_structures) < total:
        non_empty_clusters = [k for k in cluster_dict.keys() if cluster_dict[k]]  # Get non-empty clusters
        if not non_empty_clusters:
            break  # Exit if no more non-empty clusters
        chosen_cluster = random.choice(non_empty_clusters)
        selected_structure = random.choice(cluster_dict[chosen_cluster])
        selected_structures.append(selected_structure)

    print("\nSelected structures:")
    for structure in selected_structures:
        os.makedirs(f"/project/palmer/hadi/metagga/kmeans_one_eighth_dataset/refinement/deform{j}/{os.path.basename(os.path.dirname(structure))}", exist_ok=True)
        shutil.copy(structure, f"/project/palmer/hadi/metagga/kmeans_one_eighth_dataset/refinement/deform{j}/{os.path.basename(os.path.dirname(structure))}")
        print(structure)
