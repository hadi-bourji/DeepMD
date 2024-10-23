#!/usr/bin/env python

#SBATCH -p batch
#SBATCH -o myMPI.o%j    # Output and error file name (%j expands to jobID)
#SBATCH -N 1 --ntasks-per-node=24
#SBATCH -t 13-1:00:00
#SBATCH --mem 60G
#SBATCH --error=quip80.err --output=quip80.out
#SBATCH -J quip80

from quippy.descriptors import Descriptor
import ase.io
import numpy as np
import glob
import matplotlib.pyplot as plt
from ase.io import read

desc = Descriptor("soap l_max=6 n_max=12 cutoff =5 atom_sigma=.5") #set soap parameters quip
soap_vectors = []
for i in glob.glob("exyzFiles/*"): 
    print(i)
    zeolite = ase.io.read(filename=i,format='vasp')
    D2=desc.calc(zeolite)['data']
    soap_vectors.append(D2.flatten())
    
soap_vectors = np.array(soap_vectors)  
