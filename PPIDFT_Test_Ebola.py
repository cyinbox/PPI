# -*- coding: utf-8 -*-
'''
# Test program of protein-protein interaction on Ebola genomes
#
# Changchuan Yin, Ph.D.
# Dept. of Mathematics, Statistics and Computer Science
# University of Illinois at Chicago
# Chicago, IL 60607
# USA
#
# Email cyin1@uic.edu, cyinbox@gmail.com
# Last update 02/06/2017
#
# Citation
# Yin, C. & Yau, Stephen S.-T. (2017). A coevolution analysis for identifying protein-protein interactions 
# by Fourier transform. PLOS ONE.
#
'''
import numpy as np
from matplotlib import pyplot as plt
from sklearn import manifold
import PPIDFT as ppi

# Ebola proteins files, GP.fasta, NP.fasta, VP30.fasta, VP35.fasta, VP40.fasta, and L.fasta,
# shall exist in the sub-folder ./PPIData. These files contain the corresponding protein sequences from ebola genomes
proteinNames=['GP','NP','VP24','VP30','VP35','VP40','L']
n=len(proteinNames)
distM=np.zeros([n, n]) 
distV=[]

for i in range(0,n):
   nameA=proteinNames[i]
   proteinsA=ppi.getAllSequences(nameA) #The Ebola virus file: one file contain the same protein for different geneomes
   for j in range (0,n):
      nameB=proteinNames[j]
      proteinsB=ppi.getAllSequences(nameB)
      dist=1-ppi.scorePPITreesP(proteinsA,proteinsB)
      distV.append(dist)
      distM[i,j]=dist
            
'''
Multidimensional Scaling (MDS) analysis of disstances of PPIs. It maps distance matrix to 2-D
space for the PPI visualization.
'''
# Dimension is 2, the MDS fitting may be different from each run.
mds=manifold.MDS(n_components=2, max_iter=500,dissimilarity="precomputed", n_jobs=1)
pos=mds.fit(distM).embedding_

print('MDS Positions:\n',pos)
labels=['GP','NP','VP24','VP30','VP35','VP40','L']

# Plot the points
markerSizes=[676,739,251,289,341,326,1210] #Illustration for realative lengths of amino acid sequences
v=[0.02*j for j in markerSizes]
colors=['r','b','k','c','y','m','g']

for i in range(len(pos)):
    plt.plot(pos[i, 0], pos[i, 1],marker='o',c=colors[i],markersize=v[i])
    plt.text(pos[i, 0]+0.008, pos[i, 1]+0.008, proteinNames[i], fontsize=10) # a little off the point positions

plt.xlabel('Relative interaction', labelpad=5,fontsize=11)
plt.ylabel('Relative interaction',labelpad=5,fontsize=11)
plt.savefig('./PPIData/PPIDFT_Ebola.eps', dpi=200) 
plt.show()
print('Completed')

