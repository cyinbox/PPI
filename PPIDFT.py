# -*- coding: utf-8 -*-
'''
# Programs for dentifying protein-protein interaction by Fourier transform
#
# Changchuan Yin, Ph.D.
# Dept. of Mathematics, Statistics and Computer Science
# University of Illinois at Chicago
# Chicago, IL 60607
# USA
#
# Email: cyin1@uic.edu, cyinbox@gmail.com
# Last updated on 02/01/2017
#
# Citation
# Yin, C. & Yau, Stephen S.-T. (2017).A Coevolution Analysis for Identifying Protein-Protein Interactions by Fourier Transform, PLOS ONE
# 
#
'''

import numpy as np
from scipy.fftpack import fft
from scipy.spatial import distance
import math
import os.path
from Bio import SeqIO

savePath = './PPIData'
#==============================================================================
# Function to return DFT real and image coefficients of a protein sequence
# Input: seq=protein sequence
# outputs: vRI=vector DFT coefficients of the protein sequence; ps = power spectrum

# The hydrophobicity value is from:
# Kyte, J., & Doolittle, R. F. (1982). A simple method for displaying the hydropathic character of a protein.
# Journal of molecular biology, 157(1), 105-132.
#==============================================================================
def DFTHydrophobicity(seq):
  AAHydrophobicity = {'A': 1.8, 'R': -4.5, 'N': -3.5,'D': -3.5,'C':2.5, 'Q': -3.5,'E':-3.5,\
 'G':-0.4,'H': -3.2,'I':4.5,'L':3.8,'K': -3.9,'M': 1.9,'F':2.8,'P':-0.9,\
 'S': -0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2,'U':0.0,'X':0.0} #X is for padding zeros
  n=len(seq)
  listAA=np.zeros(n) #define zeros list of length n
  
  for i in range(0,n):  # use index from 0 to n
    aa=seq[i]
    listAA[i]=AAHydrophobicity[aa]
  
  fftAA=fft(listAA)
  
  ps=np.abs(fftAA[0:n])
  ps=pow(ps,2)
  

  R=fftAA.real
  I=fftAA.imag
  
  vRI=[]

  for i in range (0,n):
    vRI.append(R[i])
  for i in range(0,n):
    vRI.append(I[i])
  
  return [ps,vRI]

#==============================================================================
# Function to compute the distance of two protein sequences in Euclidean space dimension size M
# Inputs: seq1, seq2 = protein sequences; M = vector space length
# outputs: dist = Euclidean distance of two proteins
#==============================================================================
def DFTProteinDistanceSpace(seq1,seq2,M):
 n1=len(seq1)
 n2=len(seq2)
 
 if n1<M:
    n=M-n1
    seq1 = ''.join((seq1,'X'*n))
    
 if n2<M:
    n=M-n2
    seq2 = ''.join((seq2,'X'*n))
    
 [ps1,vRI1]=DFTHydrophobicity(seq1)
 [ps2,vRI2]=DFTHydrophobicity(seq2)

 dist = distance.euclidean(vRI1,vRI2)

 return dist

#==============================================================================
# Function to compute the distance of two protein sequences in Euclidean space dimension size M
# Inputs: seqList = list of one protein sequences
# outputs: distM = dissimlarity distance matrix of mutual distance among proteins in the list
#          distV = dissimlarity vector the mutual distance among proteins in the list
#==============================================================================
def mutualDistanceDFTProteins(seqList):
 n=len(seqList)
 distV=[]
 distM=np.zeros([n, n]) 

 lenM=[]
 for i in range(0,n):
  lenM.append(len(seqList[i]))
 
 M=max(lenM)
 for i in range(0,n-1):
   seqA= seqList[i]
   for j in range (i+1,n):
      seqB= seqList[j]
      dist=DFTProteinDistanceSpace(seqA,seqB,M)
      distV.append(dist)
      distM[i,j]=dist
      
 return [distM,distV]

#==============================================================================
# Function to compute the DFT of a list of protein sequences of dimension size M
# Inputs: seqList = list of one protein sequences
# outputs: dftVectors = list of the DFT vectors for all the protein sequences in the list
#==============================================================================
def DFTProteinList(seqList):
 n=len(seqList)
 
 lenM=[]
 for i in range(0,n):
  lenM.append(len(seqList[i]))
 
 M=max(lenM)
 
 dftVectors=[]
 seqP=''
 for i in range(0,n):
    seq= seqList[i]
    s=len(seq)
    p=M-s
    seqP = ''.join((seq,'X'*p))
   
    [ps,vRI]=DFTHydrophobicity(seqP)
    dftVectors.append(vRI)
   
 return dftVectors

#==============================================================================
# Function to remove zeros from two lists if zeros are in the same positions. 
# Inputs: distA,distB = two lists of distances
# Output: distA, distB = two lists that were removed zeros if these zeros are in the same positions
#==============================================================================
def removeZerosSamePositions(distA,distB):
 n=len(distA)
 indexes=[]
 for i in range (0,n):
    if distA[i]==0 and distB[i]==0:
         indexes.append(i)
         
 for offset, index in enumerate(indexes):
   index -= offset
   del distA[index]
   del distB[index]

 return [distA, distB]

#==============================================================================
# Function to compute the score of PPI of two proteins 
# Inputs: seqListA, seqListB= list of proteinA sequences, list of proteinA sequences; 
# output: p = Peserson coefficient of the distance vectors of two proteins
#==============================================================================
def scorePPITreesP(seqListA,seqListB):
 [distMA,distVA]=mutualDistanceDFTProteins(seqListA)
 [distMA,distVB]=mutualDistanceDFTProteins(seqListB)
  
 [distVA, distVB]=removeZerosSamePositions(distVA,distVB)
 p=np.corrcoef(distVA, distVB)[0, 1] #Pearson coefficient

 if math.isnan(p):  
     p=0
 return p

#------------------------------------------------------------------------------
# Function to retrieve a specific gene/sequence from a fasta file
# Input: geneName=gene name,fastaFile=fasta file
# outputs: name,sequence
#------------------------------------------------------------------------------
def getAllSequences(proteinName):
 fastaFile=proteinName+'.fasta'
 fastaFile = os.path.join(savePath, fastaFile) 

 fasta_sequences = SeqIO.parse(open(fastaFile),'fasta')
 name=''
 sequence=''
 listSeq=[]

 for fasta in fasta_sequences:
  [name, sequence] = fasta.id, str(fasta.seq)
  listSeq.append(sequence)

 #print('List',listSeq)
 return listSeq

