# # Here try to do Striatalcortical parcellation in python
# 
# The process works like this. 
# 
# * You will need a set of ROI that consists of: (1) the cortex, (2) a subcortical mask that you wish to parcellate based on its connectivity profile with the cortex, in this case the Striatal. 
# * The cortical ROIs should be labeled as integer in a NIFTI file. The subcortical mask should also have each voxel labeled as different integer values. Note the values cannot overlap with corticla ROIs.
# * This combined ROI set will then be used to generate a subcortical+cortical full correlation matrix (adjacency matrix). Essentially BOLD timesires from each ROI and subcortical voxel will be correlated with each other using 3dNetCorr.
# * The output adjancency matrix from each subject can then be averaged. Alternatively it is possible to do this parcellation on a subject by subject basis, but its not implemented yet.
# * For each subcortical voxel, its connectivity strength with different cortical ROIs will then be sorted (ranked). 
# * You can then classify the subcortical voxel based which cortical region (network, module) its most strongly connected to.

# ROI used in the following example
# * Total of 3836 ROIs in the Striatal + cortex parcellation file
# * 3539 voxels in thalammic parcellation files
# * 297 ROIs in cortical ROI file. 

# In[86]:

import numpy as np
import scipy as sp
import scipy.io as sio
import matplotlib.pyplot as plt
import pickle
import csv
import glob
from __future__ import division, print_function


# Averaging
AdjMat_Files = glob.glob('/home/despoB/kaihwang/Rest/AdjMatrices/*Ses1_FIX_striatalcortical_corrmat')

M_Sum =[];
for f in AdjMat_Files:
    M = np.loadtxt(f)
    M_Sum += [M]
    
    
AdjMat = M_Sum/len(AdjMat_Files)
np.savetxt('/home/despoB/kaihwang/bin/FuncParcel/StriatalCorticalAveMat', AdjMat)


# load the averaged Striatalcortical adj matrix
#AdjMat = np.loadtxt('/Users/kaihwang/Google Drive/Projects/Striatal-Rest/StriatalcorticalAveMat')

# load the ROI name vector of the full Striatalcortical adj matrix
Striatalcortical_ROIs = np.loadtxt('/home/despoB/kaihwang/bin/FuncParcel/Striatalcortical_ROIs_index')

# load the name vector of the Striatal voxels
Striatal_voxels = np.loadtxt('/home/despoB/kaihwang/bin/FuncParcel/striatal_voxel_index')

# load the name vector of the cortical ROIs 
Cortical_ROIs = np.loadtxt('/home/despoB/kaihwang/bin/FuncParcel/Cortical_ROI_index')


# Double check the ROI numbers match
print(len(Striatalcortical_ROIs) == len(Cortical_ROIs) + len(Striatal_voxels))


# Each ROI and Striatal voxel has a unique integer. So we need to extract the matrix indeces of those numbers.

#create cortical ROI index
Cortical_indices = np.array([])
for i in range(0,len(Cortical_ROIs)):
    Cortical_indices = np.append(Cortical_indices, np.where(Striatalcortical_ROIs==Cortical_ROIs[i]))

Cortical_indices = Cortical_indices.astype(int)

#Create Striatal voxel index
Strital_indices = np.array([])
for i in range(0,len(Striatal_voxels)):
    Strital_indices = np.append(Strital_indices, np.where(Striatalcortical_ROIs==Striatal_voxels[i]))
    
Striatal_indices = Striatal_indices.astype(int)  # I don't understand why I can't change the datatype itself using copy=false


# Then for each Striatal voxel rank the cortical ROIs based on its Striatalcortical connectivity strength. 
#initialize output as dictionary. 
Striatalcortical_Targets = {}

for i_Striatal in Striatal_indices:
    r_vec = AdjMat[i_Striatal,Cortical_indices] #extract r values bewteen a given Striatal voxel and cortical ROIs
    index = np.argsort(r_vec) # sort r values, save index       
    index = index[::-1] # reverse index (highest to lowest r values)
    
    # rank the cortical ROI number by highest to lowerst r value, then save to a dictionary where the voxel number is the key 
    Striatalcortical_Targets[str(int(Striatalcortical_ROIs[i_Striatal]))] = Cortical_ROIs[index]



# load lesion voxel numbers from patients

Patients = ['b117', 'b122', 'b138', 'b143', 'b153']

for sub in Patients:
    fn = "/home/despoB/kaihwang/bin/FuncParcel/%s_lesioned_voxels" %sub
    lesioned_vox = np.loadtxt(fn)
    
    Cortical_Targets = np.array([])
    
    for vox in lesioned_vox:
        Cortical_Targets = np.append(Cortical_Targets,Striatalcortical_Targets[str(int(vox))][0])
        Cortical_NonTargets = np.append(Cortical_NonTargets,Striatalcortical_Targets[str(int(vox))][-1])
    
    
    Cortical_Targets = np.unique(Cortical_Targets)
    Cortical_NonTargets = np.unique(Cortical_NonTargets)
    #Cortical_NonTargets = Striatalcortical_Targets[str(int(vox))][-1*len(Cortical_Targets):]
    
    fn = "/home/despoB/kaihwang/bin/FuncParcel/%s_cortical_target" %sub
    np.savetxt(fn, Cortical_Targets.astype(int), fmt='%3.d')
    
    fn = "/home/despoB/kaihwang/bin/FuncParcel/%s_cortical_nontarget" %sub
    np.savetxt(fn, Cortical_NonTargets.astype(int), fmt='%3.d')






