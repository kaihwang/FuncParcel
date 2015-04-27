# # Here try to do thalamocortical parcellation in python
# 
# The process works like this. 
# 
# * You will need a set of ROI that consists of: (1) the cortex, (2) a subcortical mask that you wish to parcellate based on its connectivity profile with the cortex, in this case the thalamus. 
# * The cortical ROIs should be labeled as integer in a NIFTI file. The subcortical mask should also have each voxel labeled as different integer values. Note the values cannot overlap with corticla ROIs.
# * This combined ROI set will then be used to generate a subcortical+cortical full correlation matrix (adjacency matrix). Essentially BOLD timesires from each ROI and subcortical voxel will be correlated with each other using 3dNetCorr.
# * The output adjancency matrix from each subject can then be averaged. Alternatively it is possible to do this parcellation on a subject by subject basis, but its not implemented yet.
# * For each subcortical voxel, its connectivity strength with different cortical ROIs will then be sorted (ranked). 
# * You can then classify the subcortical voxel based which cortical region (network, module) its most strongly connected to.

# ROI used in the following example
# * Total of 3836 ROIs in the Thalamus + cortex parcellation file
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


# In[ ]:




# In[3]:

# load the averaged thalamocortical adj matrix
#AdjMat = np.loadtxt('/Users/kaihwang/Google Drive/Projects/Thalamus-Rest/ThalamoCorticalAveMat')
#(I'm lazy, that averaged matrices was created in matlab.... will write a separate block code sometime)

# load the ROI name vector of the full thalamocortical adj matrix
Thalamocortical_ROIs = np.loadtxt('./Thalamocortical_ROIs_index')

# load the name vector of the thalamic voxels
Thalamic_voxels = np.loadtxt('./Thalamic_voxel_index')

# load the name vector of the cortical ROIs 
Cortical_ROIs = np.loadtxt('./Cortical_ROI_index')


# In[4]:

# Double check the ROI numbers match
print(len(Thalamocortical_ROIs) == len(Cortical_ROIs) + len(Thalamic_voxels))


# Each ROI and thalamic voxel has a unique integer. So we need to extract the matrix indeces of those numbers.

# In[5]:

#create cortical ROI index
Cortical_indices = np.array([])
for i in range(0,len(Cortical_ROIs)):
    Cortical_indices = np.append(Cortical_indices, np.where(Thalamocortical_ROIs==Cortical_ROIs[i]))

Cortical_indices = Cortical_indices.astype(int)

#Create thalamus voxel index
Thalamus_indices = np.array([])
for i in range(0,len(Thalamic_voxels)):
    Thalamus_indices = np.append(Thalamus_indices, np.where(Thalamocortical_ROIs==Thalamic_voxels[i]))
    
Thalamus_indices = Thalamus_indices.astype(int)  # I don't understand why I can't change the datatype itself using copy=false


# Then for each thalamus voxel rank the cortical ROIs based on its thalamocortical connectivity strength. 

# In[6]:

#initialize output as dictionary. 
Thalamocortical_Targets = {}

for i_thalamus in Thalamus_indices:
    r_vec = AdjMat[i_thalamus,Cortical_indices] #extract r values bewteen a given thalamic voxel and cortical ROIs
    index = np.argsort(r_vec) # sort r values, save index       
    index = index[::-1] # reverse index (highest to lowest r values)
    
    # rank the cortical ROI number by highest to lowerst r value, then save to a dictionary where the voxel number is the key 
    Thalamocortical_Targets[str(int(Thalamocortical_ROIs[i_thalamus]))] = Cortical_ROIs[index]


# Write out the target dictionary to a csv file and pickle

# In[7]:

with open('/Users/kaihwang/Google Drive/Projects/Thalamus-Rest/thalamocortical_target.csv', 'wb') as f: 
    w = csv.DictWriter(f, Thalamocortical_Targets.keys())
    w.writeheader()
    w.writerow(Thalamocortical_Targets)
    
pickle.dump( Thalamocortical_Targets, open( "/Users/kaihwang/Google Drive/Projects/Thalamus-Rest/Thalamocortical_Targets.p", "wb" ) )


# Read the pickle output

# In[8]:

Thalamocortical_Targets = pickle.load( open( "/Users/kaihwang/Google Drive/Projects/Thalamus-Rest/Thalamocortical_Targets.p", "rb" ) )


# In[9]:

Thalamocortical_Targets


# In[10]:

# write to text

def saveDict(fn,dict_rap):
    f=open(fn, "wb")
    w = csv.writer(f)
    for key, val in dict_rap.items():
        w.writerow([key, val])
    f.close()


# In[11]:

saveDict('/Users/kaihwang/Google Drive/Projects/Thalamus-Rest/thalamocortical_target.test', Thalamocortical_Targets)


# In[82]:

# load lesion voxel numbers from patients

Patients = [128, 162, 163, 168, 176]

for sub in Patients:
    fn = "/Users/kaihwang/bin/FuncParcel/%s_lesioned_voxels" %sub
    lesioned_vox = np.loadtxt(fn)
    
    Cortical_Targets = np.array([])
    
    for vox in lesioned_vox:
        Cortical_Targets = np.append(Cortical_Targets,Thalamocortical_Targets[str(int(vox))][0])
        Cortical_NonTargets = np.append(Cortical_NonTargets,Thalamocortical_Targets[str(int(vox))][-1])
    
    
    Cortical_Targets = np.unique(Cortical_Targets)
    Cortical_NonTargets = np.unique(Cortical_NonTargets)
    #Cortical_NonTargets = Thalamocortical_Targets[str(int(vox))][-1*len(Cortical_Targets):]
    
    fn = "/Users/kaihwang/bin/FuncParcel/%s_cortical_target" %sub
    np.savetxt(fn, Cortical_Targets.astype(int), fmt='%3.d')
    
    fn = "/Users/kaihwang/bin/FuncParcel/%s_cortical_nontarget" %sub
    np.savetxt(fn, Cortical_NonTargets.astype(int), fmt='%3.d')






