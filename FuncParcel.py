from __future__ import division, print_function
import numpy as np
import scipy as sp
import scipy.io as sio
import matplotlib.pyplot as plt
import pickle
import csv
import glob



def average_corrmat(file_path):
	'''
	This is a short function to average correlation coefficient matrices.
	example file_path = '/home/despoB/kaihwang/Rest/AdjMatrices/*Ses1_FIX_Striatalcortical_corrmat'
	'''
	# Averaging
	AdjMat_Files = glob.glob(file_path)

	M_Sum =[];
	for f in AdjMat_Files:
	    M = np.loadtxt(f)
	    M_Sum += [M]
	    
	    
	AdjMat = sum(M_Sum)/len(AdjMat_Files)
	np.savetxt('/home/despoB/kaihwang/Rest/Striatum_parcel/StriatalCorticalAveMat', AdjMat)
	return AdjMat


def parcel_subcortical_network(path_to_adjmat, path_to_list_of_subcorticalcortical_ROIs, path_to_list_of_subcortical_voxels, path_to_list_of_cortical_ROIs):
	'''#Here try to do Striatalcortical parcellation in python

	The process works like this. 

	* You will need a set of ROI that consists of: (1) the cortex, (2) a subcortical mask that you wish to parcellate based on its connectivity profile with the cortex, in this case the Striatal. 
	* The cortical ROIs should be labeled as integer in a NIFTI file. The subcortical mask should also have each voxel labeled as different integer values. Note the values cannot overlap with corticla ROIs.
	* This combined ROI set will then be used to generate a subcortical+cortical full correlation matrix (adjacency matrix). Essentially BOLD timesires from each ROI and subcortical voxel will be correlated with each other using 3dNetCorr.
	* The output adjancency matrix from each subject can then be averaged. Alternatively it is possible to do this parcellation on a subject by subject basis, but its not implemented yet.
	* For each subcortical voxel, its connectivity strength with different cortical ROIs will then be sorted (ranked). 
	* You can then classify the subcortical voxel based which cortical region (network, module) its most strongly connected to.

	ROI used in the following example
	* Total of 3836 ROIs in the Striatal + cortex parcellation file
	* 3539 voxels in thalammic parcellation files
	* 297 ROIs in cortical ROI file. 

	Example paths
		path_to_adjmat = '/home/despoB/kaihwang/Rest/Striatum_parcel/StriatalCorticalAveMat'
		path_to_list_of_subcorticalcortical_ROIs = '/home/despoB/kaihwang/bin/FuncParcel/Striatalcortical_ROIs_index'
		path_to_list_of_subcortical_voxels = '/home/despoB/kaihwang/bin/FuncParcel/striatal_voxel_index'
		path_to_list_of_cortical_ROIs = '/home/despoB/kaihwang/bin/FuncParcel/Cortical_ROI_index'
	'''



	# load the averaged adj matrix
	AdjMat = np.loadtxt(path_to_adjmat)
	# load the ROI name vector of the full Subcorticalcortical adj matrix
	Subcorticalcortical_ROIs = np.loadtxt(path_to_list_of_subcorticalcortical_ROIs)
	# load the name vector of the Subcortical voxels
	Subcortical_voxels = np.loadtxt(path_to_list_of_subcortical_voxels)
	# load the name vector of the cortical ROIs 
	Cortical_ROIs = np.loadtxt(path_to_list_of_cortical_ROIs)

	# Double check the ROI numbers match
	print(len(Subcorticalcortical_ROIs) == len(Cortical_ROIs) + len(Subcortical_voxels))

	# Each ROI and Subcortical voxel has a unique integer. So we need to extract the matrix indeces of those numbers.
	#create cortical ROI index
	Cortical_indices = np.array([])
	for i in range(0,len(Cortical_ROIs)):
	    Cortical_indices = np.append(Cortical_indices, np.where(Subcorticalcortical_ROIs==Cortical_ROIs[i]))

	Cortical_indices = Cortical_indices.astype(int)
	#Create Subcortical voxel index
	Subcortical_indices = np.array([])
	for i in range(0,len(Subcortical_voxels)):
	    Subcortical_indices = np.append(Subcortical_indices, np.where(Subcorticalcortical_ROIs==Subcortical_voxels[i]))
	    
	Subcortical_indices = Subcortical_indices.astype(int)  # I don't understand why I can't change the datatype itself using copy=false

	# Then for each Subcortical voxel rank the cortical ROIs based on its Subcorticalcortical connectivity strength. 
	#initialize output as dictionary. 
	Subcorticalcortical_Targets = {}
	for i_Subcortical in Subcortical_indices:
	    r_vec = AdjMat[i_Subcortical,Cortical_indices] #extract r values bewteen a given Subcortical voxel and cortical ROIs
	    index = np.argsort(r_vec) # sort r values, save index       
	    index = index[::-1] # reverse index (highest to lowest r values)
	    # rank the cortical ROI number by highest to lowerst r value, then save to a dictionary where the voxel number is the key 
	    Subcorticalcortical_Targets[str(int(Subcorticalcortical_ROIs[i_Subcortical]))] = Cortical_ROIs[index]

	return Subcorticalcortical_Targets


def subcortical_patients_cortical_target_change(Subcorticalcortical_Targets, Patients):
	'''
	example patient vector: Patients = ['b117', 'b122', 'b138', 'b143', 'b153']
	'''

	# write out list of target and nontarget ROIs
	for sub in Patients:
	    print(sub)
	    fn = "/home/despoB/kaihwang/bin/FuncParcel/%s_lesioned_voxels" %sub
	    lesioned_vox = np.loadtxt(fn)
	    Cortical_Targets = np.array([])
	    Cortical_NonTargets = np.array([])
	    for vox in lesioned_vox:
	        Cortical_Targets = np.append(Cortical_Targets,Subcorticalcortical_Targets[str(int(vox))][0])
	        Cortical_NonTargets = np.append(Cortical_NonTargets,Subcorticalcortical_Targets[str(int(vox))][-1])

	    # take out repetitions
	    Cortical_Targets = np.unique(Cortical_Targets)
	    Cortical_NonTargets = np.unique(Cortical_NonTargets)
	    #Cortical_NonTargets = Subcorticalcortical_Targets[str(int(vox))][-1*len(Cortical_Targets):]
	    # save data
	    fn = "/home/despoB/kaihwang/bin/FuncParcel/%s_cortical_target" %sub
	    np.savetxt(fn, Cortical_Targets.astype(int), fmt='%3.d')
	    fn = "/home/despoB/kaihwang/bin/FuncParcel/%s_cortical_nontarget" %sub
	    np.savetxt(fn, Cortical_NonTargets.astype(int), fmt='%3.d')






