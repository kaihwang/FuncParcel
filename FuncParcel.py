from __future__ import division, print_function
import numpy as np
import scipy as sp
import scipy.io as sio
import matplotlib.pyplot as plt
import pickle
import csv
import glob
import pandas as pd



def average_corrmat(file_path):
	'''
	This is a short function to average correlation coefficient matrices.
	
	------
    Parameters
    ------
    file_path:
		file_path: = '/home/despoB/kaihwang/Rest/AdjMatrices/*Ses1_FIX_Striatalcortical_corrmat'

	usage: average_corrmat(file_path)

	------
    Return
    ------
	AdjMat: it will return an averaged adjmat
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
	'''#Here try to do subcortical parcellation in python
	The process works like this. 
	* You will need a set of ROI that consists of: (1) the cortex, (2) a subcortical mask that you wish to parcellate based 
	on its connectivity profile with the cortex.
	* The cortical ROIs should be labeled as integer in a NIFTI file. The subcortical mask should also have each voxel labeled 
	as different integer values. Note the values cannot overlap with corticla ROIs.
	* This combined ROI set will then be used to generate a subcortical+cortical full correlation matrix (adjacency matrix). 
	Essentially BOLD timesires from each ROI and subcortical voxel will be correlated with each other using 3dNetCorr.
	* The output adjancency matrix from each subject can then be averaged. Alternatively it is possible to do this parcellation 
	on a subject by subject basis, but its not implemented yet.
	* For each subcortical voxel, its connectivity strength with different cortical ROIs will then be sorted (ranked). 
	* You can then classify the subcortical voxel based which cortical region (network, module) its most strongly connected to.

	ROI used in the study
	* Total of 3836 ROIs in the Striatal + cortex parcellation file
	* 3539 voxels in thalammic parcellation files
	* 297 ROIs in cortical ROI file. 

	------
    Parameters
    ------
    path_to_adjmat:
		path_to_adjmat = '/home/despoB/kaihwang/Rest/Striatum_parcel/StriatalCorticalAveMat'
	path_to_list_of_subcorticalcortical_ROIs:
		path_to_list_of_subcorticalcortical_ROIs = '/home/despoB/kaihwang/bin/FuncParcel/Striatalcortical_ROIs_index'
	path_to_list_of_subcortical_voxels:	
		path_to_list_of_subcortical_voxels = '/home/despoB/kaihwang/bin/FuncParcel/striatal_voxel_index'
	path_to_list_of_cortical_ROIs:
		path_to_list_of_cortical_ROIs = '/home/despoB/kaihwang/bin/FuncParcel/Cortical_ROI_index'

    ------
    Return
    ------
    Subcorticalcortical_Targets: dict
    	function will return a dictionary "Subcorticalcortical_Targets". Where the key is the subcortical voxel index, and 
    	values are cortical ROIs ranked by strength of conenction with the subcortical voxel

    ------
    Example usage
    ------
	To run this program give the following parameters:
		parcel_subcortical_network(path_to_adjmat, path_to_list_of_subcorticalcortical_ROIs, ...
			path_to_list_of_subcortical_voxels, path_to_list_of_cortical_ROIs)
		

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


def subcortical_patients_cortical_target(Subcorticalcortical_Targets, Patients):
	'''
	function to identify lesioned subcortical voxels' cortical targets and non-targets
	targets defined as top connected ROI
	non-target defined as the least connected ROI
	
	------
    Parameters
    ------
    Subcorticalcortical_Targets: dict
    	Subcorticalcortical_Targets is a dictionary created by function parcel_subcortical_network

    Patients: vector	
		example patient vector: Patients = ['b117', 'b122', 'b138', 'b143', 'b153']

	------
	Return
	------
	Patients_Cortical_Targets: dict
		A dictionary where key is patient ID, values are the targeted ROI indices of lesioned voxel

	Patients_Cortical_NonTargets: dict
		A dictionary where key is paeitnt ID, values are non-targeted ROI indices

	------	
	Usage 
	------
	Patients_Cortical_Targets, Patients_Cortical_non Targets = subcortical_patients_cortical_target(Subcorticalcortical_Targets, Patients)
	

	'''

	# write out list of target and nontarget ROIs
	Patients_Cortical_Targets = {}
	Patients_Cortical_NonTargets = {}
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
	    Patients_Cortical_Targets[sub] = Cortical_Targets
	    Patients_Cortical_NonTargets[sub] = Cortical_NonTargets

	return Patients_Cortical_Targets, Patients_Cortical_NonTargets     
	    #Cortical_NonTargets = Subcorticalcortical_Targets[str(int(vox))][-1*len(Cortical_Targets):]
	    # save data
	    #fn = "/home/despoB/kaihwang/bin/FuncParcel/%s_cortical_target" %sub
	    #np.savetxt(fn, Cortical_Targets.astype(int), fmt='%3.d')
	    #fn = "/home/despoB/kaihwang/bin/FuncParcel/%s_cortical_nontarget" %sub
	    #np.savetxt(fn, Cortical_NonTargets.astype(int), fmt='%3.d')



def convert_matlab_graph_str(graph_path, SubID):
	'''
	Function to load matlab graph.mat structure, convert its global graph metric content into a panda's dataframe 
	
	------
    Parameters
    ------
    graph_path: path to where matlb graph output (in matlab structure) is saved.
    	example graph_path = '/home/despoB/kaihwang/Rest/Graph/gsetCI_176.mat'
    
    SubID: Subjet ID in string
    	example SubID = '176'
	
    ------
    Return
    ------
    GlobalDataframe: Returns a panda's dataframe object for global network properties
	NodalDataframe: Returns a panda's dataframe object for nodal graph metrics

	------
	Usage
    ------
    usage: GlobalDataframe, NodalDataframe = convert_matlab_graph_str(graph_path)
	
	'''
	mat_struct = sp.io.loadmat(graph_path)

	# check density vector 
	density_vector = np.arange(0.01, 0.255, 0.005)
	len(mat_struct['Graph']['Full_Q'][0,0][0,0][0]) == len(density_vector)

	rep_num = len(density_vector)
	
	# extract global graph metric from matlab struct, put in dictionary
	tmp_dict = {}
	tmp_dict['Density'] = density_vector
	tmp_dict['Q'] = mat_struct['Graph']['Full_Q'][0,0][0,0][0]
	tmp_dict['left_Q'] = mat_struct['Graph']['Left_Q'][0,0][0,0][0]
	tmp_dict['right_Q'] = mat_struct['Graph']['Right_Q'][0,0][0,0][0]
	tmp_dict['CC'] = np.mean(mat_struct['Graph']['Full_CC'][0,0][0,0],1)
	tmp_dict['left_CC'] = np.mean(mat_struct['Graph']['Left_CC'][0,0][0,0],1)
	tmp_dict['right_CC'] = np.mean(mat_struct['Graph']['Right_CC'][0,0][0,0],1)
	tmp_dict['Subject'] = [SubID] * rep_num
	
	# conert dict to panda's dataframe
	GlobalDataframe = pd.DataFrame(tmp_dict, columns=['Subject', 'Density', 'Q', 'left_Q', 'right_Q','CC', 'left_CC','right_CC'])
	
	# extract nodal properties
	tmp_dict = {}
	tmp = mat_struct['Graph']['Full_CC'][0,0][0,0]
	tmp_dict['CC'] = tmp.ravel()
	tmp_dict['Density'] = np.tile(density_vector,len(Cortical_ROIs))
	tmp_dict['ROI'] = np.repeat(Cortical_ROIs,len(density_vector))
	NodalDataframe = pd.DataFrame(tmp_dict, columns=['ROI', 'Density','CC'])
	#tmp_dict['Density'] = density_vector
	#tmp_dict['CC'] = mat_struct['Graph']['Full_CC'][0,0][0,0][:,1]
	#pd.DataFrame(tmp_dict, columns=['Subject', 'Density', 'Q', 'left_Q', 'right_Q','CC', 'left_CC','right_CC'])
	#return
	return GlobalDataframe NodalDataframe






