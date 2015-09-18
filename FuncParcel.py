from __future__ import division, print_function
import numpy as np
import scipy.io as sio
from scipy import stats, linalg
#import matplotlib.pyplot as plt
import glob
import pandas as pd
import networkx as nx
from brainx import weighted_modularity
import bct
import nibabel as nib
from collections import Counter
import cPickle as pickle


def average_corrmat(file_path):
	'''
	This is a short function to average correlation coefficient matrices.
	
	------
    Parameters
    ------
    file_path:
		file_path = '/home/despoB/kaihwang/Rest/AdjMatrices/*Ses1_Full_WashU333_corrmat'

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
	#np.savetxt('/home/despoB/kaihwang/Rest/Striatum_parcel/StriatalCorticalAveMat', AdjMat)
	return AdjMat


def parcel_subcortical_network(path_to_adjmat = '/home/despoB/kaihwang/Rest/bin/FuncParcel/Data/ThalamoCorticalAveMat', \
	path_to_list_of_subcorticalcortical_ROIs = '/home/despoB/kaihwang/bin/FuncParcel/Data/Thalamocortical_ROIs_index', \
	path_to_list_of_subcortical_voxels = '/home/despoB/kaihwang/bin/FuncParcel/Data/Thalamic_voxel_index', \
	path_to_list_of_cortical_ROIs ='/home/despoB/kaihwang/bin/FuncParcel/Data/Cortical_ROI_index', \
	path_to_cortical_CI = '/home/despoB/kaihwang/bin/FuncParcel/Data/MGH_CI'):
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
		path_to_list_of_subcorticalcortical_ROIs = '/home/despoB/kaihwang/bin/FuncParcel/Data/Striatalcortical_ROIs_index'
	path_to_list_of_subcortical_voxels:	
		path_to_list_of_subcortical_voxels = '/home/despoB/kaihwang/bin/FuncParcel/Data/striatal_voxel_index'
	path_to_list_of_cortical_ROIs:
		path_to_list_of_cortical_ROIs = '/home/despoB/kaihwang/bin/FuncParcel/Data/Cortical_ROI_index'
	Cortical_CI:  path to where a list of cortical ROI parittion is saved	
		Cortical_CI = '/home/despoB/kaihwang/bin/FuncParcel/Data/MGH_CI'
    ------
    Return
    ------
    Subcorticalcortical_Targets: dict
    	function will return a dictionary "Subcorticalcortical_Targets". Where the key is the subcortical voxel index, and 
    	values are cortical ROIs ranked by strength of conenction with the subcortical voxel
    Subcortical_ParcelCI: 
    	 key is the subcortical voxel index, and values are arrays of cortical paritions ranked by strength of conenction with 
    	 the subcortical voxel
    Subcorticalcortical_Targets_corr: 
    	key is the subcortical voxel index, and values of arrays of correlatin coefficeint between subcortical voxel and cortical ROI,
    	ranked by the strength	 
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
    # load the CI vector
	Cortical_CI = np.loadtxt(path_to_cortical_CI)
	Cortical_CI = Cortical_CI.astype(int) # for some reason using dtype = int doesn't work...
	
	# Double check the ROI numbers match
	assert len(Subcorticalcortical_ROIs) == len(Cortical_ROIs) + len(Subcortical_voxels)

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
	Subcortical_ParcelCI = {}
	Subcorticalcortical_Targets_corr = {}
	for i_Subcortical in Subcortical_indices:
	    r_vec = AdjMat[i_Subcortical,Cortical_indices] #extract r values bewteen a given Subcortical voxel and cortical ROIs
	    index = np.argsort(r_vec) # sort r values, save index       
	    index = index[::-1] # reverse index (highest to lowest r values)
	    # rank the cortical ROI number by highest to lowerst r value, then save to a dictionary where the voxel number is the key 
	    Subcorticalcortical_Targets[int(Subcorticalcortical_ROIs[i_Subcortical])] = Cortical_ROIs[index]
	    Subcortical_ParcelCI[int(Subcorticalcortical_ROIs[i_Subcortical])] = Cortical_CI[index]
	    Subcorticalcortical_Targets_corr[int(Subcorticalcortical_ROIs[i_Subcortical])]  = r_vec[index]

	return Subcorticalcortical_Targets, Subcortical_ParcelCI, Subcorticalcortical_Targets_corr


def subcortical_patients_cortical_target(Subcorticalcortical_Targets, Patients, numROI):
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

	numROI: integer
		number of ROIs you want to associate with each lesioned voxel.	

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
	Patients_Cortical_Targets, Patients_Cortical_non Targets = subcortical_patients_cortical_target(Subcorticalcortical_Targets, Patients, 3)
	

	'''

	# write out list of target and nontarget ROIs
	Patients_Cortical_Targets = {}
	Patients_Cortical_NonTargets = {}
	for sub in Patients:
	    print(sub)
	    fn = "/home/despoB/kaihwang/bin/FuncParcel/Data/%s_lesioned_voxels" %sub # at some point need to update this...
	    lesioned_vox = np.loadtxt(fn, dtype = int)
	    Cortical_Targets = np.array([])
	    Cortical_NonTargets = np.array([])
	    for vox in lesioned_vox:
	        Cortical_Targets = np.append(Cortical_Targets,Subcorticalcortical_Targets[int(vox)][0:numROI])
	        Cortical_NonTargets = np.append(Cortical_NonTargets,Subcorticalcortical_Targets[int(vox)][::-1][0:numROI])

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



def convert_matlab_graph_str(graph_path, SubID, Cortical_ROIs):
	'''
	Function to load matlab graph.mat structure, convert its global graph metric and nodal properties into a panda's dataframe format
	
	------
    Parameters
    ------
    graph_path: path to where matlb graph output (in matlab structure) is saved.
    	example graph_path = '/home/despoB/kaihwang/Rest/Graph/gsetCI_176.mat'
    
    SubID: Subjet ID in string
    	example SubID = '176'

    Cortical_ROIs: list of array indexing the ROI number	
	
    ------
    Return
    ------
    GlobalDataframe: Returns a panda's dataframe object for global network properties
	NodalDataframe: Returns a panda's dataframe object for nodal graph metrics

	------
	Usage
    ------
    usage: GlobalDataframe, NodalDataframe = convert_matlab_graph_str(graph_path, SubID, Cortical_ROIs)
	
	'''
	mat_struct = sio.loadmat(graph_path)

	# check density vector is the same as what we think it is....
	density_vector = np.linspace(0.01, 0.25, num = 49)
	assert len(mat_struct['Graph']['Full_Q'][0,0][0,0][0]) == len(density_vector)
	
	# extract global graph metric from matlab struct, put in dictionary
	tmp_dict = {}
	tmp_dict['Density'] = density_vector
	tmp_dict['Q'] = mat_struct['Graph']['Full_Q'][0,0][0,0][0]
	tmp_dict['left_Q'] = mat_struct['Graph']['Left_Q'][0,0][0,0][0]
	tmp_dict['right_Q'] = mat_struct['Graph']['Right_Q'][0,0][0,0][0]
	tmp_dict['CC'] = np.mean(mat_struct['Graph']['Full_CC'][0,0][0,0],1)
	tmp_dict['left_CC'] = np.mean(mat_struct['Graph']['Left_CC'][0,0][0,0],1)
	tmp_dict['right_CC'] = np.mean(mat_struct['Graph']['Right_CC'][0,0][0,0],1)
	tmp_dict['Subject'] = [SubID] * len(density_vector)
	
	# conert dict to panda's dataframe
	GlobalDataframe = pd.DataFrame(tmp_dict, columns=['Subject', 'Density', 'Q', 'left_Q', 'right_Q','CC', 'left_CC','right_CC'])
	
	# extract nodal properties
	tmp_dict = {}

	tmp = mat_struct['Graph']['Full_CC'][0,0][0,0]
	tmp_dict['CC'] = tmp.ravel('F')

	tmp = mat_struct['Graph']['Within_Module_Weight'][0,0][0,0]
	tmp_dict['Within_Module_Weight'] = tmp.ravel('F')

	tmp = mat_struct['Graph']['Out_Module_Weight'][0,0][0,0]
	tmp_dict['Between_Module_Weight'] = tmp.ravel('F')

	tmp = mat_struct['Graph']['Within_Module_Degree'][0,0][0,0]
	tmp_dict['Within_Module_Degree'] = tmp.ravel('F')

	tmp = mat_struct['Graph']['Out_Module_Degree'][0,0][0,0]
	tmp_dict['Between_Module_Degree'] = tmp.ravel('F')

	tmp = mat_struct['Graph']['P'][0,0][0,0]
	tmp_dict['PC'] = tmp.ravel('F')

	tmp = mat_struct['Graph']['WMD'][0,0][0,0]
	tmp_dict['WMD'] = tmp.ravel('F')

	tmp = mat_struct['Graph']['Full_locE'][0,0][0,0]
	tmp_dict['localE'] = tmp.ravel('F')

	tmp = mat_struct['Graph']['Full_Ci'][0,0][0,0]
	tmp_dict['Ci'] = tmp.ravel('F')

	tmp_dict['Density'] = np.tile(density_vector,len(Cortical_ROIs))
	tmp_dict['ROI'] = np.repeat(Cortical_ROIs,len(density_vector))
	tmp_dict['Subject'] = [SubID] * len(tmp.ravel('F'))
	NodalDataframe = pd.DataFrame(tmp_dict, columns=['Subject', 'ROI', 'Density','CC', 'PC', 'WMD', 'Ci', 'localE', 'Between_Module_Weight','Within_Module_Weight', 'Between_Module_Degree','Within_Module_Degree'])

	return GlobalDataframe, NodalDataframe


def cal_graph_z_score(PatientDataframe, ControlDataframe, rois, metric):
	''' 
	Function to load pateint's graph metric, convert that metric into z-score (number of stds relative to control subjects).

	----
	Parameters
	----
	PatientDataframe: Dataframe of graph metrics from a particular patient. Should be output from 'convert_matlab_graph_str'

	ControlDataframe: Dataframe of graph metrics from control subjects. Note should be "ALL" control subjects, can use 'convert_matlab_graph_str' and append
		results from each subject

	rois: np array of list of ROIs, empty if for gobal data

	metric: string of the graph metric you want to convert to z-score

	----
	Return
	----
	z_score of the graph metric	
	'''
	#GlobalData, NodalData = convert_matlab_graph_str(patient_graph_path, patient_ID, Cortical_ROIs)
	if not rois.any():
		z_score = (PatientDataframe[metric] -ControlDataframe.groupby(['Density']).aggregate(np.nanmean).reset_index()[metric]) / ControlDataframe.groupby(['Density']).aggregate(np.nanstd).reset_index()[metric]
	else:
		#z_score = (PatientDataframe[metric] -ControlDataframe.groupby(['ROI','Density']).aggregate(np.nanmean).reset_index()[metric]) / ControlDataframe.groupby(['ROI','Density']).aggregate(np.nanstd).reset_index()[metric]
		x = PatientDataframe.loc[PatientDataframe['ROI'].isin(rois)].groupby(['Density']).aggregate(np.nanmean).reset_index()[metric]
		m = ControlDataframe.loc[ControlDataframe['ROI'].isin(rois)].groupby(['Density']).aggregate(np.nanmean).reset_index()[metric]
		sd = ControlDataframe.loc[ControlDataframe['ROI'].isin(rois)].groupby(['Density']).aggregate(np.nanstd).reset_index()[metric]
		z_score = (x-m) / sd
	return z_score	


def convert_partition_to_dict(input_partition):
    '''
    Function to convert graph partitin output from brainx's louvain function into a dictionary. 

    ----
    Parameters
    ----
    input_partition: partition output from LouvainCommunityDetection


    ----
    Return
    -----
    Parition in dicitionary, each key is a community, value is a list of ROIs that belong to this community
    '''
    partition_dict = {}
    Communities = 0
    for ROIs in input_partition:
        if len(list(ROIs)) != 1:
        	partition_dict[Communities] = list(ROIs)
        	Communities = Communities + 1
    return partition_dict



def convert_partition_dict_to_array(d, roi_num):
	''' To convert brainx's network partition dictionary into an array
	usage: ci = convert_partition_dict_to_array(d, roi_num)

	----
	Parameters
	----
	d: Dictionary of module partitions from brainx
	roi_num: int number indicating number of rois

	----
	Returns
	----
	ci : vector of community partition assignments for each ROI. 
	'''
	ci = np.zeros(roi_num)
	for key in d:
		for i in xrange(len(d[key])):
			ci[d[key][i]] = key

	return ci

def convert_graph_metric_dict_to_array(d, roi_num):
	''' To convert brainx's network graph metric dictionary into an array
	usage: output = convert_graph_metric_dict_to_array(d, roi_num)

	----
	Parameters
	----
	d: Dictionary of graph metrics from brainx
	roi_num: int number indicating number of rois

	----
	Returns
	----
	result : array
	'''

	output = np.zeros(roi_num)
	for key in d:
		output[key] = d[key]

	return output


def iterate_modularity_partition(subject, iter):
	''' Function to run iterations of Louvain modularity detection, return the one with the highest modularity (Q)
	usage: ci, q = iterate_modularity_partition(subject, iter)

	----
	Parameters
	----
	subject: subject number, will look for its correlation corrmat file with a hard coded path... (should prob change that)
	iter: int number indicating number of iterations

	----
	Returns
	----
	ci : community partition
	q  : modularity
	'''

	fn = '/home/despoB/kaihwang/Rest/AdjMatrices/t%s_Full_WashU333_corrmat' %subject
	AveMat = np.loadtxt(fn)
	graph = nx.from_numpy_matrix(bct.binarize(bct.threshold_proportional(AveMat, 0.05)))
	q = 0
	for i in xrange(0,iter):
		print(i)
		louvain = weighted_modularity.LouvainCommunityDetection(graph)
		weighted_partitions = louvain.run()
		if weighted_partitions[0].modularity() > q:
			q = weighted_partitions[0].modularity()
			weighted_partition = weighted_partitions[0]
			ci = convert_partition_dict_to_array(convert_partition_to_dict(weighted_partition.communities), len(AveMat))
	return ci, q

def make_image(atlas_path,image_path,ROI_list,values):
	''' Function to write out modularity community assignments to a nifit image
	usage: make_image(atlas_path,image_path,ROI_list,values)

	----
	Parameters
	----
	atlas_path : the ROI template used. Path to a nifit file
	image_path : the output path
	ROI_list : a text file of list of ROI indices in the ROI template
	values : a vector of community assignment, has to be the same length as ROI_list
	'''
	image = nib.load(atlas_path)
	image_data = image.get_data()
	ROIs = np.loadtxt(ROI_list, dtype = int)

	# check ROI number and CI are the same length
	assert len(ROIs) == len(values)

	value_data = image_data.copy()	
	partition_count = Counter(values)

	for ix,i in enumerate(values):

		if partition_count[i] > 1:
			value_data[image_data==ROIs[ix]] = i
		else: 
			value_data[image_data==ROIs[ix]] = 0

	image_data[:,:,:,] = value_data[:,:,:,]
	nib.save(image,image_path)


def save_object(obj, filename):
	''' Simple function to write out objects into a pickle file
	usage: save_object(obj, filename)
	'''
	with open(filename, 'wb') as output:
		pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)



def pcorr_subcortico_cortical_connectivity(subcortical_ts, cortical_ts):
	''' function to do partial correlation bewteen subcortical and cortical ROI timeseries. 
	Cortical signals (not subcortical) will be removed from subcortical and cortical ROIs,
	and then pairwise correlation will be calculated bewteen subcortical and cortical ROIs 
	(but not between subcortico-subcortical or cortico-cortical ROIs).
	This partial correlation/regression approach is for  cleaning subcortico-cortical 
	conectivity, which seems to be heavily influenced by a global noise.


	usage: pcorr_mat = pcorr_subcortico-cortical(subcortical_ts, cortical_ts)

	----
	Parameters
	----
	subcortical_ts: txt file of timeseries data from subcortical ROIs/voxels, each roi is an ROI
	cortical_ts: txt file of timeseries data from cortical ROIs, each roi is an ROI
	pcorr_mat: output partial correlation matrix
	'''

	# transpose so that column is ROI, this is because output from 3dNetcorr is row-based.
	subcortical_ts = subcortical_ts.T
	cortical_ts = cortical_ts.T

	#first check that the dimension is appropriate
	num_cort = cortical_ts.shape[1]
	num_subcor = subcortical_ts.shape[1]
	num_total = num_cort + num_subcor
	pcorr_mat = np.zeros((num_total, num_total), dtype=np.float)

	
	for j in range(num_cort):	
		k = np.ones(num_cort, dtype=np.bool)
		k[j] = False
		# fit cortical signal to cortical ROI TS, get betas
		beta_cortical = linalg.lstsq(cortical_ts[:,k], cortical_ts[:,j])[0]

		#get residuals
		res_cortical = cortical_ts[:, j] - cortical_ts[:, k].dot(beta_cortical)
		
		for i in range(num_subcor):
			# fit cortical signal to subcortical ROI TS, get betas
			beta_subcortical = linalg.lstsq(cortical_ts[:,k], subcortical_ts[:,i])[0]

			#get residuals
			res_subcortical = subcortical_ts[:, i] - cortical_ts[:, k].dot(beta_subcortical)

			#partial correlation
			pcorr_mat[i+num_cort, j] = stats.pearsonr(res_cortical, res_subcortical)[0]
			pcorr_mat[j,i+num_cort ] = pcorr_mat[i+num_cort, j]  

	return pcorr_mat		






