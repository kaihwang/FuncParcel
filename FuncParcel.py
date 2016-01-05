from __future__ import division, print_function
import numpy as np
import scipy.io as sio
from scipy import stats, linalg
#import matplotlib.pyplot as plt
import glob
import pandas as pd
import networkx as nx
#from brainx import weighted_modularity
import bct
import nibabel as nib
from collections import Counter
import pickle
from sklearn.decomposition import PCA


def average_corrmat(file_path, np_txt=False, pickle_object=True):
	'''
	This is a short function to average correlation coefficient matrices.
	
	------
    Parameters
    ------
    file_path:
		file_path = '/home/despoB/kaihwang/Rest/AdjMatrices/*Ses1_Full_WashU333_corrmat'

	usage: average_corrmat(file_path, np_txt=False, pickle_object=True)
	if the saved matrices are txt files from np.savetxt, set np_txt =True
	if in pickle object format, set pckle_object=True (default)
	------
    Return
    ------
	AdjMat: it will return an averaged adjmat
	'''
	# Averaging
	AdjMat_Files = glob.glob(file_path)

	M_Sum =[]
	for f in AdjMat_Files:
		if np_txt:
			M = np.loadtxt(f)
		if pickle_object:
			M = pickle.load(open(f, "rb"))
		M_Sum += [M]

	AdjMat = sum(M_Sum)/len(AdjMat_Files)
	#np.savetxt('/home/despoB/kaihwang/Rest/Striatum_parcel/StriatalCorticalAveMat', AdjMat)
	return AdjMat


def parcel_subcortical_network(AdjMat, Subcorticalcortical_ROIs, Subcortical_voxels, Cortical_ROIs, Cortical_CI):
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
	#AdjMat = np.loadtxt(path_to_adjmat)
	# load the ROI name vector of the full Subcorticalcortical adj matrix
	#Subcorticalcortical_ROIs = np.loadtxt(path_to_list_of_subcorticalcortical_ROIs)
	# load the name vector of the Subcortical voxels
	#Subcortical_voxels = np.loadtxt(path_to_list_of_subcortical_voxels)
	# load the name vector of the cortical ROIs 
	#Cortical_ROIs = np.loadtxt(path_to_list_of_cortical_ROIs)
    # load the CI vector
	#Cortical_CI = np.loadtxt(path_to_cortical_CI)
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
	NodalDataframe = pd.DataFrame(tmp_dict, columns=['Subject', 'ROI', 'Density','CC', 'PC', 'WMD', 'Ci', 'localE', 
		'Between_Module_Weight','Within_Module_Weight', 'Between_Module_Degree','Within_Module_Degree'])

	return GlobalDataframe, NodalDataframe


def cal_graph_z_score(PatientDataframe, ControlDataframe, rois, metric):
	''' 
	Function to load pateint's graph metric, convert that metric into z-score (number of stds relative to control subjects).

	----
	Parameters
	----
	PatientDataframe: Dataframe of graph metrics from a particular patient. Should be output from 'convert_matlab_graph_str'

	ControlDataframe: Dataframe of graph metrics from control subjects. Note should be "ALL" control subjects, 
	can use 'convert_matlab_graph_str' and append
		results from each subject

	rois: np array of list of ROIs, empty if for global data

	metric: string of the graph metric you want to convert to z-score

	----
	Return
	----
	z_score of the graph metric	
	'''
	#GlobalData, NodalData = convert_matlab_graph_str(patient_graph_path, patient_ID, Cortical_ROIs)
	if not rois.any():
		z_score = (PatientDataframe[metric] -ControlDataframe.groupby(['Density']).aggregate(np.nanmean).reset_index()
			[metric]) / ControlDataframe.groupby(['Density']).aggregate(np.nanstd).reset_index()[metric]
	else:
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
	header = image.get_header()
	header.set_data_dtype(np.float)

	# check ROI number and CI are the same length
	assert len(ROIs) == len(values)

	value_data = image_data.copy()	
	#partition_count = Counter(values)

	for ix,i in enumerate(values):
		value_data[image_data==ROIs[ix]] = i
		#if partition_count[i] > 1:
		#	value_data[image_data==ROIs[ix]] = i
		#else: 
	    #		value_data[image_data==ROIs[ix]] = 0

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

	# check length of data
	assert cortical_ts.shape[0] == subcortical_ts.shape[0]
	num_vol = cortical_ts.shape[0]

	#first check that the dimension is appropriate
	num_cort = cortical_ts.shape[1]
	num_subcor = subcortical_ts.shape[1]
	num_total = num_cort + num_subcor

	#maximum number of regressors that we can use
	max_num_components = int(num_vol/20)
	if max_num_components > num_cort:
		max_num_components = num_cort-1 

	pcorr_mat = np.zeros((num_total, num_total), dtype=np.float)

	for j in range(num_cort):	
		k = np.ones(num_cort, dtype=np.bool)
		k[j] = False

		#use PCA to reduce cortical data dimensionality
		pca = PCA(n_components=max_num_components)
		pca.fit(cortical_ts[:,k])
		reduced_cortical_ts = pca.fit_transform(cortical_ts[:,k])
		
		#print("Amount of varaince explanined after PCA: %s" %np.sum(pca.explained_variance_ratio_))  
		
		# fit cortical signal to cortical ROI TS, get betas
		beta_cortical = linalg.lstsq(reduced_cortical_ts, cortical_ts[:,j])[0]

		#get residuals
		res_cortical = cortical_ts[:, j] - reduced_cortical_ts.dot(beta_cortical)
		
		for i in range(num_subcor):
			# fit cortical signal to subcortical ROI TS, get betas
			beta_subcortical = linalg.lstsq(reduced_cortical_ts, subcortical_ts[:,i])[0]

			#get residuals
			res_subcortical = subcortical_ts[:, i] - reduced_cortical_ts.dot(beta_subcortical)

			#partial correlation
			pcorr_mat[i+num_cort, j] = stats.pearsonr(res_cortical, res_subcortical)[0]
			pcorr_mat[j,i+num_cort ] = pcorr_mat[i+num_cort, j]  

	return pcorr_mat


def par_pcorr_subcortico_cortical_connectivity(idx, subcortical_ts, cortical_ts):
	''' function to do partial correlation bewteen subcortical and cortical ROI timeseries. 
	Cortical signals (not subcortical) will be removed from subcortical and cortical ROIs,
	and then pairwise correlation will be calculated bewteen subcortical and cortical ROIs 
	(but not between subcortico-subcortical or cortico-cortical ROIs).
	This partial correlation/regression approach is for  cleaning subcortico-cortical 
	conectivity, which seems to be heavily influenced by a global noise.

	THIS VERSION IS CREATED TO DO PARALLEL FOR LOOPS


	usage: pcorr_coef = pcorr_subcortico-cortical(idx, subcortical_ts, cortical_ts)

	----
	Parameters
	----
	subcortical_ts: txt file of timeseries data from subcortical ROIs/voxels, each roi is an ROI
	cortical_ts: txt file of timeseries data from cortical ROIs, each roi is an ROI
	idx a tuple of indices like this (i,j):
		i = index for subcortical
		j = index for cortical
	pcorr_coef: return partial correlation between subcortical_ts(i) and cortical_ts(j)
	'''

	# transpose so that column is ROI, this is because output from 3dNetcorr is row-based.
	i = idx[0]
	j = idx[1]

	#subcortical_ts = np.loadtxt(subcortical_ts)
	#cortical_ts = np.loadtxt(cortical_ts)

	s_ts = subcortical_ts.T
	c_ts = cortical_ts.T

	#first check number of volumnes/timepoints/observations are the same
	assert c_ts.shape[0] == s_ts.shape[0]
	num_vol = c_ts.shape[0]

	#first check that the dintimension is appropriate
	num_cort = c_ts.shape[1]
	num_subcor = s_ts.shape[1]
	num_total = num_cort + num_subcor

	#maximum number of regressors that we can use
	max_num_components = int(num_vol/20)
	if max_num_components > num_cort:
			max_num_components = num_cort-1
	
	#use PCA to reduce cortical timeseries data dimension
	k = np.ones(num_cort, dtype=np.bool)
	k[j] = False

	pca = PCA(n_components=max_num_components)
	pca.fit(c_ts[:,k])
	reduced_c_ts = pca.fit_transform(c_ts[:,k])
	#print("Amount of varaince explanined by after PCA: %s" %np.sum(pca.explained_variance_ratio_))  

	# fit cortical signal to cortical ROI TS, get betas
	beta_cortical = linalg.lstsq(reduced_c_ts, c_ts[:,j])[0]

	#get residuals
	res_cortical = c_ts[:, j] - reduced_c_ts.dot(beta_cortical)
		
	# fit cortical signal to subcortical ROI TS, get betas
	beta_subcortical = linalg.lstsq(reduced_c_ts, s_ts[:,i])[0]

	#get residuals
	res_subcortical = s_ts[:, i] - reduced_c_ts.dot(beta_subcortical)

	#partial correlation
	pcorr_coef = stats.pearsonr(res_cortical, res_subcortical)[0]

	return pcorr_coef	


def map_subcortical_cortical_targets(corrmat, Cortical_ROIs, Subcortical_voxels):
	''' function to generate a dictionary listing the cortical ROIs connected with each subcortical voxel, at different costs

	usage: cortical_targets, thresholds = map_subcortical_cortical_targets(corrmat, Cortical_ROIs, Subcortical_voxels):
	----
	Parameters
	----
	corrmat: adj matrix in numpy format, the cortical ROIs first, then concatnated by subcortical voxels
	Cortical_ROIs : list of cortical ROIs
	Subcortical_voxels : list of subcortical voxels

	It will return a dictionary of cortical_targets
	where the key is (cost, voxel_number), and the value is a list of cortical ROIs

	It will also return the threshold for each cost, from .15 to .01
	'''

	#do the annoying conversion to integer
	Subcortical_voxels = Subcortical_voxels.astype('int', copy=True)
	Cortical_ROIs = Cortical_ROIs.astype('int', copy=True)

	#check data length
	assert corrmat.shape[0] == Cortical_ROIs.shape[0]+ Subcortical_voxels.shape[0]
	assert corrmat.shape[1] == Cortical_ROIs.shape[0]+ Subcortical_voxels.shape[0]

	#determine size
	num_cortical_rois =  len(Cortical_ROIs)
	num_subcortical_voxels =  len(Subcortical_voxels)

	#extract subcortical_cortical_matrix, ignoring whithn cortical and within subcortical weights
	sc_mat = corrmat[range(num_cortical_rois, num_subcortical_voxels+num_cortical_rois),:][:, range(num_cortical_rois)].copy()

	#determine density thresholds, cost at .15 to .01, right now this range is hard coded
	thresholds = np.percentile(sc_mat.flatten(), range(85,100))
	#thresholds_lb = np.percentile(sc_mat.flatten(), range(1,16)) #lower bound threshold for neg weights

	costs = 100 -np.array(range(85,100))
	#save target in a dictionary, stepping by cost
	cortical_targets={}
	for p, th in enumerate(thresholds):
		for i in range(sc_mat.shape[0]):
			cortical_targets[costs[p],Subcortical_voxels[i]] = Cortical_ROIs[sc_mat[i]>th]
	
	cortical_nontargets={}
	for p, th in enumerate(thresholds):
		for i in range(sc_mat.shape[0]):
			rvec=abs(sc_mat[i]-0) #find the ones close to zero
			cortical_nontargets[costs[p],Subcortical_voxels[i]] \
			= Cortical_ROIs[rvec.argsort()[:len(cortical_targets[costs[p],Subcortical_voxels[i]])]]
			
	return cortical_targets, cortical_nontargets, thresholds


def pcorr_cortico_cortical_connectivity(cortical_ts):
    ''' function to do partial correlation within different ROIs of a cortical ROI
    timeseries. Based off of pcorr_subcortico_cortical_connectivity

    usage: pcorr_mat = pcorr_subcortico-cortical(cortical_ts)
    
    ----
    Parameters
    ----
    cortical_ts: txt file of timeseries data from cortical ROIs, each roi is an ROI
    pcorr_mat: output partial correlation matrix
    '''

    # transpose so that column is ROI, this is because output from 3dNetcorr is row-based.
    # subcortical_ts = subcortical_ts.T
    cortical_ts = cortical_ts.T

    # check length of data
    # assert cortical_ts.shape[0] == subcortical_ts.shape[0]
    num_vol = cortical_ts.shape[0]

    #first check that the dimension is appropriate
    num_cort = cortical_ts.shape[1]
    # num_subcor = subcortical_ts.shape[1]
    # num_cort = num_cort + num_subcor

    #maximum number of regressors that we can use
    max_num_components = int(num_vol/20)
    if max_num_components > num_cort:
        max_num_components = num_cort-1 

    pcorr_mat = np.zeros((num_cort, num_cort), dtype=np.float)
    for i in range(num_cort):
        pcorr_mat[i][i] = 1.0

    for j in range(num_cort):
        for i in range(j+1, num_cort):
            k = np.ones(num_cort, dtype=np.bool)
            k[j] = False
            k[i] = False

            #use PCA to reduce cortical data dimensionality
            pca = PCA(n_components=max_num_components)
            pca.fit(cortical_ts[:,k])
            reduced_cortical_ts = pca.fit_transform(cortical_ts[:,k])
            
            #print("Amount of varaince explanined after PCA: %s" %np.sum(pca.explained_variance_ratio_))  
            
            # fit cortical signal to cortical ROI TS, get betas
            beta_cortical_j = linalg.lstsq(reduced_cortical_ts, cortical_ts[:,j])[0]
            beta_cortical_i = linalg.lstsq(reduced_cortical_ts, cortical_ts[:,i])[0]

            #get residuals
            res_cortical_j = cortical_ts[:, j] - reduced_cortical_ts.dot(beta_cortical_j)
            res_cortical_i = cortical_ts[:, i] - reduced_cortical_ts.dot(beta_cortical_i)
            pcorr_mat[i, j] = stats.pearsonr(res_cortical_i, res_cortical_j)[0]
            pcorr_mat[j,i] = stats.pearsonr(res_cortical_i, res_cortical_j)[0]

    return pcorr_mat	


def cal_thalamus_and_cortical_ROIs_nodal_properties(Thalamocortical_corrmat, Cortical_adj, \
	Cortical_plus_thalamus_CI, Thalamus_CIs, Cortical_CI, Cortical_ROIs_positions, Thalamus_voxel_positions, cost_thresholds):
	'''Function to calculate voxel-wise nodal properties of the thalamus, and nodal properties of cortical ROIs. 
	Metrics to be calculated include:
	
	Participation Coefficient (PC)
	Between network connectivity weiight (BNWR)
		Ratio of connection weight devoted to between network interactions
	Number of network/modules/components connected (NNC)
	Within module degree zscore (WMD)
		For WMD, matrices will be binarzied, and normalized to corticocortical connections' mean and SD

	usage: PC, BNWR, NNC, WMD = cal_thalamus_and_cortical_ROIs_nodal_properties(Thalamocor_adj,
                Cortical_adj,
                Cortical_plus_thalamus_CI,
                Thalamus_CIs,
                Cortical_CI,
                Cortical_ROIs_positions,
                Thalamus_voxel_positions,
                cost_thresholds)
    
    ----
    Parameters
    ----
    Thalamocor_adj: Thalamocortical adj matrix
    Cortical_adj: corticocortical adj matrix
    Cortical_plus_thalamus_CI: A vector of community/module/network assignment of all nodes, cortical ROIs + thalamic voxels
    Thalamus_CIs: A vector of network assignements for thalamic voxels
    Cortical_CI: A vector of network assignments for cortical ROIs 
    Cortical_ROIs_positions: a position vector indicating in the thalamocortical adj matrix which rows/columns are cortical ROIs
    Thalamus_voxel_posistions: a position vector indicating in the thalamocortical adj matrix which rows/columns are thalamic voxels
    cost_thresholds: the thoresholds that can threshold the thalamocortical edges at density .01 to .15. For now this is hard coded
    '''

	##Thalamus nodal roles
	Thalamocortical_corrmat[np.isnan(Thalamocortical_corrmat)] = 0
	#PC
	PCs = np.zeros(Cortical_plus_thalamus_CI.size)
	#BNWR between network connectivity weight
	BNWRs = np.zeros(Cortical_plus_thalamus_CI.size)
	#get number of networks/communities connected
	NNCs = np.zeros(Cortical_plus_thalamus_CI.size)
		
	#loop through costs
	for c in cost_thresholds:
		#copy adj matrix and then threshold
		Par_adj = Thalamocortical_corrmat.copy()
		#remove weights connected to low SNR communities (CI==0, orbital frontal, inferior temporal)
		Par_adj[Cortical_ROIs_positions[Cortical_CI==0],:]=0
		Par_adj[:,Cortical_ROIs_positions[Cortical_CI==0]]=0
		Par_adj[Par_adj<c]=0

		#PC
		PCs += bct.participation_coef(Par_adj, Cortical_plus_thalamus_CI)

		#BNWR and NNCs
		Tha_BNWR = np.zeros(Cortical_plus_thalamus_CI.size)
		Tha_NNCs = np.zeros(Cortical_plus_thalamus_CI.size)
		for ix, i in enumerate(Thalamus_voxel_positions):
			sum_between_weight = np.nansum(Par_adj[i,Cortical_plus_thalamus_CI!=Thalamus_CIs[ix]])
			sum_total = np.nansum(Par_adj[i,:])
			Tha_BNWR[i] = sum_between_weight / sum_total
			Tha_BNWR[i] = np.nan_to_num(Tha_BNWR[i])

			Tha_NNCs[i] = len(np.unique(Cortical_plus_thalamus_CI[Par_adj[i,]!=0]))
		BNWRs += Tha_BNWR
		NNCs += Tha_NNCs

	##Cortical nodal roles
	#remove weights connected to low SNR communities (CI==0, orbital frontal, inferior temporal)
	Cortical_adj[np.isnan(Cortical_adj)] = 0
	Cortical_adj[Cortical_CI==0,:]=0
	Cortical_adj[:,Cortical_CI==0]=0

	Cortical_PCs = np.zeros(Cortical_CI.size)
	Cortical_BNWR = np.zeros(Cortical_CI.size)
	Cortical_NNCs = np.zeros(Cortical_plus_thalamus_CI.size)
	for ix, c in enumerate(np.arange(0.01,0.16, 0.01)):
		M = bct.threshold_proportional(Cortical_adj, c, copy=True)
		#PC
		Cortical_PCs += bct.participation_coef(M, Cortical_CI)

		#BNWR and NNC
		BNWR = np.zeros(Cortical_CI.size)
		Cor_NNCs = np.zeros(Cortical_plus_thalamus_CI.size)
		for i in range(len(Cortical_CI)):
			sum_between_weight = np.nansum(M[i,Cortical_CI!=Cortical_CI[i]])
			sum_total = np.nansum(M[i,:])
			BNWR[i] = sum_between_weight / sum_total
			BNWR[i] = np.nan_to_num(BNWR[i])

			Cor_NNCs[i] = len(np.unique(Cortical_CI[M[i,]!=0]))
		Cortical_BNWR += BNWR	
		Cortical_NNCs += Cor_NNCs

	#do WMD, first convert matrices to binary, then calcuate z score using mean and std of "corticocortical degrees"
	Cortical_wm_mean = {}
	Cortical_wm_std = {}
	Cortical_WMDs = np.zeros(Cortical_CI.size)
	WMDs = np.zeros(Cortical_plus_thalamus_CI.size)
	for ix, c in enumerate(np.arange(0.01,0.16, 0.01)):
		
		#threshold by density 
		M = bct.weight_conversion(bct.threshold_proportional(Cortical_adj, c, copy=True), 'binarize')
		Cortical_WMDs += bct.module_degree_zscore(M, Cortical_CI)

		#return mean and degree 
		for CI in np.unique(Cortical_CI):
			Cortical_wm_mean[ix+1, CI] = np.nanmean(np.sum(M[Cortical_CI==CI,:],1))
			Cortical_wm_std[ix+1, CI] = np.nanstd(np.sum(M[Cortical_CI==CI,:],1))

		#thalamic WMD, threshold by density
		M = bct.weight_conversion(bct.threshold_absolute(Thalamocortical_corrmat, cost_thresholds[ix], copy=True), 'binarize')	

		tha_wmd = np.zeros(Cortical_plus_thalamus_CI.size)
		for i in np.unique(Cortical_CI):
			tha_wmd[Cortical_plus_thalamus_CI==i] = (np.sum(M[Cortical_plus_thalamus_CI==i][:, Cortical_ROIs_positions],1)\
			- Cortical_wm_mean[ix+1,i])/Cortical_wm_std[ix+1,i]
		tha_wmd = np.nan_to_num(tha_wmd)
		WMDs += tha_wmd

	# organize output	
	NNCs[Cortical_ROIs_positions] = Cortical_NNCs[Cortical_ROIs_positions]
	BNWRs[Cortical_ROIs_positions] = Cortical_BNWR[Cortical_ROIs_positions]
	PCs[Cortical_ROIs_positions] = Cortical_PCs[Cortical_ROIs_positions]	
	WMDs[Cortical_ROIs_positions] = Cortical_WMDs[Cortical_ROIs_positions]
	# average across thresholds, convert into percentage
	NNCs = (NNCs/15.0) * 100
	BNWRs = (BNWRs/15.0) * 100
	PCs = (PCs/13.5) * 100 #this is the thoretical upperbound
	WMDs = (WMDs/15.0) * 100
 	
	return PCs, BNWRs, NNCs, WMDs


def cal_within_thalamus_nodal_roles(Thalamus_corrmat, Thalamus_CIs, Thalamus_voxel_positions):
	'''Thalamus's nodal properties based on within thalamus weight

	usage: PC, WMD = al_within_thalamus_nodal_roles(Thalamus_corrmat, Thalamus_CIs, Thalamus_voxel_positions)

	'''
	#get PC and WMD
	within_Tha_PCs = np.zeros(Thalamus_voxel_positions.size)
	within_Tha_WMDs = np.zeros(Thalamus_voxel_positions.size)

	for c in np.arange(0.01,0.16,0.01):
		
		graph = matrix_to_igraph(Thalamus_corrmat.copy(),cost=c)
		nodal_graph = brain_graph(VertexClustering(graph, membership=Thalamus_CIs))
		tmp_pc = np.nan_to_num(nodal_graph.pc)
		tmp_wmd = np.nan_to_num(nodal_graph.wmd)
		within_Tha_PCs += tmp_pc
		within_Tha_WMDs += tmp_wmd

	within_Tha_PCs_percentage = (within_Tha_PCs/13.5)*100 
	within_Tha_WMDs_percentage = (within_Tha_WMDs/15)*100

	return within_Tha_PCs_percentage, within_Tha_WMDs_percentage

def sort_CIs(Thalamo_ParcelCIs, Thalamus_voxels):
	Thalamus_CIs = np.zeros(len(Thalamus_voxels))
	for i, thalamus_voxel_index in enumerate(Thalamus_voxels):
		Thalamus_CIs[i] = Thalamo_ParcelCIs[thalamus_voxel_index][0]
	Thalamus_CIs = Thalamus_CIs.astype(int)
	return Thalamus_CIs
