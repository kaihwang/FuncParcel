#scrip to analyze cortical ROI nodal proerties, the other script taking to long (too many thalamic voxel)
#from brain_graphs import *
from FuncParcel import *

#### setup
Parcel_path = '/home/despoB/connectome-thalamus/Thalamic_parcel'
Partialcorr_path = '/home/despoB/connectome-thalamus/Partial_CorrMats'
Adjmat_path = '/home/despoB/connectome-thalamus/NotBackedUp/AdjMatrices'
path_to_ROIs = '/home/despoB/connectome-thalamus/ROIs'
path_to_data_folder = '/home/despoB/kaihwang/Rest/Graph'

Thalamus_CIs = np.loadtxt(Parcel_path + '/MGH_Thalamus_WTA_CIs')
Thalamus_Parcels = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
Cortical_CI = np.loadtxt(path_to_ROIs + '/Gordon_consensus_CI')
Cortical_plus_thalamus_CI = np.append(Cortical_CI, Thalamus_CIs)
Cortical_plus_thalamus_parcel_CI = np.append(Cortical_CI, Thalamus_Parcels)
Cortical_ROIs_positions = np.arange(0,len(Cortical_CI),1)
Thalamus_voxel_positions = np.arange(len(Cortical_CI),len(Cortical_plus_thalamus_CI),1)
Thalamus_parcel_positions = np.arange(len(Cortical_CI),len(Cortical_plus_thalamus_parcel_CI),1)
Thalamus_voxels = np.loadtxt(path_to_ROIs+'/thalamus_voxel_indices', dtype = int)
Cortical_ROIs = np.loadtxt(path_to_ROIs+'/Gordon_333', dtype = int)
ROIs = np.append(Cortical_ROIs, Thalamus_voxels)

Cortical_adj  = average_corrmat('/home/despoB/kaihwang/Rest/NotBackedUp/ParMatrices/MGH*_Gordon_333_cortical_corrmat', np_txt=True, pickle_object=False)

Cortical_adj[np.isnan(Cortical_adj)] = 0  
Cortical_PCs = []#np.zeros(Cortical_CI.size)
Cortical_bPCs = []#np.zeros(Cortical_CI.size)
Cortical_BNWRs = []#np.zeros(Cortical_CI.size)
Cortical_NNCs = []#np.zeros(Cortical_plus_thalamus_CI.size)

for ix, c in enumerate(np.arange(0.01,0.16, 0.01)):
	M = bct.threshold_proportional(Cortical_adj, c, copy=True)
	bM = bct.weight_conversion(M, 'binarize', copy=True);

	#PC
	Cortical_PCs += [bct.participation_coef(M, Cortical_CI)]
	Cortical_bPCs += [bct.participation_coef(bM, Cortical_CI)]

	#BNWR and NNC
	BNWR = np.zeros(Cortical_CI.size)
	Cor_NNCs = np.zeros(Cortical_plus_thalamus_CI.size)
	for i in range(len(Cortical_CI)):
		sum_between_weight = np.nansum(M[i,Cortical_CI!=Cortical_CI[i]])
		sum_total = np.nansum(M[i,:])
		BNWR[i] = sum_between_weight / sum_total
		BNWR[i] = np.nan_to_num(BNWR[i])

		Cor_NNCs[i] = len(np.unique(Cortical_CI[M[i,]!=0]))
	Cortical_BNWRs += [BNWR]	
	Cortical_NNCs += [Cor_NNCs]

Cortical_WMDs = []#np.zeros(Cortical_CI.size)
WMDs = []#np.zeros(Cortical_plus_thalamus_CI.size)
for ix, c in enumerate(np.arange(0.01,0.16, 0.01)):
		
	#threshold by density 
	bM = bct.weight_conversion(bct.threshold_proportional(Cortical_adj, c, copy=True), 'binarize')
	Cortical_WMDs += [bct.module_degree_zscore(bM, Cortical_CI)]

mean_Cortical_PC = (np.sum(Cortical_PCs,axis=0)/13.5) * 100
mean_Cortical_bPC = (np.sum(Cortical_bPCs,axis=0)/13.5) * 100
mean_Cortical_WMD = (np.sum(Cortical_WMDs,axis=0)/15.0) * 100
mean_Cortical_BNWR = (np.sum(Cortical_BNWRs, axis=0)/15.0) * 100
mean_Cortical_NNC = (np.sum(Cortical_NNCs,axis=0)/15.0) * 100


path_to_data_folder = '/home/despoB/kaihwang/Rest/Graph'
file_path = path_to_data_folder + '/MGH_avemat_cortical_nodal_corr_PCs' 
save_object(Cortical_PCs, file_path)
file_path = path_to_data_folder + '/MGH_avemat_cortical_nodal_corr_WMDs' 
save_object(Cortical_WMDs, file_path)
file_path = path_to_data_folder + '/MGH_avemat_cortical_nodal_corr_meanPCs' 
save_object(mean_Cortical_PC, file_path)
file_path = path_to_data_folder + '/MGH_avemat_cortical_nodal_corr_meanWMDs' 
save_object(mean_Cortical_WMD, file_path)











