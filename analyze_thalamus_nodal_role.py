#scrip to analyze thalamus nodal proerties
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

def subject_nodal_role(subject):	
	print('\n'+subject)
	
	#### run nodal role function
	Thalamocor_adj = pickle.load(open(Partialcorr_path + '/MGH_' + subject + '_Ses1_Gordon_333_cortical_pcorr_mat', "rb"))
	Cortical_adj = np.loadtxt(Adjmat_path + '/MGH_' + subject + '_Ses1_Gordon_333_cortical_corrmat')
	print('\ndata succesfully loaded')

	_, _, cost_thresholds = map_subcortical_cortical_targets(Thalamocor_adj, Cortical_ROIs, Thalamus_voxels)
	
	PCs, BNWRs, NNCs, WMDs = cal_thalamus_and_cortical_ROIs_nodal_properties(Thalamocor_adj, \
		Cortical_adj, \
		Cortical_plus_thalamus_CI, \
		Thalamus_CIs, \
		Cortical_CI, \
		Cortical_ROIs_positions, \
		Thalamus_voxel_positions, \
		cost_thresholds)

	#### save output
	file_path = path_to_data_folder + '/MGH_%s_tha_nodal_pcorr_PCs' %subject
	save_object(PCs, file_path)

	file_path = path_to_data_folder + '/MGH_%s_tha_nodal_pcorr_BNWRs' %subject
	save_object(BNWRs, file_path)

	file_path = path_to_data_folder + '/MGH_%s_tha_nodal_pcorr_NNCs' %subject
	save_object(NNCs, file_path)

	file_path = path_to_data_folder + '/MGH_%s_tha_nodal_pcorr_WMDs' %subject
	save_object(WMDs, file_path)


def ave_subject_nodal_role(dset, pattern):
	global path_to_data_folder
	'''average individual subject's nodal role output'''
	
	files = glob.glob(path_to_data_folder+'/'+dset+pattern)
	
	M_Sum =[]
	for fn in files:
		M = pickle.load(open(fn, "rb"))
		M_Sum += [M]
	aveM = np.array(sum(M_Sum)/len(files))


def group_parcel_role(Datasets):
	'''calculate nodal metrics on group averaged matrices, unit of analysis at the thalamus is "parcel"'''

	for dat in Datasets:
		
		Thalamocor_adj = np.loadtxt(Parcel_path + '/' + dat + 'Gordon_ThaWTA_avemat')
		Cortical_adj = np.loadtxt(Parcel_path + '/' + dat + 'Gordon_333_cortical_avemat')

		_, _, cost_thresholds = map_subcortical_cortical_targets(Thalamocor_adj, Cortical_ROIs, Thalamus_parcel_positions)

		PCs, BNWRs, NNCs, WMDs, bPCs, mean_NNC, mean_BNWR, mean_PC, mean_bPC, mean_WMD = cal_thalamus_and_cortical_ROIs_nodal_properties(Thalamocor_adj, \
				Cortical_adj, \
				Cortical_plus_thalamus_parcel_CI, \
				Thalamus_Parcels, \
				Cortical_CI, \
				Cortical_ROIs_positions, \
				Thalamus_parcel_positions, \
				cost_thresholds)

		file_path = path_to_data_folder + '/%savemat_tha_parcel_pcorr_PCs' %dat
		save_object(PCs, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_parcel_pcorr_BNWRs' %dat
		save_object(BNWRs, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_parcel_pcorr_NNCs' %dat
		save_object(NNCs, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_parcel_pcorr_WMDs' %dat
		save_object(WMDs, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_parcel_pcorr_bPCs' %dat
		save_object(bPCs, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_parcel_pcorr_meanPC' %dat
		save_object(mean_PC, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_parcel_pcorr_meanBNWR' %dat
		save_object(mean_BNWR, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_parcel_pcorr_meanNNC' %dat
		save_object(mean_NNC, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_parcel_pcorr_meanWMD' %dat
		save_object(mean_WMD, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_parcel_pcorr_meanbPC' %dat
		save_object(mean_bPC, file_path)


def group_voxel_role(Datasets):
	'''calculate nodal metrics on group averaged matrices for "voxels" at the thalamus level'''

	for dat in Datasets:
		Thalamocor_adj = np.loadtxt(Parcel_path + '/' + dat + 'Gordon_333_thalamocortical_pcorr_avemat')
		Cortical_adj = np.loadtxt(Parcel_path + '/' + dat + 'Gordon_333_cortical_avemat')

		_, _, cost_thresholds = map_subcortical_cortical_targets(Thalamocor_adj, Cortical_ROIs, Thalamus_voxels)

		PCs, BNWRs, NNCs, WMDs, bPCs, mean_NNC, mean_BNWR, mean_PC, mean_bPC, mean_WMD = cal_thalamus_and_cortical_ROIs_nodal_properties(Thalamocor_adj, \
				Cortical_adj, \
				Cortical_plus_thalamus_CI, \
				Thalamus_CIs, \
				Cortical_CI, \
				Cortical_ROIs_positions, \
				Thalamus_voxel_positions, \
				cost_thresholds)

		file_path = path_to_data_folder + '/%savemat_tha_nodal_pcorr_PCs' %dat
		save_object(PCs, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_nodal_pcorr_BNWRs' %dat
		save_object(BNWRs, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_nodal_pcorr_NNCs' %dat
		save_object(NNCs, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_nodal_pcorr_WMDs' %dat
		save_object(WMDs, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_nodal_pcorr_bPCs' %dat
		save_object(bPCs, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_nodal_pcorr_meanPC' %dat
		save_object(mean_PC, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_nodal_pcorr_meanBNWR' %dat
		save_object(mean_BNWR, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_nodal_pcorr_meanNNC' %dat
		save_object(mean_NNC, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_nodal_pcorr_meanWMD' %dat
		save_object(mean_WMD, file_path)

		file_path = path_to_data_folder + '/%savemat_tha_nodal_pcorr_meanbPC' %dat
		save_object(mean_bPCs, file_path)


if __name__ == "__main__":
#subject = sys.stdin.read().strip('\n')	
#subject_nodal_role(subject)

### group averaged
	#Datasets =	['MGH_', 'NKI_645_', 'NKI_1400_'] 
	Datasets =	['MGH_'] 
	group_voxel_role(Datasets)
	#group_parcel_role(Datasets)

