#from brain_graphs import *
from FuncParcel import *
from scipy.stats.mstats import zscore as zscore

#### paths and data
AvgMat_path = '/home/despoB/connectome-thalamus/AvgMatrices'
Parcel_path = '/home/despoB/connectome-thalamus/Thalamic_parcel'
path_to_ROIs = '/home/despoB/connectome-thalamus/ROIs'
path_to_data_folder = '/home/despoB/kaihwang/bin/FuncParcel/Data'

ROIs = np.loadtxt(path_to_ROIs+'/Gordon_333_plus_thalamus', dtype = int)
Thalamus_voxels = np.loadtxt(path_to_ROIs+'/thalamus_voxel_indices', dtype = int)
Cortical_ROIs = np.loadtxt(path_to_ROIs+'/Gordon_333', dtype = int)
Cortical_ROIs_positions = np.arange(0,333,1)

MGH_thalamocor_adj = np.loadtxt(Parcel_path+'/MGH_Gordon_333_thalamocortical_pcorr_avemat')
#NKI_thalamocor_adj = np.loadtxt(Parcel_path+'/NKI_mx_645_Craddock_300_cortical_plus_thalamus_parcorrmatavg')

#### call function and map targets for MGH and NKI data, compare them
MGH_Cortical_targets, MGH_Cortical_nontargets, MGH_cost_thresholds = map_subcortical_cortical_targets(MGH_thalamocor_adj, Cortical_ROIs, Thalamus_voxels)
#NKI_Cortical_targets, NKI_Cortical_nontargets, NKI_cost_thresholds = map_subcortical_cortical_targets(NKI_thalamocor_adj, Cortical_ROIs, Thalamus_voxels)


#### get AN, MD, PuM, intralaminar cortical targets
AN_vox = np.loadtxt(path_to_ROIs + '/Morel_AN_indices')
MD_vox = np.loadtxt(path_to_ROIs + '/Morel_MD_indices')
PuM_vox = np.loadtxt(path_to_ROIs + '/Morel_PuM_indices')
Intralaminar_vox = np.loadtxt(path_to_ROIs + '/Morel_intralaminar_indices')
LP_vox = np.loadtxt(path_to_ROIs + '/Morel_LP_indices')
VL_vox = np.loadtxt(path_to_ROIs + '/Morel_VL_indices')
LGN_vox = np.loadtxt(path_to_ROIs + '/Morel_LGN_indices')
VA_vox = np.loadtxt(path_to_ROIs + '/Morel_VA_indices')
VM_vox = np.loadtxt(path_to_ROIs + '/Morel_VM_indices')

def map_vox_target(voxlist, targetlist):
	targets = np.zeros(0, dtype='int')
	for i, vox in enumerate(voxlist):
		if vox !=0:
			targets = np.concatenate((MGH_Cortical_targets[2,vox], targets))
	targets = np.unique(targets)
	return targets

# targets = map_vox_target(AN_vox, MGH_Cortical_targets)
# np.savetxt(path_to_ROIs + '/Morel_AN_targets', targets)
# targets = map_vox_target(MD_vox, MGH_Cortical_targets)
# np.savetxt(path_to_ROIs + '/Morel_MD_targets', targets)
# targets = map_vox_target(PuM_vox, MGH_Cortical_targets)
# np.savetxt(path_to_ROIs + '/Morel_PuM_targets', targets)
# targets = map_vox_target(Intralaminar_vox, MGH_Cortical_targets)
# np.savetxt(path_to_ROIs + '/Morel_intralaminar_targets', targets)
# targets = map_vox_target(CL_vox, MGH_Cortical_targets)
# np.savetxt(path_to_ROIs + '/Morel_CL_targets', targets)
# targets = map_vox_target(VL_vox, MGH_Cortical_targets)
# np.savetxt(path_to_ROIs + '/Morel_VL_targets', targets)


### cal average con strength between nuclei and cortical networks

MGHadjmat = np.loadtxt('/home/despoB/connectome-thalamus/Thalamic_parcel/MGH_Gordon_ConsesnsusCI_thalamocortical_pcorr_avemat')
# load subcortical voxel info
Thalamus_voxel_coordinate = np.loadtxt('/home/despoB/connectome-thalamus/ROIs/thalamus_voxels_ijk_indices', dtype = int)

subcorticalcortical_ROIs = np.loadtxt('/home/despoB/connectome-thalamus/ROIs/Gordon_CI_plus_thalamus_ROIs')
subcortical_voxels = np.loadtxt('/home/despoB/connectome-thalamus/ROIs/thalamus_voxel_indices')
cortical_ROIs =np.loadtxt('/home/despoB/connectome-thalamus/ROIs/Gordon_Network_consensusCI')
Cortical_CI = np.loadtxt('/home/despoB/connectome-thalamus/ROIs/Gordon_Network_consensusCI')

#call function
_, _, tha_cor_str, = parcel_subcortical_network(MGHadjmat, subcorticalcortical_ROIs, \
            subcortical_voxels, cortical_ROIs, Cortical_CI)

#loop through voxels
def get_vox_str(voxlist, targetlist):
	targets = np.zeros(9, dtype='int')
	for i, vox in enumerate(voxlist):
		if vox !=0:
			targets = targets + zscore(targetlist[vox])
	targets = targets/len(voxlist)		
	return targets


