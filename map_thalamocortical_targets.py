from brain_graphs import *
from FuncParcel import *

#### paths and data
AvgMat_path = '/home/despoB/connectome-thalamus/AvgMatrices'
Parcel_path = '/home/despoB/connectome-thalamus/Thalamic_parcel'
roi_template = 'Craddock_300_plus_thalamus_ROIs_ncsreg' 
path_to_ROIs = '/home/despoB/connectome-thalamus/ROIs'
path_to_data_folder = '/home/despoB/kaihwang/bin/FuncParcel/Data'

ROIs = np.loadtxt(path_to_ROIs+'/Craddock_300_cortical_plus_thalamus_ROIs', dtype = int)
Thalamus_voxels = np.loadtxt(path_to_ROIs+'/thalamus_voxel_indices', dtype = int)
Cortical_ROIs = np.loadtxt(path_to_ROIs+'/Craddock_300_cortical_ROIs', dtype = int)
Cortical_ROIs_positions = np.arange(0,320,1)

MGH_thalamocor_adj = np.loadtxt(Parcel_path+'/MGH_Craddock_300_cortical_plus_thalamus_parcorrmatavg')
NKI_thalamocor_adj = np.loadtxt(Parcel_path+'/NKI_mx_645_Craddock_300_cortical_plus_thalamus_parcorrmatavg')

#### call function and map targets for MGH and NKI data, compare them
MGH_Cortical_targets, MGH_Cortical_nontargets, MGH_cost_thresholds = map_subcortical_cortical_targets(MGH_thalamocor_adj, Cortical_ROIs, Thalamus_voxels)
NKI_Cortical_targets, NKI_Cortical_nontargets, NKI_cost_thresholds = map_subcortical_cortical_targets(NKI_thalamocor_adj, Cortical_ROIs, Thalamus_voxels)

# lets intersect the two dataset to find reliable targets... and non-targets....
Cortical_targets={}
Cortical_nontargets ={}
for i in Thalamus_voxels:
	Cortical_targets[i] = np.intersect1d(NKI_Cortical_targets[15,i], MGH_Cortical_targets[15,i])
	Cortical_nontargets[i] = np.intersect1d(NKI_Cortical_nontargets[15,i], MGH_Cortical_nontargets[15,i])

save_object(Cortical_targets, path_to_data_folder +'/Cortical_targets')
save_object(Cortical_nontargets, path_to_data_folder +'/Cortical_nontargets')

#### do a count for each ROI, how many thalamic voxels it connects, save result to nifti
targets = np.zeros(0, dtype='int')
for v in Cortical_targets.itervalues():
	targets=np.concatenate((targets,v))
np.savetxt(path_to_data_folder+'/targets_across_costs.txt', targets)

# save to nifti
Cortical_ROI_target_count = np.zeros(len(Cortical_ROIs), dtype='int')
for i, roi in enumerate(Cortical_ROIs):
	Cortical_ROI_target_count[i] = np.sum(targets == roi)

atlas_path = path_to_ROIs+'/Craddock_300_cortical.nii.gz' 
image_path = Parcel_path +'/MGH+NKI_cortical_roi_target_count.nii.gz' 
make_image(atlas_path, image_path, Cortical_ROIs, Cortical_ROI_target_count)