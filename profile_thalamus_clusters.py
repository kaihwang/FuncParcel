### script to profile the distribution of cortical components each thalamus cluster is connected to, "pooled" across subjects
import sys
from brain_graphs import *
from FuncParcel import *
import glob
from scipy.stats.mstats import zscore as zscore

### setup paths and variables needed
corrmat_path = '/home/despoB/connectome-thalamus/Partial_CorrMats/'
path_to_ROIs = '/home/despoB/connectome-thalamus/ROIs'
output_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/Indiv_consensus/'
Parcel_path = '/home/despoB/connectome-thalamus/Thalamic_parcel'

### load data needed
Cortical_ROIs_CI = np.loadtxt(path_to_ROIs + '/Gordon_consensus_CI')
Cortical_ROI_positions = np.array(range(0, len(Cortical_ROIs_CI)))
Thalamocortical_aveMat = np.loadtxt(Parcel_path + '/MGH_Gordon_333_thalamocortical_pcorr_avemat')
Thalamus_voxels = np.loadtxt(path_to_ROIs+'/thalamus_voxel_indices', dtype = int)
Thalamus_voxel_positions = np.array(range(len(Cortical_ROIs_CI),len(Thalamocortical_aveMat)))

#Cortical_Network_CIs = np.loadtxt(path_to_ROIs + '/Gordon_Network_CIs')
#Cortical_Network_CIs_plus_thalamus = np.loadtxt(path_to_ROIs +'/Gordon_CI_plus_thalamus_ROIs')
# Thalamus_voxel_coordinate = np.loadtxt(path_to_ROIs +'/thalamus_voxels_ijk_indices', dtype = int)

#### load the subject's partial corr matrix (bertween the thalamus and the cortex)
# the file pattern, MGH subjects (303 of them), and then with this ROI template: MGH_Gordon_333_consensus_CI*
# MGH*MGH_Gordon_333_consensus_CI*

# file_pattern = 'NKI*Gordon*consensus*'
# pcorrMat_Files = glob.glob(corrmat_path + file_pattern)

# Indiv_Subject_Thalamus_Winner_CIs = np.empty((len(Thalamus_voxels), len(pcorrMat_Files)))

# for i, fn in enumerate(pcorrMat_Files):
# 	M = pickle.load(open(fn, "rb"))

#### step 2:, do the winner take all of thalamus voxel mapping for each subject. 
# 	_, Thalamo_rankCIs, _, = parcel_subcortical_network(M, Cortical_Network_CIs_plus_thalamus, Thalamus_voxels, Cortical_Network_CIs, Cortical_Network_CIs)
# 	# then sourt through the ranking output, and keep the "winner"
# 	Indiv_Subject_Thalamus_Winner_CIs[:,i] = sort_CIs(Thalamo_rankCIs)


#### Get a histogram  for each thalamus cluster
Thalamus_clusters = np.loadtxt(Parcel_path + '/NKI_thalamus_clusters_c07-15')
Cluster_profile = {}
for i in np.unique(Thalamus_clusters):
	tmp = np.zeros(len(np.unique(Cortical_ROIs_CI)))
	for j, n in enumerate(np.unique(Cortical_ROIs_CI)):
		tmp[j] = Thalamocortical_aveMat[Thalamus_voxel_positions[Thalamus_clusters==i],:][:,Cortical_ROI_positions[Cortical_ROIs_CI==n]].mean()
	Cluster_profile[i] = tmp	
save_object(Cluster_profile, Parcel_path + '/NKI_c07-15_Cluster_profiles')


Thalamus_parcel_CI = np.zeros(len(Thalamus_voxels))
for i in np.unique(Thalamus_clusters):
	Thalamus_parcel_CI[Thalamus_clusters==i] = np.argmax(zscore(Cluster_profile[i]))

np.savetxt(Parcel_path + '/NKI_Thalamus_clusters_sort_by_max_corticalCI', Thalamus_parcel_CI)

atlas_path = path_to_ROIs+'/Thalamus_indices.nii.gz' 
ROI_list = path_to_ROIs + '/thalamus_voxel_indices' 
image_path = Parcel_path + '/NKI_gordon_consensus_based_thalamus_parcels_sort_by_max_corticalCIs.nii.gz' 
make_image(atlas_path, image_path, ROI_list, Thalamus_parcel_CI)



#### new profile summary dictionary

#### look at segmentation confidence
os.chdir(corrmat_path)
#files = glob.glob('MGH*MGH_Gordon_333_consensus_CI_pcorr_mat')
#seg_confidence = np.zeros(3539)
#for fn in files:
#	M = pickle.load(open(fn, "rb"))
#	seg_confidence += np.sort(M[10:,0:10])[:,-2:][:,0] / np.sort(M[10:,0:10])[:,-2:][:,1]

M = average_corrmat('MGH*MGH_Gordon_333_consensus_CI_pcorr_mat')
seg_confidence = np.sort(M[10:,0:10])[:,-2:][:,0] / np.sort(M[10:,0:10])[:,-2:][:,1]
#seg_confidence = seg_confidence/ len(files)
Thalamus_parcel_CI = np.loadtxt(Parcel_path + '/Thalamus_clusters_corticalCI')