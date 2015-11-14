### script to profile the distribution of cortical components each thalamus cluster is connected to, "pooled" across subjects
import sys
from brain_graphs import *
from FuncParcel import *
import glob

### setup paths and variables needed
corrmat_path = '/home/despoB/connectome-thalamus/Partial_CorrMats/'
path_to_ROIs = '/home/despoB/connectome-thalamus/ROIs'
output_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/Indiv_consensus/'
Parcel_path = '/home/despoB/connectome-thalamus/Thalamic_parcel'
Cortical_Network_CIs = np.loadtxt(path_to_ROIs + '/Gordon_Network_CIs')
Cortical_Network_CIs_plus_thalamus = np.loadtxt(path_to_ROIs +'/Gordon_CI_plus_thalamus_ROIs')
Thalamus_voxels = np.loadtxt(path_to_ROIs+'/thalamus_voxel_indices', dtype = int)
Thalamus_voxel_positions = np.arange(320,3859,1)
Thalamus_voxel_coordinate = np.loadtxt(path_to_ROIs +'/thalamus_voxels_ijk_indices', dtype = int)

def sort_CIs(Thalamo_ParcelCIs):
	global Thalamus_voxel_coordinate
	Thalamus_CIs = np.zeros(len(Thalamus_voxel_coordinate))
	for i, thalamus_voxel_index in enumerate(Thalamus_voxel_coordinate[:,3]):
		Thalamus_CIs[i] = Thalamo_ParcelCIs[thalamus_voxel_index][0]
	Thalamus_CIs = Thalamus_CIs.astype(int)
	return Thalamus_CIs

#step 1: load the subject's partial corr matrix (bertween the thalamus and the cortex)
# the file pattern, MGH subjects (303 of them), and then with this ROI template: MGH_Gordon_333_consensus_CI*
# MGH*MGH_Gordon_333_consensus_CI*

file_pattern = 'NKI*Gordon*consensus*'
pcorrMat_Files = glob.glob(corrmat_path + file_pattern)

Indiv_Subject_Thalamus_Winner_CIs = np.empty((len(Thalamus_voxels), len(pcorrMat_Files)))

for i, fn in enumerate(pcorrMat_Files):
	M = pickle.load(open(fn, "rb"))

	# step 2:, do the winner take all of thalamus voxel mapping for each subject. 
	_, Thalamo_rankCIs, _, = parcel_subcortical_network(M, Cortical_Network_CIs_plus_thalamus, Thalamus_voxels, Cortical_Network_CIs, Cortical_Network_CIs)
	# then sourt through the ranking output, and keep the "winner"
	Indiv_Subject_Thalamus_Winner_CIs[:,i] = sort_CIs(Thalamo_rankCIs)


# step3, get a histogram counter for each thalamus cluster
#load the group based clusters
Thalamus_clusters = np.loadtxt(Parcel_path + '/Thalamus_clusters')

# get the distribution of each voxel's connected cortical component
for i in np.unique(Thalamus_clusters):
	print(Counter(Indiv_Subject_Thalamus_Winner_CIs[Thalamus_clusters ==i].flatten()))
	print(Counter(stats.mode(Indiv_Subject_Thalamus_Winner_CIs[Thalamus_clusters ==i])[0][0]))