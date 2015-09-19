#scrip to analyze thalamus nodal proerty

from brain_graphs import *
from FuncParcel import *
import matplotlib.pylab as plt

#### setup
AvgMat_path = '/home/despoB/connectome-thalamus/AvgMatrices/'
roi_template = 'Craddock_300_plus_thalamus_ROIs_ncsreg' 
path_to_ROIs = '/home/despoB/connectome-thalamus/ROIs'
path_to_data_folder = '/home/despoB/kaihwang/bin/FuncParcel/Data'

ROIs = np.loadtxt(path_to_ROIs+'/Craddock_300_cortical_plus_thalamus_ROIs', dtype = int)
Cortical_ROIs_positions = np.arange(0,320,1)
Thalamus_voxel_positions = np.arange(320,3859,1)	
Thalamus_voxel_coordinate = np.loadtxt(path_to_ROIs +'/thalamus_voxels_ijk_indices', dtype = int)


## do network based parcellation
path_to_adjmat = '/home/despoB/connectome-thalamus/Thalamic_parcel/NKI_mx_645_cortical_network_plus_thalamus_corrmatavg'
path_to_list_of_subcorticalcortical_ROIs = '/home/despoB/connectome-thalamus/ROIs/Cortical_CI_plus_thalamus_ROIs'
path_to_list_of_subcortical_voxels = '/home/despoB/connectome-thalamus/ROIs/thalamus_voxel_indices'
path_to_list_of_cortical_ROIs ='/home/despoB/connectome-thalamus/ROIs/Cortical_Network_CIs'
path_to_Cortical_CI = '/home/despoB/connectome-thalamus/ROIs/Cortical_Network_CIs'

#call parcel function
_, Thalamo_ParcelCIs, _, = parcel_subcortical_network(path_to_adjmat, \
            path_to_list_of_subcorticalcortical_ROIs, \
            path_to_list_of_subcortical_voxels, path_to_list_of_cortical_ROIs, path_to_Cortical_CI)

# sort CI vector
Thalamus_CIs = np.zeros(len(Thalamus_voxel_coordinate))
for i, thalamus_voxel_index in enumerate(Thalamus_voxel_coordinate[:,3]):
	Thalamus_CIs[i]= Thalamo_ParcelCIs[thalamus_voxel_index][0]

# show histogram of CI distribution
Thalamus_CIs = Thalamus_CIs.astype(int)
plt.hist(Thalamus_CIs, bins=np.unique(Thalamus_CIs))
#plt.show()

#write out parcel info into nifti
atlas_path = path_to_ROIs+'/Thalamus_indices.nii.gz' 
ROI_list = path_to_ROIs + '/thalamus_voxel_indices' 
image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/MGH_network_based_thalamus_parcel.nii.gz' 
make_image(atlas_path, image_path, ROI_list, Thalamus_CIs)

##Here analyze each thalamic voxel's nodal property
#first load thalamus voxel + cortical ROI's matrix
path = AvgMat_path + 'MGH/_%s_corrmatavg' %(roi_template)
adj = np.loadtxt(path)
adj[adj<0] = 0.0
np.fill_diagonal(adj,0)

#clip within thalamus edges
adj[320:3859,320:3859] = 0.0

# load cortical CI vectors
partition = pickle.load(open(path_to_data_folder+"/Graph/graph_MGH_Craddock_300_cortical_min_cost0.02_min_weight2.0_min_size5", "rb"))
Cortical_CI=np.array(partition.community.membership) +1
Cortical_CI[Cortical_CI>10] = 0

#clip edges of nodes with no good CI assignment (CI==0)
Cortical_plus_thalamus_CI = np.concatenate((Cortical_CI, Thalamus_CIs), 1) 
adj[Cortical_plus_thalamus_CI==0,:]=0
adj[:,Cortical_plus_thalamus_CI==0]=0

# calculate cortical node metrics
#cortical_adj = adj[0:320,0:320].copy()
#graph = matrix_to_igraph(cortical_adj,cost=0.15)
#cortical_graph = brain_graph(VertexClustering(graph, membership=Cortical_CI))

# calculate thalamic voxel nodal metrics
adj[0:320,0:320] = 0.0 #only thalamo-cortical edges are preserved
graph = matrix_to_igraph(adj.copy(),cost=0.0087)
nodal_graph = brain_graph(VertexClustering(graph, membership=Cortical_plus_thalamus_CI))
save_path = '/home/despoB/kaihwang/bin/FuncParcel/Data/Graph/graph_MGH_%s' %roi_template
save_object(nodal_graph, save_path)

#get cortical PC and wmd
cortical_pc = partition.pc.copy()
cortical_pc[np.isnan(cortical_pc)]=0.0
cortical_wmd = partition.wmd.copy()
cortical_wmd[np.isnan(cortical_wmd)]=0.0

atlas_path = path_to_ROIs+'/Craddock_300_plus_thalamus_ROIs.nii.gz' 
ROI_list = path_to_ROIs+'/Craddock_300_cortical_plus_thalamus_ROIs' 

pc = nodal_graph.pc
wmd = nodal_graph.wmd
pc[0:320] = cortical_pc
wmd[0:320] = cortical_wmd
Thalamus_pc = pc[320::]

pc = pc *100
wmd = wmd * 100
image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/MGH_%s_nodal_role_pc.nii.gz' %roi_template
make_image(atlas_path, image_path, ROI_list, pc.astype(int))

image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/MGH_%s_nodal_role_wmd.nii.gz' %roi_template
make_image(atlas_path, image_path, ROI_list, wmd.astype(int))

#get descriptive stats of PC within parcel
mean_PC_per_thalamus_parcel = np.zeros(max(Thalamus_CIs))
for ix, i in enumerate(np.arange(1,max(Thalamus_CIs)+1,1)):
	mean_PC_per_thalamus_parcel[ix] = np.nanmean(Thalamus_pc[Thalamus_CIs==i])

#get ratio of between/total connectivity weight
adjt = bct.threshold_proportional(adj,0.015)

thalamus_between_ratio = np.zeros(len(Thalamus_voxel_positions))
for ix, i in enumerate(Thalamus_voxel_positions):
	sum_between_weight = np.nansum(adjt[i,Cortical_plus_thalamus_CI!=Thalamus_CIs[ix]])
	sum_total = np.nansum(adjt[i,:])
	thalamus_between_ratio[ix] = sum_between_weight / sum_total

mean_ratio_per_thalamus_parcel = np.zeros(max(Thalamus_CIs))	
for ix, i in enumerate(np.arange(1,max(Thalamus_CIs)+1,1)):
	mean_ratio_per_thalamus_parcel[ix] = np.nanmean(thalamus_between_ratio[Thalamus_CIs==i])

#write to nifti
thalamus_between_ratio = thalamus_between_ratio*100
atlas_path = path_to_ROIs+'/Thalamus_indices.nii.gz' 
image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/MGH_%s_nodal_role_bcr.nii.gz' %roi_template
make_image(atlas_path, image_path, ROI_list, pc.astype(int))




