#scrip to analyze thalamus nodal proerty

from brain_graphs import *
from FuncParcel import *
import matplotlib.pylab as plt

###### setup
AvgMat_path = '/home/despoB/connectome-thalamus/AvgMatrices'
Parcel_path = '/home/despoB/connectome-thalamus/Thalamic_parcel'
roi_template = 'Craddock_300_plus_thalamus_ROIs_ncsreg' 
path_to_ROIs = '/home/despoB/connectome-thalamus/ROIs'
path_to_data_folder = '/home/despoB/kaihwang/bin/FuncParcel/Data'

ROIs = np.loadtxt(path_to_ROIs+'/Craddock_300_cortical_plus_thalamus_ROIs', dtype = int)
Thalamus_voxels = np.loadtxt(path_to_ROIs+'/thalamus_voxel_indices', dtype = int)
Cortical_ROIs = np.loadtxt(path_to_ROIs+'/Craddock_300_cortical_ROIs', dtype = int)
Cortical_ROIs_positions = np.arange(0,320,1)
Thalamus_voxel_positions = np.arange(320,3859,1)	
Thalamus_voxel_coordinate = np.loadtxt(path_to_ROIs +'/thalamus_voxels_ijk_indices', dtype = int)
Thalamocortical_corrmat = np.loadtxt(Parcel_path+'/MGH_Craddock_300_cortical_plus_thalamus_parcorrmatavg')

Cortical_CI = np.loadtxt(path_to_ROIs+'/Cortical_CI', dtype='int')

###### find cortical target and nontarget, save output 
Cortical_targets, Cortical_nontargets, cost_thresholds = map_subcortical_cortical_targets(Thalamocortical_corrmat, Cortical_ROIs, Thalamus_voxels)
save_object(Cortical_targets, path_to_data_folder +'/Cortical_targets')
save_object(Cortical_nontargets, path_to_data_folder +'/Cortical_nontargets')

targets = np.empty(0, dtype='int')
for v in Cortical_targets.itervalues():
	targets=np.concatenate((targets,v))

np.savetxt(path_to_data_folder+'/targets_across_costs.txt', targets)

#save to nifti
Cortical_ROI_target_count = np.zeros(len(Cortical_ROIs), dtype='int')
for i, roi in enumerate(Cortical_ROIs):
	Cortical_ROI_target_count[i] = np.sum(targets == roi)

atlas_path = path_to_ROIs+'/Craddock_300_cortical.nii.gz' 
image_path = Parcel_path +'/MGH_cortical_roi_target_count.nii.gz' 
make_image(atlas_path, image_path, Cortical_ROIs, Cortical_ROI_target_count)

###### do network based parcellation


#call function
Netadjmat = np.loadtxt('/home/despoB/connectome-thalamus/Thalamic_parcel/MGH_cortical_network_plus_thalamus_parcorrmatavg')
Cortical_Network_CIs = np.loadtxt('/home/despoB/connectome-thalamus/ROIs/Cortical_Network_CIs')
Cortical_Network_CIs_plus_thalamus = np.loadtxt('/home/despoB/connectome-thalamus/ROIs/Cortical_CI_plus_thalamus_ROIs')
#path_to_Cortical_CI = '/home/despoB/connectome-thalamus/ROIs/Cortical_Network_CIs'


#clip out bad SNR networks (orbital frontal and inf temporal)
Netadjmat[9,:]=0
Netadjmat[:,9]=0
Netadjmat[10,:]=0
Netadjmat[:,10]=0

#call function
_, Thalamo_ParcelCIs, _ = parcel_subcortical_network(Netadjmat, Cortical_Network_CIs_plus_thalamus, Thalamus_voxels, Cortical_Network_CIs, Cortical_Network_CIs)
save_object(Thalamo_ParcelCIs, path_to_data_folder +'/Thalamo_ParcelCIs')


#sort CI vector
Thalamus_CIs = np.zeros(len(Thalamus_voxel_coordinate))
for i, thalamus_voxel_index in enumerate(Thalamus_voxel_coordinate[:,3]):
	Thalamus_CIs[i]= Thalamo_ParcelCIs[thalamus_voxel_index][0]

# show histogram of CI distribution
Thalamus_CIs = Thalamus_CIs.astype(int)
plt.hist(Thalamus_CIs, bins=9)
#plt.show()

#write out parcel info into nifti
atlas_path = path_to_ROIs+'/Thalamus_indices.nii.gz' 
ROI_list = path_to_ROIs + '/thalamus_voxel_indices' 
image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/MGH_network_based_thalamus_parcel.nii.gz' 
make_image(atlas_path, image_path, ROI_list, Thalamus_CIs)

###### Here analyze each thalamic voxel's nodal property
#first load thalamus voxel + cortical ROI's matrix, there are two. One partial, one full correlation
#path = AvgMat_path + 'MGH/_Craddock_300_plus_thalamus_ROIs_ncsreg_corrmatavg' 
#Full_adj = np.loadtxt(path)
#Full_adj[Full_adj<0] = 0.0
#np.fill_diagonal(Full_adj,0)


#clip within thalamus edges
#adj[320:3859,320:3859] = 0.0

# load cortical CI vectors
#partition = pickle.load(open(path_to_data_folder+"/Graph/graph_MGH_Craddock_300_cortical_min_cost0.02_min_weight2.0_min_size5", "rb"))
#Cortical_CI=np.array(partition.community.membership) +1
#Cortical_CI[Cortical_CI>10] = 0

#Thalamus_CIs = np.zeros(len(Thalamus_voxel_coordinate))
Cortical_plus_thalamus_CI = np.concatenate((Cortical_CI, Thalamus_CIs), 1) 
Cortical_plus_thalamus_CI = Cortical_plus_thalamus_CI.astype(int)

## calculate PC, BNWR, WMD, iterate through costs 
#PC
Tha_PCs = np.empty(Cortical_plus_thalamus_CI.size)
for c in cost_thresholds:
	Par_adj = Thalamocortical_corrmat.copy()
	Par_adj[Cortical_ROIs_positions[Cortical_CI==0],:]=0
	Par_adj[:,Cortical_ROIs_positions[Cortical_CI==0]]=0
	Par_adj[Par_adj<c]=0
	Tha_PCs += bct.participation_coef(Par_adj, Cortical_plus_thalamus_CI)
	#graph = matrix_to_igraph(Par_adj.copy(),cost=1)
	#nodal_graph = brain_graph(VertexClustering(graph, membership=Cortical_plus_thalamus_CI))
	#save_path = '/home/despoB/kaihwang/bin/FuncParcel/Data/Graph/graph_MGH_%s' %roi_template
	#save_object(nodal_graph, save_path)

Tha_PCs[Cortical_ROIs_positions.astype('int')]=0
Tha_PCs_percentage = (Tha_PCs/13.5)*100 #13.5 is the highest PC possible  (upper bound is .9 for each cost, .9*15 = 13.5)
atlas_path = path_to_ROIs+'/Craddock_300_plus_thalamus_ROIs.nii.gz' 
image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/MGH_tha_nodal_role_pc.nii.gz'
make_image(atlas_path, image_path, ROIs, Tha_PCs_percentage)
save_object(Tha_PCs, path_to_data_folder +'/Tha_PCs')

#BNWR
Tha_BNWR = np.empty(Cortical_plus_thalamus_CI.size)
for c in cost_thresholds:
	Par_adj = Thalamocortical_corrmat.copy()
	#Par_adj[Cortical_ROIs_positions[Cortical_CI==0],:]=0
	#Par_adj[:,Cortical_ROIs_positions[Cortical_CI==0]]=0
	Par_adj[Par_adj<c]=0
	
	BNWR = np.zeros(Cortical_plus_thalamus_CI.size)
	for ix, i in enumerate(Thalamus_voxel_positions):
		sum_between_weight = np.nansum(Par_adj[i,Cortical_plus_thalamus_CI!=Thalamus_CIs[ix]])
		sum_total = np.nansum(Par_adj[i,:])
		BNWR[i] = sum_between_weight / sum_total
		BNWR[i] = np.nan_to_num(BNWR[i])
		#sBNWR = np.concatenate(BNWR,sBNWR)
	Tha_BNWR += BNWR

Tha_BNWR[Cortical_ROIs_positions.astype('int')]=0
Tha_BNWR_percentage = (Tha_BNWR/15)*100 
atlas_path = path_to_ROIs+'/Craddock_300_plus_thalamus_ROIs.nii.gz' 
image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/MGH_tha_nodal_role_bnwr.nii.gz'
make_image(atlas_path, image_path, ROIs, Tha_BNWR_percentage)
save_object(Tha_BNWR, path_to_data_folder +'/Tha_BNWR')


#get cortical PC and BNEW and WMD



#get thlamus WMD, using mean and STD from cortex. 



###### calculate cortical node metrics
#cortical_adj = adj[0:320,0:320].copy()
#graph = matrix_to_igraph(cortical_adj,cost=0.15)
#cortical_graph = brain_graph(VertexClustering(graph, membership=Cortical_CI))

#### get descriptive stats of PC within parcel
mean_PC_per_thalamus_parcel = np.zeros(max(Thalamus_CIs))
for ix, i in enumerate(np.arange(1,max(Thalamus_CIs)+1,1)):
	mean_PC_per_thalamus_parcel[ix] = np.nanmean(Thalamus_pc[Thalamus_CIs==i])

#get ratio of between/total connectivity weight
adjt = bct.threshold_proportional(adj,0.015)



mean_ratio_per_thalamus_parcel = np.zeros(max(Thalamus_CIs))	
for ix, i in enumerate(np.arange(1,max(Thalamus_CIs)+1,1)):
	mean_ratio_per_thalamus_parcel[ix] = np.nanmean(thalamus_between_ratio[Thalamus_CIs==i])

#write to nifti
thalamus_between_ratio = thalamus_between_ratio*100
atlas_path = path_to_ROIs+'/Thalamus_indices.nii.gz' 
image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/MGH_%s_nodal_role_bcr.nii.gz' %roi_template
make_image(atlas_path, image_path, ROI_list, pc.astype(int))




