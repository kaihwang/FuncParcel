#scrip to analyze thalamus nodal proerty

from brain_graphs import *
from FuncParcel import *
import matplotlib.pylab as plt


################################################################
###### setup
################################################################

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



################################################################
###### find cortical target and nontarget, save output 
################################################################

# load MGH and NKI
MGH_thalamocor_adj = np.loadtxt(Parcel_path+'/MGH_Craddock_300_cortical_plus_thalamus_parcorrmatavg')
NKI_thalamocor_adj = np.loadtxt(Parcel_path+'/NKI_mx_645_Craddock_300_cortical_plus_thalamus_parcorrmatavg')

MGH_Cortical_targets, MGH_Cortical_nontargets, cost_thresholds = map_subcortical_cortical_targets(MGH_thalamocor_adj, Cortical_ROIs, Thalamus_voxels)
NKI_Cortical_targets, NKI_Cortical_nontargets, cost_thresholds = map_subcortical_cortical_targets(NKI_thalamocor_adj, Cortical_ROIs, Thalamus_voxels)

# lets intersect the two dataset to find reliable targets... and non-targets....
Cortical_targets={}
Cortical_nontargets ={}
for i in Thalamus_voxels:
	Cortical_targets[i] = np.intersect1d(NKI_Cortical_targets[15,i], MGH_Cortical_targets[15,i])
	Cortical_nontargets[i] = np.intersect1d(NKI_Cortical_nontargets[15,i], MGH_Cortical_nontargets[15,i])


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
image_path = Parcel_path +'/MGH+NKI_cortical_roi_target_count.nii.gz' 
make_image(atlas_path, image_path, Cortical_ROIs, Cortical_ROI_target_count)


################################################################
###### do network based parcellation
################################################################

#call function
Netadjmat = np.loadtxt('/home/despoB/connectome-thalamus/Thalamic_parcel/NKI_mx_645_cortical_network_plus_thalamus_parcorrmatavg')
Cortical_Network_CIs = np.loadtxt('/home/despoB/connectome-thalamus/ROIs/Cortical_Network_CIs')
Cortical_Network_CIs_plus_thalamus = np.loadtxt('/home/despoB/connectome-thalamus/ROIs/Cortical_CI_plus_thalamus_ROIs')
#path_to_Cortical_CI = '/home/despoB/connectome-thalamus/ROIs/Cortical_Network_CIs'


#clip out bad SNR networks (orbital frontal and inf temporal)
Netadjmat[8,:]=0
Netadjmat[:,8]=0
Netadjmat[9,:]=0
Netadjmat[:,9]=0

#call function
_, Thalamo_ParcelCIs, _, = parcel_subcortical_network(Netadjmat, Cortical_Network_CIs_plus_thalamus, \
	Thalamus_voxels, Cortical_Network_CIs, Cortical_Network_CIs)
save_object(Thalamo_ParcelCIs, path_to_data_folder +'/Thalamo_ParcelCIs')


#sort CI vector
Thalamus_CIs = np.zeros(len(Thalamus_voxel_coordinate))
for i, thalamus_voxel_index in enumerate(Thalamus_voxel_coordinate[:,3]):
	Thalamus_CIs[i] = Thalamo_ParcelCIs[thalamus_voxel_index][0]

# show histogram of CI distribution
Thalamus_CIs = Thalamus_CIs.astype(int)
plt.hist(Thalamus_CIs, 11)
plt.show()

#write out parcel info into nifti
atlas_path = path_to_ROIs+'/Thalamus_indices.nii.gz' 
ROI_list = path_to_ROIs + '/thalamus_voxel_indices' 
image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/MGH_network_based_thalamus_parcel.nii.gz' 
make_image(atlas_path, image_path, ROI_list, Thalamus_CIs)

save_object(Thalamus_CIs, path_to_data_folder +'/Thalamus_CIs')


################################################################
###### Here analyze each thalamic voxel's nodal property
################################################################
Cortical_plus_thalamus_CI = np.concatenate((Cortical_CI, Thalamus_CIs), 1) 
Cortical_plus_thalamus_CI = Cortical_plus_thalamus_CI.astype(int)

save_object(Cortical_plus_thalamus_CI, path_to_data_folder +'/Cortical_plus_thalamus_CI')
save_object(Thalamus_CIs, path_to_data_folder +'/Thalamus_CIs')

# do nodal roles on MGH and NKI
MGH_thalamocor_adj = np.loadtxt(Parcel_path+'/MGH_Craddock_300_cortical_plus_thalamus_parcorrmatavg')
NKI_thalamocor_adj = np.loadtxt(Parcel_path+'/NKI_mx_645_Craddock_300_cortical_plus_thalamus_parcorrmatavg')

MGH_cor_adj = np.loadtxt('/home/despoB/connectome-thalamus/AvgMatrices/MGH/_Craddock_300_cortical_corrmatavg')
NKI_cor_adj = np.loadtxt('/home/despoB/connectome-thalamus/AvgMatrices/NKI/_Craddock_300_cortical__mx_645_corrmatavg')

##cal_thalamus_and_cortical_ROIs_nodal_properties(MGH_thalamocor_adj, MGH_cor_adj, 'MGH')
#cal_thalamus_and_cortical_ROIs_nodal_properties(NKI_thalamocor_adj, NKI_cor_adj, 'NKI')
def cal_thalamus_and_cortical_ROIs_nodal_properties(Thalamocortical_corrmat, Cortical_adj, dset):
	global Cortical_plus_thalamus_CI, cost_thresholds, Cortical_ROIs_positions, Cortical_CI, path_to_ROIs, path_to_data_folder
	global Thalamus_voxel_positions, Thalamus_CIs, ROIs, Cortical_ROIs 

	## setup
	#Thalamocortical_corrmat = MGH_thalamocor_adj
	#Cortical_adj = MGH_cor_adj
	#dset = 'MGH'
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
	image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/%s_tha_nodal_role_pc.nii.gz' %dset
	make_image(atlas_path, image_path, ROIs, Tha_PCs_percentage)
	file_path = path_to_data_folder +'/%s_Tha_PCs' %dset
	save_object(Tha_PCs, file_path) 

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
	image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/%s_tha_nodal_role_bnwr.nii.gz' %dset
	make_image(atlas_path, image_path, ROIs, Tha_BNWR_percentage)
	file_path = path_to_data_folder +'/%s_Tha_BNWR' %dset
	save_object(Tha_BNWR, file_path) %dset


	#get cortical PC and BNEW and WMD
	Cortical_adj[np.isnan(Cortical_adj)] = 0

	Cortical_wm_mean = {}
	Cortical_wm_std = {}
	Cortical_PCs = np.zeros(Cortical_CI.size)
	Cortical_WMDs = np.zeros(Cortical_CI.size)
	Cortical_BNWR = np.empty(Cortical_CI.size)
	for ix, c in enumerate(np.arange(0.01,0.16, 0.01)):

		#pc
		#graph = matrix_to_igraph(Cortical_adj.copy(),cost=c)
		#nodal_graph = brain_graph(VertexClustering(graph, membership=Cortical_CI))
		Cortical_PCs += bct.participation_coef(bct.threshold_proportional(Cortical_adj, c, copy=True), Cortical_CI)

		#wmd
		M = bct.weight_conversion(bct.threshold_proportional(Cortical_adj, c, copy=True), 'binarize')
		Cortical_WMDs += bct.module_degree_zscore(M, Cortical_CI)
		#return mean and degree 
		for CI in np.unique(Cortical_CI):
			Cortical_wm_mean[ix+1, CI] = np.nanmean(np.sum(M[Cortical_CI==CI,:],1))
			Cortical_wm_std[ix+1, CI] = np.nanstd(np.sum(M[Cortical_CI==CI,:],1))

		#BNWR
		M = bct.threshold_proportional(Cortical_adj, c, copy=True)
		BNWR = np.zeros(Cortical_CI.size)
		for i in range(len(Cortical_CI)):
			sum_between_weight = np.nansum(M[i,Cortical_CI!=Cortical_CI[i]])
			sum_total = np.nansum(M[i,:])
			BNWR[i] = sum_between_weight / sum_total
			BNWR[i] = np.nan_to_num(BNWR[i])
			#sBNWR = np.concatenate(BNWR,sBNWR)
		Cortical_BNWR += BNWR	


	Cortical_PCs_percentage = (Cortical_PCs/13.5)*100 #13.5 is the highest PC possible  (upper bound is .9 for each cost, .9*15 = 13.5)
	atlas_path = path_to_ROIs+'/Craddock_300_cortical.nii.gz' 
	image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/%s_cortical_nodal_role_pc.nii.gz' %dset
	make_image(atlas_path, image_path, Cortical_ROIs, Cortical_PCs_percentage)
	file_path = path_to_data_folder +'/%s_Cortical_PCs' %dset
	save_object(Cortical_PCs, file_path) 

	Cortical_WMDs_percentage = (Cortical_WMDs/15)*100 
	image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/%s_cortical_nodal_role_wmd.nii.gz' %dset
	make_image(atlas_path, image_path, Cortical_ROIs, Cortical_WMDs_percentage)
	file_path = path_to_data_folder +'/%s_Cortical_WMDs' %dset
	save_object(Cortical_WMDs, file_path) 

	Cortical_BNWR_percentage = (Cortical_BNWR/15)*100 
	image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/%s_cortical_nodal_role_BNWR.nii.gz' %dset
	make_image(atlas_path, image_path, Cortical_ROIs, Cortical_BNWR_percentage)
	file_path = path_to_data_folder +'/%s_Cortical_BNWR' %dset
	save_object(Cortical_BNWR, file_path) 

	#get thlamus WMD, using mean and STD from cortex. 
	Tha_WMDs = np.empty(Cortical_plus_thalamus_CI.size)
	for ix, c in enumerate(cost_thresholds):
		Par_adj = Thalamocortical_corrmat.copy()
		Par_adj[Cortical_ROIs_positions[Cortical_CI==0],:]=0
		Par_adj[:,Cortical_ROIs_positions[Cortical_CI==0]]=0
		Par_adj[Par_adj<c]=0
		M = bct.weight_conversion(Par_adj, 'binarize')

		tha_wmd = np.zeros(Cortical_plus_thalamus_CI.size)
		for i in np.unique(Cortical_CI):
			tha_wmd[Cortical_plus_thalamus_CI==i] = (np.sum(M[Cortical_plus_thalamus_CI==i,0:320],1) - Cortical_wm_mean[ix+1,i])/Cortical_wm_std[ix+1,i]
		tha_wmd = np.nan_to_num(tha_wmd)
		Tha_WMDs += tha_wmd
		#graph = matrix_to_igraph(Par_adj.copy(),cost=1)
		#nodal_graph = brain_graph(VertexClustering(graph, membership=Cortical_plus_thalamus_CI))
		#save_path = '/home/despoB/kaihwang/bin/FuncParcel/Data/Graph/graph_MGH_%s' %roi_template
		#save_object(nodal_graph, save_path)

	Tha_WMDs[Cortical_ROIs_positions.astype('int')]=0
	Tha_WMDs_percentage = (Tha_WMDs/15)*100 #13.5 is the highest PC possible  (upper bound is .9 for each cost, .9*15 = 13.5)
	atlas_path = path_to_ROIs+'/Craddock_300_plus_thalamus_ROIs.nii.gz' 
	image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/%s_tha_nodal_role_WMD.nii.gz' %dset
	make_image(atlas_path, image_path, ROIs, Tha_WMDs_percentage)
	file_path = path_to_data_folder + '/%s_Tha_WMDs' %dset
	save_object(Tha_WMDs, file_path)


	#get number of networks/communities connected
	NNCs = np.empty(Cortical_plus_thalamus_CI.size)
	for ic, c in enumerate(cost_thresholds):
		Par_adj = Thalamocortical_corrmat.copy()
		#Par_adj[Cortical_ROIs_positions[Cortical_CI==0],:]=0
		#Par_adj[:,Cortical_ROIs_positions[Cortical_CI==0]]=0
		Par_adj[Par_adj<c]=0

		Tha_NNCs = np.zeros(Cortical_plus_thalamus_CI.size)
		for ix, i in enumerate(Thalamus_voxel_positions):
			Tha_NNCs[i] = len(np.unique(Cortical_plus_thalamus_CI[Par_adj[i,]!=0]))

		NNCs += Tha_NNCs


	Cortical_NNCs = np.empty(Cortical_plus_thalamus_CI.size)
	for ic, c in enumerate(np.arange(0.01,0.16, 0.01)):
		M= bct.threshold_proportional(Cortical_adj, c, copy=True)

		Cor_NNCs = np.zeros(Cortical_plus_thalamus_CI.size)
		for ix, i in enumerate(Cortical_ROIs_positions):
			Cor_NNCs[i] = len(np.unique(Cortical_plus_thalamus_CI[M[i,]!=0]))

		Cortical_NNCs += Cor_NNCs
	NNCs[Cortical_ROIs_positions] = Cortical_NNCs[Cortical_ROIs_positions] 

	NNCs_percentage = (NNCs/15)*100 #13.5 is the highest PC possible  (upper bound is .9 for each cost, .9*15 = 13.5)
	atlas_path = path_to_ROIs+'/Craddock_300_plus_thalamus_ROIs.nii.gz' 
	image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/%s_cortex_tha_nodal_role_NNC.nii.gz' %dset
	make_image(atlas_path, image_path, ROIs, NNCs_percentage)
	file_path = path_to_data_folder + '/%s_NNCs' %dset
	save_object(NNCs, file_path) 



################################################################
#### Thalamus's nodal properties based on within thalamus weight
################################################################

Thalamus_corrmat = np.loadtxt('/home/despoB/connectome-thalamus/AvgMatrices/MGH/_Craddock_300_plus_thalamus_ROIs_ncsreg_corrmatavg')
Thalamus_corrmat = Thalamus_corrmat[Thalamus_voxel_positions,:][:,Thalamus_voxel_positions]

#get PC and WMD
within_Tha_PCs = np.empty(Thalamus_voxel_positions.size)
within_Tha_WMDs = np.empty(Thalamus_voxel_positions.size)

for c in np.arange(0.01,0.16,0.01):
	
	graph = matrix_to_igraph(Thalamus_corrmat.copy(),cost=c)
	nodal_graph = brain_graph(VertexClustering(graph, membership=Thalamus_CIs))
	tmp_pc = np.nan_to_num(nodal_graph.pc)
	tmp_wmd = np.nan_to_num(nodal_graph.wmd)
	within_Tha_PCs += tmp_pc
	within_Tha_WMDs += tmp_wmd
	#save_path = '/home/despoB/kaihwang/bin/FuncParcel/Data/Graph/graph_MGH_%s' %roi_template
	#save_object(nodal_graph, save_path)

within_Tha_PCs_percentage = (within_Tha_PCs/13.5)*100 
atlas_path = path_to_ROIs+'/Thalamus_indices.nii.gz' 
image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/MGH_within_tha_nodal_role_pc.nii.gz'
make_image(atlas_path, image_path, Thalamus_voxels, within_Tha_PCs_percentage)
save_object(within_Tha_PCs, path_to_data_folder +'/within_Tha_PCs')

within_Tha_WMDs_percentage = (within_Tha_WMDs/13.5)*100 #13.5 is the highest PC possible  (upper bound is .9 for each cost, .9*15 = 13.5)
atlas_path = path_to_ROIs+'/Thalamus_indices.nii.gz' 
image_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/MGH_within_tha_nodal_role_wmd.nii.gz'
make_image(atlas_path, image_path, Thalamus_voxels, within_Tha_WMDs_percentage)
save_object(within_Tha_WMDs, path_to_data_folder +'/within_Tha_WMDs')



################################################################
#### get descriptive stats of nodal rols per parcel
################################################################
mean_PC_per_thalamus_parcel = np.zeros(max(Thalamus_CIs))
for ix, i in enumerate(np.unique(Thalamus_CIs)):
	mean_PC_per_thalamus_parcel[ix] = np.nanmean(Tha_PCs_percentage[Thalamus_CIs==i])

mean_WMD_per_thalamus_parcel = np.zeros(max(Thalamus_CIs))
for ix, i in enumerate(np.unique(Thalamus_CIs)):
	mean_WMD_per_thalamus_parcel[ix] = np.nanmean(Tha_WMDs_percentage[Thalamus_CIs==i])	

mean_BNWR_per_thalamus_parcel = np.zeros(max(Thalamus_CIs))
for ix, i in enumerate(np.unique(Thalamus_CIs)):
	mean_BNWR_per_thalamus_parcel[ix] = np.nanmean(Tha_BNWR_percentage[Thalamus_CIs==i])

Tha_NNCs = NNCs_percentage[Thalamus_voxel_positions]/100		
mean_NNCs_per_thalamus_parcel = np.zeros(max(Thalamus_CIs))
for ix, i in enumerate(np.unique(Thalamus_CIs)):
	mean_NNCs_per_thalamus_parcel[ix] = np.nanmean(Tha_NNCs[Thalamus_CIs==i])


mean_PC_per_Cortical_parcel = np.zeros(len(np.unique(Cortical_CI)))
for ix, i in enumerate(np.unique(Cortical_CI)):
	mean_PC_per_Cortical_parcel[ix] = np.nanmean(Cortical_PCs_percentage[Cortical_CI==i])

mean_WMD_per_Cortical_parcel = np.zeros(len(np.unique(Cortical_CI)))
for ix, i in enumerate(np.unique(Cortical_CI)):
	mean_WMD_per_Cortical_parcel[ix] = np.nanmean(Cortical_WMDs_percentage[Cortical_CI==i])

mean_BNWR_per_Cortical_parcel = np.zeros(len(np.unique(Cortical_CI)))
for ix, i in enumerate(np.unique(Cortical_CI)):
	mean_BNWR_per_Cortical_parcel[ix] = np.nanmean(Cortical_BNWR_percentage[Cortical_CI==i])

Cortical_NNCs = NNCs_percentage[Cortical_ROIs_positions]/100
mean_NNCs_per_Cortical_parcel = np.zeros(len(np.unique(Cortical_CI)))
for ix, i in enumerate(np.unique(Cortical_CI)):
	mean_NNCs_per_Cortical_parcel[ix] = np.nanmean(Cortical_NNCs[Cortical_CI==i])




