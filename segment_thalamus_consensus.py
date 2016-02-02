# segment thalamus using consensus matrices across subjects
import sys
from brain_graphs import *
from FuncParcel import *
import glob


#path to partial correlations
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

def indiv_thalamus_parcel_consensus(subject):
	''' generate each subject's individual thalamus parcel, save consensus matrix'''
	global corrmat_path, Cortical_Network_CIs_plus_thalamus, Thalamus_voxel, Cortical_Network_CIs
	fn = corrmat_path + subject + 'MGH_Gordon_333_consensus_CI_pcorr_mat'
	M = pickle.load(open(fn, "rb"))

	# winner take out ranking
	_, Thalamo_voxCIs, _, = parcel_subcortical_network(M, Cortical_Network_CIs_plus_thalamus, \
		Thalamus_voxels, Cortical_Network_CIs, Cortical_Network_CIs)

	# sort through, keep the strongest connected CI
	Tha_CIs = sort_CIs(Thalamo_voxCIs)

	# generate thalamus consensus matrix
	concensus = community_matrix(Tha_CIs ,100)
	save_path = output_path + subject +'_thalamus_consensus'
	np.savetxt(save_path, np.nan_to_num(concensus), fmt='%2.1f')


def ave_consensus(file_pattern):
	''' load all indiv concensus matrices and average across subjects'''
	global output_path
	file_path = output_path + file_pattern
	Ave_consensus = average_corrmat(file_path, np_txt=True, pickle_object=False)
	return Ave_consensus	

def group_thalamus_consensus(file_pattern):
	''' average consensus matrices across subjects, run partition on averaged matrix'''
	AveMat = ave_consensus(file_pattern)
	ave_path = output_path + 'Group_thalamus_consensus'
	np.savetxt(ave_path, AveMat)
	graph = recursive_network_partition(matrix=AveMat, min_cost=.01, max_cost=0.07, min_community_size=100 ,min_weight=0.5)
	#graph = partition_avg_costs(matrix=AveMat, costs = np.arange(0.01,0.16, 0.01), min_community_size=75 ,graph_cost=0.05)
	return graph


if __name__ == "__main__":
	##run subject consensus
	#subject = sys.stdin.read().strip('\n')
	#indiv_thalamus_parcel_consensus(subject)

	## cluster group results
	#MGH
	#file_pattern = 'MGH*consensus'
	#AveMat = ave_consensus(file_pattern)
	#np.savetxt(Parcel_path + '/Group_thalamus_consensus', AveMat)

	#NKI
	#file_pattern = 'NKI*consensus'
	#AveMat = ave_consensus(file_pattern)
	#np.savetxt(Parcel_path + '/NKI_thalamus_consensus', AveMat)

	# write out clustering results
	#for c in np.arange(0.05, 0.16, 0.1):
	AveMat = np.loadtxt(output_path +'Group_thalamus_consensus')
	#graph = matrix_to_igraph(AveMat,cost=c)
	graph = recursive_network_partition(matrix=AveMat, min_cost=.07, max_cost=0.15, min_community_size=40 ,min_weight=0.1)
	#graph = partition_avg_costs(matrix=AveMat, costs = np.arange(0.01,0.15, 0.01), min_community_size=50 ,graph_cost=0.05)
	save_object(graph, Parcel_path +'/MGH_tha_consensus_parcel_graph')
	Tha_CI = np.array(graph.community.membership)+1
	# for i in Counter(Tha_CI):
	# 	if Counter(Tha_CI)[i] < 40:
	# 		Tha_CI[Tha_CI==i] = 0

	print(Counter(Tha_CI))
	np.savetxt(Parcel_path + '/MGH_thalamus_clusters_c07-15', Tha_CI)
	atlas_path = path_to_ROIs+'/Thalamus_indices.nii.gz' 
	ROI_list = path_to_ROIs + '/thalamus_voxel_indices' 
	image_path = Parcel_path + '/MGH_gordon_consensus_based_thalamus_parcels_c07-15.nii.gz' 
	make_image(atlas_path, image_path, ROI_list, Tha_CI)

	#repeat for NKI
	AveMat = np.loadtxt(output_path +'NKI_thalamus_consensus')
	#graph = matrix_to_igraph(AveMat,cost=c)
	graph = recursive_network_partition(matrix=AveMat, min_cost=.07, max_cost=0.15, min_community_size=40 ,min_weight=0.1)
	save_object(graph, Parcel_path +'/NKI_tha_consensus_parcel_graph')
	#graph = partition_avg_costs(matrix=AveMat, costs = np.arange(0.01,0.15, 0.01), min_community_size=50 ,graph_cost=0.05)
	#save_object(graph, Parcel_path +'/Tha_consensus_parc_graph')
	Tha_CI = np.array(graph.community.membership)+1
	# for i in Counter(Tha_CI):
	# 	if Counter(Tha_CI)[i] < 40:
	# 		Tha_CI[Tha_CI==i] = 0

	print(Counter(Tha_CI))
	np.savetxt(Parcel_path + '/NKI_thalamus_clusters_c07-15', Tha_CI)
	atlas_path = path_to_ROIs+'/Thalamus_indices.nii.gz' 
	ROI_list = path_to_ROIs + '/thalamus_voxel_indices' 
	image_path = Parcel_path + '/NKI_gordon_consensus_based_thalamus_parcels_c07-15.nii.gz' 
	make_image(atlas_path, image_path, ROI_list, Tha_CI)

	image_path = Parcel_path + '/MGH_WTA_parcel.nii.gz' 

