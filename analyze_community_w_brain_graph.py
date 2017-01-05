#script to run Maxwell's brain graph algorithms, for finding community structures

from brain_graphs import *
#from sklearn.metrics import normalized_mutual_info_score as nmi 
from FuncParcel import *
import glob

#data_path='/home/despoB/connectome-thalamus/AvgMatrices/'
data_path = '/home/despoB/kaihwang/Rest/NotBackedUp/AdjMatrices'
#data_path = '/home/despoB/kaihwang/Rest/NotBackedUp/WILLIAM/PartialCorrelationMatrices'
#templates = ['Craddock_300_cortical', 'Craddock_900_cortical']
templates = ['Gordon_333_cortical']
#sequences = ['mx_1400', 'mx_645', 'std_2500']
#min_costs = [0.02, 0.05]
#min_comm_sizes = [5, 10]
#min_weights = [0., 0.5, 1.0, 1.5, 2.0, 2.5 ]
output_path = '/home/despoB/connectome-thalamus/Graph/'


def cluster_on_grouped_ave_matrices(templates, min_costs, min_comm_sizes, min_weights, data_path):
	#this function explore a range of thresholds for clustering group averaged matrices
	for roi_template in templates:

		for mc in min_costs:
			for ms in min_comm_sizes:
				for mw in min_weights:
					cost_string = str(mc)
					size_string = str(ms)
					weight_string = str(mw)
					path = data_path + 'MGH/_%s_corrmatavg' %(roi_template)
					adj = np.loadtxt(path)
					
					# run brain graph partition
					#def partition_avg_costs(matrix,costs,min_community_size,graph_cost):
					#graph_NKI645 = recursive_network_partition(matrix=adj645, min_cost=0.05)
					#partition_avg_costs(matrix = adj, costs = np.arange(0.02, 0.11, 0.005), min_community_size=5, graph_cost=0.1)
					graph = recursive_network_partition(matrix=adj, min_cost=mc, min_community_size=ms,min_weight=mw)

					#save graph objects
					save_path = '/home/despoB/kaihwang/bin/FuncParcel/Data/graph_MGH_%s_min_cost%s_min_weight%s_min_size%s' %(roi_template, cost_string, weight_string, size_string)
					save_object(graph, save_path)
					
					#save output nifti
					atlas_path = '/home/despoB/connectome-thalamus/ROIs/%s.nii.gz' %roi_template
					ROI_list = '/home/despoB/connectome-thalamus/ROIs/%s_ROIs' %roi_template
					image_path = '/home/despoB/connectome-thalamus/ROIs/MGH_%s_CI_min_cost%s_min_weight%s_min_size%s.nii.gz' %(roi_template, cost_string, weight_string, size_string)
					make_image(atlas_path, image_path, ROI_list, graph.community.membership+1)

					for sequence in sequences:
						#load adj
						path = data_path + 'NKI/_%s__%s_corrmatavg' %(roi_template, sequence)
						adj = np.loadtxt(path)
						
						graph = recursive_network_partition(matrix=adj, min_cost=mc, min_community_size=ms,min_weight=mw)

						#save graph objects
						save_path = '/home/despoB/kaihwang/bin/FuncParcel/Data/graph_NKI_%s_%s_min_cost%s_min_weight%s_min_size%s' %(roi_template, sequence, cost_string, weight_string, size_string)
						save_object(graph, save_path)
						
						#save output nifti
						atlas_path = '/home/despoB/connectome-thalamus/ROIs/%s.nii.gz' %roi_template
						ROI_list = '/home/despoB/connectome-thalamus/ROIs/%s_ROIs' %roi_template
						image_path = '/home/despoB/connectome-thalamus/ROIs/NKI_%s_%s_CI_min_cost%s_min_weight%s_min_size%s.nii.gz' %(sequence, roi_template, cost_string, weight_string, size_string)
						make_image(atlas_path, image_path, ROI_list, graph.community.membership+1)


def gen_indiv_graph(file_pattern, output_path, prefix):
	''' generate individual graph outputs, and save consensus matrix'''

	for roi_template in templates:
		AdjMat_Files = glob.glob(file_pattern)

		for s, f in enumerate(AdjMat_Files):
			M = np.loadtxt(f)
			graph = recursive_network_partition(matrix=M, min_cost=.01,max_cost=0.15, min_community_size=5 ,min_weight=2)
			save_path = output_path + prefix + str(s) +'_graph.output'
			save_object(graph, save_path)

			#concensus = community_matrix(graph.community.membership,5)
			#save_path = output_path + 'MGH_Craddock_300_full_' + str(s) +'_consensus'
			#np.savetxt(save_path, np.nan_to_num(concensus), fmt='%2.1f')


def ave_consensus(file_pattern):
	''' load all indiv concensus matrices and average across subjects'''
	global output_path
	file_path = output_path + file_pattern
	Ave_consensus = average_corrmat(file_path, np_txt=True, pickle_object=False)
	return Ave_consensus


def compare_full_v_partial_corr(output_path):
	''' compare partial and full correlation outputs'''
	file_path = output_path + '/*Craddock_300_full*consensus'
	Full_Ave_consensus = average_corrmat(file_path, np_txt=True, pickle_object=False)

	file_path = output_path + '/*Craddock_300_partial*consensus'
	Partial_Ave_consensus = average_corrmat(file_path, np_txt=True, pickle_object=False)

	Full_graph = recursive_network_partition(matrix=Full_Ave_consensus, min_cost=.01,max_cost=0.15, min_community_size=5 ,min_weight=2)
	Partial_graph = recursive_network_partition(matrix=Partial_Ave_consensus, min_cost=.01,max_cost=0.15, min_community_size=5 ,min_weight=2)

	atlas_path = '/home/despoB/connectome-thalamus/ROIs/Craddock_300_cortical.nii.gz' 
	ROI_list = '/home/despoB/connectome-thalamus/ROIs/Craddock_300_cortical_ROIs'
	image_path = '/home/despoB/connectome-thalamus/ROIs/MGH_full_corr_CI.nii.gz' 
	make_image(atlas_path, image_path, ROI_list, np.array(Full_graph.community.membership)+1)

	image_path = '/home/despoB/connectome-thalamus/ROIs/MGH_par_corr_CI.nii.gz' 
	make_image(atlas_path, image_path, ROI_list, np.array(Partial_graph.community.membership)+1)


def cluster_group_consensus_matrices(file_pattern):
	AveMat = ave_consensus(file_pattern)
	graph = recursive_network_partition(matrix=AveMat, min_cost=.01,max_cost=0.15, min_community_size=5 ,min_weight=2)
	return graph


def write_group_consensus_cluste_results():
	''' do consensus based clustering across subjects, remove isolated clusters, write out results'''

	### cluster across averaged consensus matrices with craddock
	file_pattern = 'MGH*Craddock_300_full*consensus'
	Craddock_300_group_graph = cluster_group_consensus_matrices(file_pattern)
	Craddock_CI = np.array(Craddock_300_group_graph.community.membership)+1

	#remove isolated CIs
	Count = Counter(Craddock_CI)
	for i in Count:
		if Count[i] <5:
			Craddock_CI[Craddock_CI ==i] = 0

	#write output		
	atlas_path = '/home/despoB/connectome-thalamus/ROIs/Craddock_300_cortical.nii.gz' 
	ROI_list = '/home/despoB/connectome-thalamus/ROIs/Craddock_300_cortical_ROIs'
	image_path = '/home/despoB/connectome-thalamus/ROIs/MGH_Craddock_300_consensus_CI.nii.gz' 
	make_image(atlas_path, image_path, ROI_list, Craddock_CI)

	### cluster across averaged consensus matrices with Gordon
	file_pattern = 'MGH*Gordon*consensus'
	Gordon_group_graph = cluster_group_consensus_matrices(file_pattern)
	Gordon_CI = np.array(Gordon_group_graph.community.membership)+1
	Count = Counter(Gordon_CI)

	#remove isolated CIs
	for i in Count:
		if Count[i] <5:
			Gordon_CI[Gordon_CI ==i] = 0
	#write output
	atlas_path = '/home/despoB/connectome-thalamus/ROIs/Gordon_333_cortical.nii.gz' 
	ROI_list = '/home/despoB/connectome-thalamus/ROIs/Gordon_333'
	image_path = '/home/despoB/connectome-thalamus/ROIs/MGH_Gordon_333_consensus_CI.nii.gz' 
	make_image(atlas_path, image_path, ROI_list, Gordon_CI)


def compare_spectral_v_infomap_Q(file_pattern):
	''' compare Q values from infomap and spectral'''

	AdjMat_Files = glob.glob(file_pattern)
	Q_info = np.zeros(len(AdjMat_Files))
	Q_spectral = np.zeros(len(AdjMat_Files))
	for s, f in enumerate(AdjMat_Files):
		M = np.loadtxt(AdjMat_Files[s])
		graph = matrix_to_igraph(M.copy(), .05)
		infomap = graph.community_infomap()
		spectral = graph.community_fastgreedy().as_clustering()
		Q_info[s] = infomap.modularity
		Q_spectral[s] = spectral.modularity
	return Q_info, Q_spectral	


if __name__ == "__main__":
	#write_group_consensus_cluste_results()

	file_pattern = 'MGH*cortical*corrmat'
	Q_info, Q_spectral = compare_spectral_v_infomap_Q(file_pattern)






	


