#script to run Maxwell's brain graph algorithms
# Nov 5, here try cluster individual's, then do consensus across subjects, and cluster group results

from brain_graphs import *
#from sklearn.metrics import normalized_mutual_info_score as nmi 
from FuncParcel import *
import glob

#data_path='/home/despoB/connectome-thalamus/AvgMatrices/'
data_path = '/home/despoB/kaihwang/Rest/NotBackedUp/AdjMatrices'
#templates = ['Craddock_300_cortical', 'Craddock_900_cortical']
templates = ['Gordon_333_cortical']
#sequences = ['mx_1400', 'mx_645', 'std_2500']
#min_costs = [0.02, 0.05]
#min_comm_sizes = [5, 10]
#min_weights = [0., 0.5, 1.0, 1.5, 2.0, 2.5 ]
output_path = '/home/despoB/connectome-thalamus/Graph/'

for roi_template in templates:
	file_pattern = data_path + '/' + 'MGH*' + roi_template + '*'
	AdjMat_Files = glob.glob(file_pattern)

	for s, f in enumerate(AdjMat_Files):
		M = np.loadtxt(f)
		graph = recursive_network_partition(matrix=M, min_cost=.01,max_cost=0.15, min_community_size=5 ,min_weight=2)
		save_path = output_path + 'MGH_Gordon_333_' + str(s) +'_graph.output'
		save_object(graph, save_path)
		
		concensus = community_matrix(graph.community.membership,5)
		save_path = output_path + 'MGH_Gordon_333_' + str(s) +'_consensus'
		np.savetxt(save_path, np.nan_to_num(concensus), fmt='%2.1f')




# for roi_template in templates:

# 	for mc in min_costs:
# 		for ms in min_comm_sizes:
# 			for mw in min_weights:
# 				cost_string = str(mc)
# 				size_string = str(ms)
# 				weight_string = str(mw)
# 				path = data_path + 'MGH/_%s_corrmatavg' %(roi_template)
# 				adj = np.loadtxt(path)
				
# 				# run brain graph partition
# 				#def partition_avg_costs(matrix,costs,min_community_size,graph_cost):
# 				#graph_NKI645 = recursive_network_partition(matrix=adj645, min_cost=0.05)
# 				#partition_avg_costs(matrix = adj, costs = np.arange(0.02, 0.11, 0.005), min_community_size=5, graph_cost=0.1)
# 				graph = recursive_network_partition(matrix=adj, min_cost=mc, min_community_size=ms,min_weight=mw)

# 				#save graph objects
# 				save_path = '/home/despoB/kaihwang/bin/FuncParcel/Data/graph_MGH_%s_min_cost%s_min_weight%s_min_size%s' %(roi_template, cost_string, weight_string, size_string)
# 				save_object(graph, save_path)
				
# 				#save output nifti
# 				atlas_path = '/home/despoB/connectome-thalamus/ROIs/%s.nii.gz' %roi_template
# 				ROI_list = '/home/despoB/connectome-thalamus/ROIs/%s_ROIs' %roi_template
# 				image_path = '/home/despoB/connectome-thalamus/ROIs/MGH_%s_CI_min_cost%s_min_weight%s_min_size%s.nii.gz' %(roi_template, cost_string, weight_string, size_string)
# 				make_image(atlas_path, image_path, ROI_list, graph.community.membership+1)

# 				for sequence in sequences:
# 					#load adj
# 					path = data_path + 'NKI/_%s__%s_corrmatavg' %(roi_template, sequence)
# 					adj = np.loadtxt(path)
					
# 					graph = recursive_network_partition(matrix=adj, min_cost=mc, min_community_size=ms,min_weight=mw)

# 					#save graph objects
# 					save_path = '/home/despoB/kaihwang/bin/FuncParcel/Data/graph_NKI_%s_%s_min_cost%s_min_weight%s_min_size%s' %(roi_template, sequence, cost_string, weight_string, size_string)
# 					save_object(graph, save_path)
					
# 					#save output nifti
# 					atlas_path = '/home/despoB/connectome-thalamus/ROIs/%s.nii.gz' %roi_template
# 					ROI_list = '/home/despoB/connectome-thalamus/ROIs/%s_ROIs' %roi_template
# 					image_path = '/home/despoB/connectome-thalamus/ROIs/NKI_%s_%s_CI_min_cost%s_min_weight%s_min_size%s.nii.gz' %(sequence, roi_template, cost_string, weight_string, size_string)
# 					make_image(atlas_path, image_path, ROI_list, graph.community.membership+1)


	


