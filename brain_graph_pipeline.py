#script to run Maxwell's brain graph algorithms
from brain_graphs import *
#from sklearn.metrics import normalized_mutual_info_score as nmi 
from FuncParcel import *

data_path='/home/despoB/connectome-thalamus/AvgMatrices/'
templates = ['Craddock_300_cortical', 'Craddock_900_cortical']
sequences = ['mx_1400', 'mx_645', 'std_2500']
min_costs = [0.02, 0.05]
min_comm_sizes = [5, 10]
min_weights = [0., 0.5, 1.0, 1.5, 2.0, 2.5 ]

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


	


