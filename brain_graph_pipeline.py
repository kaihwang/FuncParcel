#script to run Maxwell's brain graph algorithms
from brain_graphs import *
from sklearn.metrics import normalized_mutual_info_score as nmi 
from FuncParcel import *

data_path='/home/despoB/connectome-thalamus/AvgMatrices/'

#load adj
adj645 = np.loadtxt(data_path+'NKI/_Craddock_300_cortical__mx_645_corrmatavg')
adj1400 = np.loadtxt(data_path+'NKI/_Craddock_300_cortical__mx_1400_corrmatavg')
adj2500 = np.loadtxt(data_path+'NKI/_Craddock_300_cortical__std_2500_corrmatavg') 
adjMGH = np.loadtxt(data_path+'MGH/MGH_Craddock_300_cortical_corrmatavg') 

# run brain graph partition
graph_NKI645 = recursive_network_partition(matrix=adj645, min_cost=0.05)
graph_NKI1400 = recursive_network_partition(matrix=adj1400, min_cost=0.05)
graph_NKI2500 = recursive_network_partition(matrix=adj2500, min_cost=0.05)
graph_MGH = recursive_network_partition(matrix=adjMGH, min_cost=0.05)

#look at the size of communities
Counter(graph_NKI645.community.membership)
Counter(graph_NKI1400.community.membership)
Counter(graph_NKI2500.community.membership)
Counter(graph_MGH.community.membership)

#get normalized mutual information bewteen datasets (MGH v NKI)
nmi(graph_NKI645.community.membership, graph_MGH.community.membership)
nmi(graph_NKI1400.community.membership, graph_MGH.community.membership)
nmi(graph_NKI2500.community.membership, graph_MGH.community.membership)

#get NMI between sequences (for NKI)
nmi(graph_NKI645.community.membership, graph_NKI2500.community.membership)
nmi(graph_NKI1400.community.membership, graph_NKI2500.community.membership)
nmi(graph_NKI1400.community.membership, graph_NKI645.community.membership)


atlas_path = '/home/despoB/connectome-thalamus/ROIs/Craddock_300_cortical.nii.gz'
ROI_list = '/home/despoB/connectome-thalamus/ROIs/Craddock_300_cortical_ROIs'

image_path = '/home/despoB/connectome-thalamus/ROIs/NKI1400_Craddock_300_CI_mincost.05.nii.gz'
make_image(atlas_path, image_path, ROI_list, graph_NKI1400.community.membership)


image_path = '/home/despoB/connectome-thalamus/ROIs/NKI645_Craddock_300_CI_mincost.05.nii.gz'
make_image(atlas_path, image_path, ROI_list, graph_NKI645.community.membership)


image_path = '/home/despoB/connectome-thalamus/ROIs/NKI2500_Craddock_300_CI_mincost.05.nii.gz'
make_image(atlas_path, image_path, ROI_list, graph_NKI2500.community.membership)


image_path = '/home/despoB/connectome-thalamus/ROIs/MGH_Craddock_300_CI_mincost.05.nii.gz'
make_image(atlas_path, image_path, ROI_list, graph_MGH.community.membership)




