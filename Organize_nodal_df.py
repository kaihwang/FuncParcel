
from FuncParcel import *
import pandas as pd
#from ggplot import *
import pickle
import matplotlib.pyplot as plt
from scipy.stats.mstats import zscore as zscore

################################################################
###### setup
################################################################
Parcel_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/'
path_to_data_folder = '/home/despoB/kaihwang/bin/FuncParcel/Data/'
path_to_graph = '/home/despoB/kaihwang/Rest/Graph/'
path_to_ROIs = '/home/despoB/connectome-thalamus/ROIs/'


#### load data
Thalamus_CIs = np.loadtxt(Parcel_path+'MGH_Thalamus_WTA_CIs')
Thalamus_voxels = np.loadtxt(path_to_ROIs+'thalamus_voxel_indices', dtype=int)
Cortical_ROIs = np.loadtxt(path_to_ROIs+'Gordon_333', dtype=int)
Cortical_CI = Cortical_CI = np.loadtxt(path_to_ROIs+'/Gordon_consensus_CI', dtype='int')

#PC
Tha_PC = pickle.load(open(path_to_graph+'MGH_avemat_tha_nodal_pcorr_PCs', "rb"))
Cortical_PC = pickle.load(open(path_to_graph+'MGH_avemat_tha_nodal_pcorr_PCs', "rb"))

#NNCs
NNCs = pickle.load(open(path_to_graph+'MGH_avemat_tha_nodal_pcorr_NNCs', "rb"))

#WMDs
Tha_WMD = pickle.load(open(path_to_graph+'MGH_avemat_tha_nodal_pcorr_WMDs', "rb"))
Cortical_WMD = pickle.load(open(path_to_graph+'MGH_avemat_tha_nodal_pcorr_WMDs', "rb"))

#BNWR
Tha_BNWR = pickle.load(open(path_to_graph+'MGH_avemat_tha_nodal_pcorr_BNWRs', "rb"))
Cortical_BNWR = pickle.load(open(path_to_graph+'MGH_avemat_tha_nodal_pcorr_BNWRs', "rb"))

#BCC
#bcc = pickle.load(open(path_to_graph+'MGH_bcc', "rb"))


################################################################
###### create Dataframe
################################################################
#thalamus
Thalamus_df = pd.DataFrame(columns=('Voxel', 'Associated System')) 
Thalamus_df['Voxel'] = Thalamus_voxels
Thalamus_df['Associated System'] = Thalamus_CIs

Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==1] = 'Default'
Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==2] = 'Cingulo-opercular'
Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==3] = 'Somatomotor'
Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==4] = 'Fronto-parietal'
Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==5] = 'Occipital-Parietal'
Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==6] = 'Visual'
Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==7] = 'Retrospinal'
Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==8] = 'Temporal'
Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==9] = 'Attention'
Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==10] = 'Temporal'


#cortical
Cortical_df = pd.DataFrame(columns=('ROI', 'Associated System'))
Cortical_df['ROI'] = Cortical_ROIs
Cortical_df['Associated System'] = Cortical_CI

Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==1] = 'Default'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==2] = 'Cingulo-opercular'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==3] = 'Somatomotor'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==4] = 'Fronto-parietal'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==5] = 'Occipital-Parietal'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==6] = 'Visual'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==7] = 'Retrospinal'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==8] = 'Temporal'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==9] = 'Attention'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==10] = 'Temporal'


##########################################
###### populate data
##########################################


#nodal roles
Thalamus_df['PC'] = Tha_PC[333:]  
Thalamus_df['WMD'] = Tha_WMD[333:] 
Thalamus_df['NNC'] = NNCs[333:] 
Thalamus_df['BNWR'] = Tha_BNWR[333:]


Cortical_df['PC'] = Cortical_PC[0:333]  
Cortical_df['WMD'] = Cortical_WMD[0:333] 
Cortical_df['NNC'] = NNCs[0:333]
Cortical_df['BNWR'] = Cortical_BNWR[0:333]



################################################################
####### thalamo-cortical target's nodal role for each thalamus voxel
################################################################

# Thalamus_target_df = pd.DataFrame(columns=('Voxel', 'Associated System', 'Cortical Target','Target PC','Target WMD' )) 

# Cortical_targets = pickle.load(open(path_to_data_folder +'/Cortical_targets', "rb"))

# for i, v in enumerate(Thalamus_voxels):

# 	tmp_df = pd.DataFrame(columns=('Voxel', 'Associated System', 'Cortical CI', 'Cortical Target','Target PC','Target WMD' ))
# 	tmp_PC = Cortical_PC[np.in1d(Cortical_ROIs, Cortical_targets[v])]
# 	tmp_WMD = Cortical_WMD[np.in1d(Cortical_ROIs, Cortical_targets[v])]
# 	tmp_targets = Cortical_ROIs[np.in1d(Cortical_ROIs, Cortical_targets[v])]
# 	tmp_CIs = Cortical_CI[np.in1d(Cortical_ROIs, Cortical_targets[v])]
# 	for u in range(len(tmp_PC)):
# 		tmp_df.set_value(u, 'Voxel', v)
# 		tmp_df.set_value(u, 'Associated System', Thalamus_CIs[i])
# 		tmp_df.set_value(u, 'Cortical CI', tmp_CIs[u])
# 		tmp_df.set_value(u, 'Cortical Target', tmp_targets[u])
# 		tmp_df.set_value(u, 'Target PC', tmp_PC[u]/13.5)
# 		tmp_df.set_value(u, 'Target WMD', tmp_WMD[u]/15)
# 	Thalamus_target_df = pd.concat([Thalamus_target_df, tmp_df], ignore_index=True)	

# Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==1] = 'Default'
# Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==2] = 'Visual'
# Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==3] = 'Somatomotor'
# Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==4] = 'Fronto-parietal'
# Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==5] = 'Attention'
# Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==6] = 'Cingulo-opercular'
# Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==7] = 'Temporal'
# Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==8] = 'Cingulo-parietal'
# Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==11] = 'Sailency'



################################################################
###### save
################################################################
Thalamus_df.to_csv(path_to_data_folder+'Thalamus_nodal_WTA.csv')
Cortical_df.to_csv(path_to_data_folder+'Cortical_nodal_WTA.csv')
#Thalamus_target_df.to_csv(path_to_data_folder+'Thalamus_target.csv')



