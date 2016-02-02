
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
Thalamus_ana_parcel = np.loadtxt(Parcel_path+'fsl_thalamus_ana_parcel') # this is diffusion based parcel from FSL
Morel_parcel = np.loadtxt(Parcel_path+'Morel_parcel') # this is histology based parcel of the Morel atlas
Cortical_ROIs = np.loadtxt(path_to_ROIs+'Craddock_300_cortical_ROIs', dtype=int)
Cortical_CI = np.loadtxt(path_to_ROIs+'/Cortical_CI', dtype='int')
#Cortical_CI = int(np.array(Cortical_CI))

#PC
Tha_PC = pickle.load(open(path_to_data_folder+'MGH_Tha_PCs', "rb"))
Cortical_PC = pickle.load(open(path_to_data_folder+'Cortical_PCs', "rb"))

#NNCs
NNCs = pickle.load(open(path_to_data_folder+'MGH_NNCs', "rb"))

#WMDs
Tha_WMD = pickle.load(open(path_to_data_folder+'MGH_Tha_WMDs', "rb"))
Cortical_WMD = pickle.load(open(path_to_data_folder+'Cortical_WMDs', "rb"))

#BNWR
Tha_BNWR = pickle.load(open(path_to_data_folder+'MGH_Tha_BNWR', "rb"))
Cortical_BNWR = pickle.load(open(path_to_data_folder+'Cortical_BNWR', "rb"))

#BCC
#bcc = pickle.load(open(path_to_graph+'MGH_bcc', "rb"))


################################################################
###### create Dataframe
################################################################

#thalamus parcellations
Thalamus_df = pd.DataFrame(columns=('Voxel', 'Associated System', 'Anatomical Parcellations')) 
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

Thalamus_df['Anatomical Parcellations'] = Thalamus_ana_parcel
Thalamus_df['Anatomical Parcellations'].loc[Thalamus_df['Anatomical Parcellations'] ==1] = 'Primary Motor'
Thalamus_df['Anatomical Parcellations'].loc[Thalamus_df['Anatomical Parcellations'] ==2] = 'Sensory'
Thalamus_df['Anatomical Parcellations'].loc[Thalamus_df['Anatomical Parcellations'] ==3] = 'Occipital'
Thalamus_df['Anatomical Parcellations'].loc[Thalamus_df['Anatomical Parcellations'] ==4] = 'Prefrontal'
Thalamus_df['Anatomical Parcellations'].loc[Thalamus_df['Anatomical Parcellations'] ==5] = 'Premotor'
Thalamus_df['Anatomical Parcellations'].loc[Thalamus_df['Anatomical Parcellations'] ==6] = 'Parietal'
Thalamus_df['Anatomical Parcellations'].loc[Thalamus_df['Anatomical Parcellations'] ==7] = 'Temporal'


Thalamus_df['Morel Parcellations'] = Morel_parcel
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==0] = 'Unclassified'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==10] = 'AD'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==11] = 'AM'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==12] = 'AV'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==13] = 'LD'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==20] = 'MD'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==21] = 'MD'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==22] = 'MV'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==23] = 'CL'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==24] = 'CeM'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==25] = 'CM'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==26] = 'Pf'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==30] = 'PuM'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==31] = 'PuI'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==32] = 'PuL'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==33] = 'PuA'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==34] = 'LP'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==35] = 'MGN'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==36] = 'SG'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==37] = 'Li'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==38] = 'Po'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==39] = 'LGN'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==40] = 'LGN'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==41] = 'VPL'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==42] = 'VPL'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==43] = 'VPM'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==44] = 'VPI'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==45] = 'VL'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==46] = 'VL'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==47] = 'VL'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==48] = 'VA'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==49] = 'VA'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==50] = 'VM'

#cortical
Cortical_df = pd.DataFrame(columns=('ROI', 'Associated System'))
#Cortical_df['ROI'] = Cortical_ROIs
Cortical_df['Associated System'] = Cortical_CI

Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==1] = 'Default'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==2] = 'Visual'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==3] = 'Somatomotor'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==4] = 'Fronto-parietal'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==5] = 'Attention'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==6] = 'Retrospinal'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==7] = 'Temporal'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==8] = 'Occipital-parietal'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==9] = 'Cingulo-opercular'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==0] = 'Other' 

##########################################
###### populate data
##########################################


#nodal roles
Thalamus_df['PC'] = Tha_PC[320:]/13.5  #take out cortical ones, normalize by max PC
Thalamus_df['WMD'] = Tha_WMD[320:]/15 
Thalamus_df['NNC'] = NNCs[320:]/15 
Thalamus_df['BNWR'] = Tha_BNWR[320:]/15

Cortical_df['PC'] = Cortical_PC[0:320]/13.5  #take out cortical ones, normalize by max PC
Cortical_df['WMD'] = Cortical_WMD[0:320]/15 
Cortical_df['NNC'] = NNCs[0:320]/15 
Cortical_df['BNWR'] = Cortical_BNWR[0:320]/15



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



