
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
Cog_component = np.loadtxt(Parcel_path+'yeo_flex')
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
Thalamus_df = pd.DataFrame(columns=('Voxel', 'Functional Network', 'Anatomical Parcellations', 'Morel Parcellation', 'Classficiation')) 
Thalamus_df['Voxel'] = Thalamus_voxels

Thalamus_df['Functional Network'] = Thalamus_CIs
Thalamus_df['Functional Network'].loc[Thalamus_df['Functional Network'] ==1] = 'DF'
Thalamus_df['Functional Network'].loc[Thalamus_df['Functional Network'] ==2] = 'CO'
Thalamus_df['Functional Network'].loc[Thalamus_df['Functional Network'] ==3] = 'SM'
Thalamus_df['Functional Network'].loc[Thalamus_df['Functional Network'] ==4] = 'FP'
Thalamus_df['Functional Network'].loc[Thalamus_df['Functional Network'] ==5] = 'latO'
Thalamus_df['Functional Network'].loc[Thalamus_df['Functional Network'] ==6] = 'mO'
Thalamus_df['Functional Network'].loc[Thalamus_df['Functional Network'] ==7] = 'mT'
Thalamus_df['Functional Network'].loc[Thalamus_df['Functional Network'] ==8] = 'T'
Thalamus_df['Functional Network'].loc[Thalamus_df['Functional Network'] ==9] = 'sFP'
Thalamus_df['Functional Network'].loc[Thalamus_df['Functional Network'] ==10] = 'T'

Thalamus_df['Anatomical Parcellations'] = Thalamus_ana_parcel
Thalamus_df['Anatomical Parcellations'].loc[Thalamus_df['Anatomical Parcellations'] ==1] = 'M'
Thalamus_df['Anatomical Parcellations'].loc[Thalamus_df['Anatomical Parcellations'] ==2] = 'S'
Thalamus_df['Anatomical Parcellations'].loc[Thalamus_df['Anatomical Parcellations'] ==3] = 'O'
Thalamus_df['Anatomical Parcellations'].loc[Thalamus_df['Anatomical Parcellations'] ==4] = 'PFC'
Thalamus_df['Anatomical Parcellations'].loc[Thalamus_df['Anatomical Parcellations'] ==5] = 'pM'
Thalamus_df['Anatomical Parcellations'].loc[Thalamus_df['Anatomical Parcellations'] ==6] = 'PL'
Thalamus_df['Anatomical Parcellations'].loc[Thalamus_df['Anatomical Parcellations'] ==7] = 'T'


Thalamus_df['Morel Parcellations'] = Morel_parcel
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==0] = 'Unclassified'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==1] = 'AN'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==2] = 'AN'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==3] = 'AN'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==4] = 'LD'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==5] = 'MD'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==6] = 'MD'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==7] = 'MV'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==8] = 'CL'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==9] = 'CeM'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==10] = 'CM'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==11] = 'Pf'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==12] = 'PuM'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==13] = 'PuI'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==14] = 'PuL'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==15] = 'PuA'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==16] = 'LP'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==17] = 'MGN'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==18] = 'SG'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==19] = 'Li'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==20] = 'Po'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==21] = 'LGN'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==22] = 'LGN'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==23] = 'VPL'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==24] = 'VPL'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==25] = 'VPM'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==26] = 'VPI'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==27] = 'VL'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==28] = 'VL'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==29] = 'VL'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==30] = 'VA'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==31] = 'VA'
Thalamus_df['Morel Parcellations'].loc[Thalamus_df['Morel Parcellations'] ==32] = 'VM'


Thalamus_df['Classification'] = Morel_parcel
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==0] = 'Unclassified'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==1] = 'Higher Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==2] = 'Higher Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==3] = 'Higher Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==4] = 'Higher Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==5] = 'Higher Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==6] = 'Higher Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==7] = 'Higher Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==8] = 'Nonspecific \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==9] = 'Nonspecific \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==10] = 'Nonspecific \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==11] = 'Nonspecific \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==12] = 'Higher Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==13] = 'Higher Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==14] = 'Higher Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==15] = 'Higher Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==16] = 'Higher Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==17] = 'First Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==18] = 'First Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==19] = 'Nonspecific \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==20] = 'Higher Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==21] = 'First Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==22] = 'First Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==23] = 'First Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==24] = 'First Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==25] = 'First Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==26] = 'Higher Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==27] = 'First Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==28] = 'First Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==29] = 'First Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==30] = 'First Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==31] = 'First Order \nThalamic Nuclei'
Thalamus_df['Classification'].loc[Thalamus_df['Classification'] ==32] = 'First Order \nThalamic Nuclei'

#cortical
Cortical_df = pd.DataFrame(columns=('ROI', 'Functional Network'))
#Cortical_df['ROI'] = Cortical_ROIs
Cortical_df['Functional Network'] = Cortical_CI

Cortical_df['Functional Network'].loc[Cortical_df['Functional Network'] ==1] = 'DF'
Cortical_df['Functional Network'].loc[Cortical_df['Functional Network'] ==2] = 'mO'
Cortical_df['Functional Network'].loc[Cortical_df['Functional Network'] ==3] = 'SM'
Cortical_df['Functional Network'].loc[Cortical_df['Functional Network'] ==4] = 'FP'
Cortical_df['Functional Network'].loc[Cortical_df['Functional Network'] ==5] = 'sFP'
Cortical_df['Functional Network'].loc[Cortical_df['Functional Network'] ==6] = 'mT'
Cortical_df['Functional Network'].loc[Cortical_df['Functional Network'] ==7] = 'T'
Cortical_df['Functional Network'].loc[Cortical_df['Functional Network'] ==8] = 'latO'
Cortical_df['Functional Network'].loc[Cortical_df['Functional Network'] ==9] = 'CO'
Cortical_df['Functional Network'].loc[Cortical_df['Functional Network'] ==0] = 'Other' 

##########################################
###### populate data
##########################################


#nodal roles
Thalamus_df['PC'] = Tha_PC[320:]/13.5  #take out cortical ones, normalize by max PC
Thalamus_df['WMD'] = Tha_WMD[320:]/15 
Thalamus_df['NNC'] = NNCs[320:]/15 
Thalamus_df['BNWR'] = Tha_BNWR[320:]/15
Thalamus_df['cog'] = Cog_component

Cortical_df['PC'] = Cortical_PC[0:320]/13.5  #take out cortical ones, normalize by max PC
Cortical_df['WMD'] = Cortical_WMD[0:320]/15 
Cortical_df['NNC'] = NNCs[0:320]/15 
Cortical_df['BNWR'] = Cortical_BNWR[0:320]/15

Cortical_df['Classification'] = 'Cortical \nNone Hubs'
Cortical_df['Classification'].loc[Cortical_df['PC'] > .63] = 'Cortical \nConnector Hubs'
Cortical_df['Classification'].loc[Cortical_df['WMD'] > .63] = 'Cortical \nProvincial Hubs'
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
#Thalamus_df.groupby(['Morel Parcellations']).mean()
total_df=pd.concat([Thalamus_df, Cortical_df])
total_df.to_csv(path_to_data_folder+'Cortical_plus_thalamus_nodal_WTA.csv')
