
from FuncParcel import *
import pandas as pd
from ggplot import *
import pickle
import matplotlib.pyplot as plt
from scipy.stats.mstats import zscore as zscore

################################################################
###### setup
################################################################
Parcel_path = '/home/despoB/connectome-thalamus/Thalamic_parcel/'
path_to_data_folder = '/home/despoB/kaihwang/bin/FuncParcel/Data/'
path_to_ROIs = '/home/despoB/connectome-thalamus/ROIs/'


#### load data
Thalamus_CIs = pickle.load(open(path_to_data_folder+'MGH_Thalamus_CIs', "rb"))
Thalamus_voxels = np.loadtxt(path_to_ROIs+'thalamus_voxel_indices', dtype=int)
Cortical_ROIs = np.loadtxt(path_to_ROIs+'Craddock_300_cortical_ROIs', dtype=int)
Cortical_CI = Cortical_CI = np.loadtxt(path_to_ROIs+'/Cortical_CI', dtype='int')

#PC
Tha_PC = pickle.load(open(path_to_data_folder+'Tha_PCs', "rb"))
Cortical_PC = pickle.load(open(path_to_data_folder+'Cortical_PCs', "rb"))

#NNCs
NNCs = pickle.load(open(path_to_data_folder+'NNCs', "rb"))

#WMDs
Tha_WMD = pickle.load(open(path_to_data_folder+'Tha_WMDs', "rb"))
Cortical_WMD = pickle.load(open(path_to_data_folder+'Cortical_WMDs', "rb"))

#BNWR
Tha_BNWR = pickle.load(open(path_to_data_folder+'Tha_BNWR', "rb"))
Cortical_BNWR = pickle.load(open(path_to_data_folder+'Cortical_BNWR', "rb"))


################################################################
###### create Dataframe
################################################################
#thalamus
Thalamus_df = pd.DataFrame(columns=('Voxel', 'Associated System')) 
Thalamus_df['Voxel'] = Thalamus_voxels
Thalamus_df['Associated System'] = Thalamus_CIs

Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==1] = 'Default'
Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==2] = 'Visual'
Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==3] = 'Somatomotor'
Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==4] = 'Fronto-parietal'
Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==5] = 'Attention'
Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==6] = 'Cingulo-opercular'
Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==7] = 'Temporal'
Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==8] = 'Cingulo-parietal'
Thalamus_df['Associated System'].loc[Thalamus_df['Associated System'] ==11] = 'Sailency'


#cortical
Cortical_df = pd.DataFrame(columns=('ROI', 'Associated System'))
Cortical_df['ROI'] = Cortical_ROIs
Cortical_df['Associated System'] = Cortical_CI

Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==1] = 'Default'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==2] = 'Visual'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==3] = 'Somatomotor'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==4] = 'Fronto-parietal'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==5] = 'Attention'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==6] = 'Cingulo-opercular'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==7] = 'Temporal'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==8] = 'Cingulo-parietal'
Cortical_df['Associated System'].loc[Cortical_df['Associated System'] ==9] = 'Sailency'
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

#within thalamus stuff
within_Tha_PC = pickle.load(open(path_to_data_folder+'within_Tha_PCs', "rb"))
within_Tha_WMD = pickle.load(open(path_to_data_folder+'within_Tha_WMDs', "rb"))

Thalamus_df['within_PC'] = within_Tha_PC/15
Thalamus_df['within_WMD'] = within_Tha_WMD/15

#correlation bewteen datasets
MGH_Tha_PC = pickle.load(open(path_to_data_folder+'MGH_Tha_PCs', "rb"))
NKI_Tha_PC = pickle.load(open(path_to_data_folder+'NKI_Tha_PCs', "rb"))

MGH_Tha_WMD = pickle.load(open(path_to_data_folder+'MGH_Tha_WMDs', "rb"))
NKI_Tha_WMD = pickle.load(open(path_to_data_folder+'NKI_Tha_WMDs', "rb"))

MGH_NNC = pickle.load(open(path_to_data_folder+'MGH_NNCs', "rb"))
NKI_NNC = pickle.load(open(path_to_data_folder+'NKI_NNCs', "rb"))

MGH_Tha_BNWR = pickle.load(open(path_to_data_folder+'MGH_Tha_BNWR', "rb"))
NKI_Tha_BNWR = pickle.load(open(path_to_data_folder+'NKI_Tha_BNWR', "rb"))

tsnr = np.loadtxt('/home/despoB/kaihwang/Rest/ROIs/tsnr')
Thalamus_df['TSNR'] = tsnr

Thalamus_df['MGH_PC'] = MGH_Tha_PC[320:]/13.5  #take out cortical ones, normalize by max PC
Thalamus_df['MGH_WMD'] = MGH_Tha_WMD[320:]/15 
Thalamus_df['MGH_NNC'] = MGH_NNC[320:]/15 
Thalamus_df['MGH_BNWR'] = MGH_Tha_BNWR[320:]/15

Thalamus_df['NKI_PC'] = NKI_Tha_PC[320:]/13.5  #take out cortical ones, normalize by max PC
Thalamus_df['NKI_WMD'] = NKI_Tha_WMD[320:]/15 
Thalamus_df['NKI_NNC'] = NKI_NNC[320:]/15 
Thalamus_df['NKI_BNWR'] = NKI_Tha_BNWR[320:]/15

################################################################
####### thalamo-cortical target's nodal role for each thalamus voxel
################################################################

Thalamus_target_df = pd.DataFrame(columns=('Voxel', 'Associated System', 'Cortical Target','Target PC','Target WMD' )) 

Cortical_targets = pickle.load(open(path_to_data_folder +'/Cortical_targets', "rb"))

for i, v in enumerate(Thalamus_voxels):

	tmp_df = pd.DataFrame(columns=('Voxel', 'Associated System', 'Cortical CI', 'Cortical Target','Target PC','Target WMD' ))
	tmp_PC = Cortical_PC[np.in1d(Cortical_ROIs, Cortical_targets[v])]
	tmp_WMD = Cortical_WMD[np.in1d(Cortical_ROIs, Cortical_targets[v])]
	tmp_targets = Cortical_ROIs[np.in1d(Cortical_ROIs, Cortical_targets[v])]
	tmp_CIs = Cortical_CI[np.in1d(Cortical_ROIs, Cortical_targets[v])]
	for u in range(len(tmp_PC)):
		tmp_df.set_value(u, 'Voxel', v)
		tmp_df.set_value(u, 'Associated System', Thalamus_CIs[i])
		tmp_df.set_value(u, 'Cortical CI', tmp_CIs[u])
		tmp_df.set_value(u, 'Cortical Target', tmp_targets[u])
		tmp_df.set_value(u, 'Target PC', tmp_PC[u]/13.5)
		tmp_df.set_value(u, 'Target WMD', tmp_WMD[u]/15)
	Thalamus_target_df = pd.concat([Thalamus_target_df, tmp_df], ignore_index=True)	

Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==1] = 'Default'
Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==2] = 'Visual'
Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==3] = 'Somatomotor'
Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==4] = 'Fronto-parietal'
Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==5] = 'Attention'
Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==6] = 'Cingulo-opercular'
Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==7] = 'Temporal'
Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==8] = 'Cingulo-parietal'
Thalamus_target_df['Associated System'].loc[Thalamus_target_df['Associated System'] ==11] = 'Sailency'



################################################################
###### save
################################################################
Thalamus_df.to_csv(path_to_data_folder+'Thalamus_nodal.csv')
Cortical_df.to_csv(path_to_data_folder+'Cortical_nodal.csv')
Thalamus_target_df.to_csv(path_to_data_folder+'Thalamus_target.csv')



