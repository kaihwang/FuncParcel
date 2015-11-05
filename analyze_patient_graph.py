#analyze patient data

from brain_graphs import *
from FuncParcel import *
import matplotlib.pylab as plt

################################################################
###### Setup
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
#Thalamocortical_corrmat = np.loadtxt(Parcel_path+'/MGH_Craddock_300_cortical_plus_thalamus_parcorrmatavg')

Cortical_CI = np.loadtxt(path_to_ROIs+'/Cortical_CI', dtype='int')
Thalamus_CIs = pickle.load(open(path_to_data_folder +'/MGH_Thalamus_CIs', "rb"))

Cortical_targets= pickle.load(open(path_to_data_folder +'/Cortical_targets', "rb"))
Cortical_nontargets= pickle.load(open(path_to_data_folder +'/Cortical_nontargets', "rb"))


path_to_lesion_masks = '/home/despoB/connectome-thalamus/Lesion_Masks/'
path_to_adjmat = '/home/despoB/connectome-thalamus/NotBackedUp/AdjMatrices/'

thalamic_patients = ['128', '162', '163', '168', '176']

Partition_CIs=np.array(range(1,12), dtype=int)
Network_names = ['Default', 'Visual', 'SM', 'FP', 'Attn', 'CO', 'Aud', 'CingP', 'OFC', 'IFT', 'Saliency']

################################################################
###### look at each subject's lesioned voxel distribution
################################################################

#load lesioned voxels indices
Lesioned_voxels = {}
for patient in thalamic_patients:
	Lesioned_voxels[patient] = np.loadtxt(path_to_lesion_masks+ patient+'_lesioned_voxels', dtype='int')


#look at ditribution of CI
Lesioned_CIs = {}
for patient in thalamic_patients:
	Lesioned_CIs[patient] = Thalamus_CIs[np.in1d(Thalamus_voxels, Lesioned_voxels[patient]).nonzero()[0]]


for patient in thalamic_patients:	
	To_Plot = np.zeros(Partition_CIs.size)
	for i, ci in enumerate(Partition_CIs):
		To_Plot[i] = sum(Lesioned_CIs[patient]==ci) 

	plt.figure()	
	plt.bar(Partition_CIs, To_Plot, align='center')
	plt.xticks(Partition_CIs, Network_names, rotation=30)
	plt.title(patient)
	plt.show()
	
#look at ditribution of nodal properties

Lesioned_values = np.array([])
for patient in thalamic_patients:
	Lesioned_values = np.concatenate((Lesioned_values, \
		Tha_WMDs_percentage[np.in1d(Thalamus_voxels, Lesioned_voxels[patient]).nonzero()[0]]))

plt.hist(Lesioned_values, 20)
plt.show()


################################################################
###### Load adj mats for all subjects, in 3d array ROIxROIxsubject
################################################################

AdjMat_Files = glob.glob(path_to_adjmat + 'MGH*Craddock_300_cortical_corrmat')

Control_AdjMats = np.loadtxt(AdjMat_Files[0])
for f in AdjMat_Files[1:]:

	M = np.loadtxt(f)
	Control_AdjMats = np.dstack((Control_AdjMats,M))

save_object(Control_AdjMats, path_to_data_folder +'/Control_AdjMats')
#np.dstack((a,a)).shape


################################################################
###### Create patient data dataframe
################################################################

Control_AdjMats = pickle.load(open('/home/despoB/connectome-thalamus/AvgMatrices/Control_AdjMats', "rb"))
Cortical_targets = pickle.load(open(path_to_data_folder +'/Cortical_targets', "rb"))
Cortical_nontargets = pickle.load(open(path_to_data_folder +'/Cortical_nontargets', "rb"))

Tha_PC = pickle.load(open(path_to_data_folder+'/Tha_PCs', "rb"))
Tha_WMD = pickle.load(open(path_to_data_folder+'/Tha_WMDs', "rb"))
Tha_BNWR = pickle.load(open(path_to_data_folder+'/Tha_BNWR', "rb"))

patient_df = pd.DataFrame()
for patient in thalamic_patients:
	Patient_adjmat = np.loadtxt(path_to_adjmat + 'Tha_' + patient + '_Craddock_300_cortical_corrmat' )

	tmp_df = pd.DataFrame()
	#tmp_df = pd.DataFrame(columns=('SubjID', 'Voxel', 'CI', 'PC', 'WMD', 'BNCR', \ #strange issue with indexing... 
	#	'Target_total_Weight', 'nonTarget_total_Weight', \
	#	'Target_total_Weight_bn', 'nonTarget_total_Weight_wn', \
	#	'Target_total_Weight_wn', 'nonTarget_total_Weight_wn', \
	#	'Target_connected_Weight', 'nonTarget_connected_Weight', \
	#	'Target_connected_Weight_bn', 'nonTarget_connected_Weight_bn', \
	#	'Target_connected_Weight_wn', 'nonTarget_connected_Weight_wn', ))

	for i, v in enumerate(Lesioned_voxels[patient]):
		tmp_df.set_value(i, 'SubjID', patient)
		tmp_df.set_value(i, 'Voxel', v)
		tmp_df.set_value(i, 'CI', Thalamus_CIs[Thalamus_voxels==v])
		tmp_df.set_value(i, 'PC', Tha_PC[Thalamus_voxels==v]/13.5)
		tmp_df.set_value(i, 'WMD', Tha_WMD[Thalamus_voxels==v]/15.)
		tmp_df.set_value(i, 'BNCR', Tha_BNWR[Thalamus_voxels==v]/15.)

		#logical position vecotros for targets and nontargets
		target_pos = np.in1d(Cortical_ROIs,Cortical_targets[v])
		nontarget_pos =  np.in1d(Cortical_ROIs,Cortical_nontargets[v])	

		#total weight
		tmp_df.set_value(i, 'Target_total_weight', \
			(np.nanmean(Patient_adjmat[target_pos,:]) - \
				np.nanmean(Control_AdjMats[target_pos,:,:])) / np.nanstd(Control_AdjMats[target_pos,:,:]))
		tmp_df.set_value(i, 'nonTarget_total_weight', \
			(np.nanmean(Patient_adjmat[nontarget_pos,:]) - \
				np.nanmean(Control_AdjMats[nontarget_pos,:,:])) / np.nanstd(Control_AdjMats[nontarget_pos,:,:]))

		#total weight bn
		tmp_df.set_value(i, 'Target_total_weight_bn', \
			(np.nanmean(Patient_adjmat[target_pos,:][:,Cortical_CI != Thalamus_CIs[Thalamus_voxels == v]]) - \
				np.nanmean(Control_AdjMats[target_pos,:,:][:,Cortical_CI != Thalamus_CIs[Thalamus_voxels == v],:])) \
			/ np.nanstd(Control_AdjMats[target_pos,:,:][:,Cortical_CI != Thalamus_CIs[Thalamus_voxels == v],:]))
		tmp_df.set_value(i, 'nonTarget_total_weight_bn', \
			(np.nanmean(Patient_adjmat[nontarget_pos,:][:,Cortical_CI != Thalamus_CIs[Thalamus_voxels == v]]) - \
				np.nanmean(Control_AdjMats[nontarget_pos,:,:][:,Cortical_CI != Thalamus_CIs[Thalamus_voxels == v],:])) \
			/ np.nanstd(Control_AdjMats[nontarget_pos,:,:][:,Cortical_CI != Thalamus_CIs[Thalamus_voxels == v],:]))

		#total weight wn
		tmp_df.set_value(i, 'Target_total_weight_wn', \
			(np.nanmean(Patient_adjmat[target_pos,:][:,Cortical_CI == Thalamus_CIs[Thalamus_voxels == v]]) - \
				np.nanmean(Control_AdjMats[target_pos,:,:][:,Cortical_CI == Thalamus_CIs[Thalamus_voxels == v],:])) \
			/ np.nanstd(Control_AdjMats[target_pos,:,:][:,Cortical_CI == Thalamus_CIs[Thalamus_voxels == v],:]))
		tmp_df.set_value(i, 'nonTarget_total_weight_wn', \
			(np.nanmean(Patient_adjmat[nontarget_pos,:][:,Cortical_CI == Thalamus_CIs[Thalamus_voxels == v]]) - \
				np.nanmean(Control_AdjMats[nontarget_pos,:,:][:,Cortical_CI == Thalamus_CIs[Thalamus_voxels == v],:])) \
			/ np.nanstd(Control_AdjMats[nontarget_pos,:,:][:,Cortical_CI == Thalamus_CIs[Thalamus_voxels == v],:]))
		
		#logical positional vectors for targets and non targets that are the same CI as the thalamic voxel
		target_wn_pos = np.in1d(Cortical_ROIs,Cortical_targets[v]) & (Cortical_CI == Thalamus_CIs[Thalamus_voxels == v])
		nontarget_wn_pos = np.in1d(Cortical_ROIs,Cortical_nontargets[v]) & (Cortical_CI == Thalamus_CIs[Thalamus_voxels == v])

		#logical positional vectors for targets and non targets that are not the same CI as the thalamic voxel
		target_bn_pos = np.in1d(Cortical_ROIs,Cortical_targets[v]) & (Cortical_CI != Thalamus_CIs[Thalamus_voxels == v])
		nontarget_bn_pos = np.in1d(Cortical_ROIs,Cortical_nontargets[v]) & (Cortical_CI != Thalamus_CIs[Thalamus_voxels == v])

		## extract the adj matrices using the positional vectors
		padj_t_tmp = Patient_adjmat[target_pos,:][:,target_pos]		
		padj_nt_tmp = Patient_adjmat[nontarget_pos,:][:,nontarget_pos]	

		cadj_t_tmp = Control_AdjMats[target_pos,:,:][:,target_pos,:]
		cadj_nt_tmp = Control_AdjMats[nontarget_pos,:,:][:,nontarget_pos,:]

		padj_t_wn_tmp = Patient_adjmat[target_wn_pos,:][:,target_wn_pos]		
		padj_nt_wn_tmp = Patient_adjmat[nontarget_wn_pos,:][:,nontarget_wn_pos]	

		cadj_t_wn_tmp = Control_AdjMats[target_wn_pos,:,:][:,target_wn_pos,:]
		cadj_nt_wn_tmp = Control_AdjMats[nontarget_wn_pos,:,:][:,nontarget_wn_pos,:]

		padj_t_bn_tmp = Patient_adjmat[target_bn_pos,:][:,target_bn_pos]		
		padj_nt_bn_tmp = Patient_adjmat[nontarget_bn_pos,:][:,nontarget_bn_pos]	

		cadj_t_bn_tmp = Control_AdjMats[target_bn_pos,:,:][:,target_bn_pos,:]
		cadj_nt_bn_tmp = Control_AdjMats[nontarget_bn_pos,:,:][:,nontarget_bn_pos,:]

		#target_connected_weight
		#only get upper triangle, get indices for that
		a = np.triu_indices(sum(target_pos),1)[0]
		b = np.triu_indices(sum(target_pos),1)[1]
		c = np.triu_indices(sum(nontarget_pos),1)[0]
		d = np.triu_indices(sum(nontarget_pos),1)[1]

		if any(a) and any(b):
			tmp_df.set_value(i, 'Target_connected_weight', (np.mean(padj_t_tmp[a,b]) - np.mean(cadj_t_tmp[a,b,:])) \
			/ np.std(cadj_t_tmp[a,b,:]))
		else:
			tmp_df.set_value(i, 'Target_connected_weight', 0)

		if any(c) and any(d):	
			tmp_df.set_value(i, 'nonTarget_connected_weight', (np.mean(padj_nt_tmp[c,d]) - np.mean(cadj_nt_tmp[c,d,:])) \
			/ np.std(cadj_nt_tmp[c,d,:]))
		else:
			tmp_df.set_value(i, 'nonTarget_connected_weight', 0)

		#target_connected_weight_bn
		a = np.triu_indices(sum(target_bn_pos),1)[0]
		b = np.triu_indices(sum(target_bn_pos),1)[1]
		c = np.triu_indices(sum(nontarget_bn_pos),1)[0]
		d = np.triu_indices(sum(nontarget_bn_pos),1)[1]

		if any(a) and any(b):
			tmp_df.set_value(i, 'Target_connected_weight_bn', (np.mean(padj_t_bn_tmp[a,b]) - np.mean(cadj_t_bn_tmp[a,b,:])) \
			/ np.std(cadj_t_bn_tmp[a,b,:]))
		else: 
			tmp_df.set_value(i, 'Target_connected_weight_bn', 0)

		if any(c) and any(d):	
			tmp_df.set_value(i, 'nonTarget_connected_weight_bn', (np.mean(padj_nt_bn_tmp[c,d]) - np.mean(cadj_nt_bn_tmp[c,d,:])) \
			/ np.std(cadj_nt_bn_tmp[c,d,:]))
		else:
			tmp_df.set_value(i, 'nonTarget_connected_weight_bn', 0)

		#target_connected_weight_wn
		a = np.triu_indices(sum(target_wn_pos),1)[0]
		b = np.triu_indices(sum(target_wn_pos),1)[1]
		c = np.triu_indices(sum(nontarget_wn_pos),1)[0]
		d = np.triu_indices(sum(nontarget_wn_pos),1)[1]

		if any(a) and any(b):
			tmp_df.set_value(i, 'Target_connected_weight_wn', (np.mean(padj_t_wn_tmp[a,b]) - np.mean(cadj_t_wn_tmp[a,b,:])) \
			/ np.std(cadj_t_wn_tmp[a,b,:]))
		else:
			tmp_df.set_value(i, 'Target_connected_weight_wn', 0)

		if any(c) and any(d):	
			tmp_df.set_value(i, 'nonTarget_connected_weight_wn', (np.mean(padj_nt_wn_tmp[c,d]) - np.mean(cadj_nt_wn_tmp[c,d,:])) \
			/ np.std(cadj_nt_wn_tmp[c,d,:]))
		else:
			tmp_df.set_value(i, 'nonTarget_connected_weight_wn', 0)

	patient_df = patient_df.append(tmp_df)		



###### save
patient_df['CI'].loc[patient_df['CI'] ==1] = 'Default'
patient_df['CI'].loc[patient_df['CI'] ==2] = 'Visual'
patient_df['CI'].loc[patient_df['CI'] ==3] = 'Somatomotor'
patient_df['CI'].loc[patient_df['CI'] ==4] = 'Fronto-parietal'
patient_df['CI'].loc[patient_df['CI'] ==5] = 'Attention'
patient_df['CI'].loc[patient_df['CI'] ==6] = 'Cingulo-opercular'
patient_df['CI'].loc[patient_df['CI'] ==7] = 'Temporal'
patient_df['CI'].loc[patient_df['CI'] ==8] = 'Cingulo-parietal'
patient_df['CI'].loc[patient_df['CI'] ==9] = 'Sailency'

patient_df.to_csv(path_to_data_folder + '/patient_df.csv', index = False)




################################################################
###### organize data not at node/roi level but at Network level
################################################################

Control_AdjMats = pickle.load(open('/home/despoB/connectome-thalamus/AvgMatrices/Control_AdjMats', "rb"))


patient_df = pd.DataFrame()
for patient in thalamic_patients:
	Patient_adjmat = np.loadtxt(path_to_adjmat + 'Tha_' + patient + '_Craddock_300_cortical_corrmat' )

	tmp_df = pd.DataFrame()
	for i, CI in enumerate(range(1,10)):

		tmp_df.set_value(i, 'SubjID', patient)
		tmp_df.set_value(i, 'CI', CI)

		p_bn = np.nanmean(Patient_adjmat[Cortical_CI==CI,:][:,Cortical_CI!=CI])
		c_bn_m = np.nanmean(Control_AdjMats[Cortical_CI==CI,:,:][:,Cortical_CI!=CI,:])
		c_bn_std = np.nanstd(Control_AdjMats[Cortical_CI==CI,:,:][:,Cortical_CI!=CI,:])

		p_wn = np.nanmean(Patient_adjmat[Cortical_CI==CI,:][:,Cortical_CI==CI])
		c_wn_m = np.nanmean(Control_AdjMats[Cortical_CI==CI,:,:][:,Cortical_CI==CI,:])
		c_wn_std = np.nanstd(Control_AdjMats[Cortical_CI==CI,:,:][:,Cortical_CI==CI,:])

		tmp_df.set_value(i, 'Between_network_connectivity_weight', p_bn)
		tmp_df.set_value(i, 'Within_network_connectivity_weight', p_wn)
	patient_df = patient_df.append(tmp_df)	


patient_df['CI'].loc[patient_df['CI'] ==1] = 'Default'
patient_df['CI'].loc[patient_df['CI'] ==2] = 'Visual'
patient_df['CI'].loc[patient_df['CI'] ==3] = 'Somatomotor'
patient_df['CI'].loc[patient_df['CI'] ==4] = 'Fronto-parietal'
patient_df['CI'].loc[patient_df['CI'] ==5] = 'Attention'
patient_df['CI'].loc[patient_df['CI'] ==6] = 'Cingulo-opercular'
patient_df['CI'].loc[patient_df['CI'] ==7] = 'Temporal'
patient_df['CI'].loc[patient_df['CI'] ==8] = 'Cingulo-parietal'
patient_df['CI'].loc[patient_df['CI'] ==9] = 'Sailency'

patient_df.to_csv(path_to_data_folder + '/patient_nework_df.csv', index = False)




################################################################
###### left v right
################################################################
Right_ROIs = np.unique(np.loadtxt(path_to_ROIs + '/Craddock_right_rois', dtype=int))
Left_ROIs = np.unique(np.loadtxt(path_to_ROIs + '/Craddock_left_rois', dtype=int))

Left_ROIs_pos = np.in1d(Cortical_ROIs, Left_ROIs)
Right_ROIs_pos = np.in1d(Cortical_ROIs, Right_ROIs)

patient_df = pd.DataFrame()
for patient in thalamic_patients:
	m = np.loadtxt(path_to_adjmat + 'Tha_' + patient + '_Craddock_300_cortical_corrmat' )
	Patient_adjmat = bct.threshold_proportional(m, 1)
	#Patient_adjmat[Patient_adjmat==0] = np.nan
	#Patient_adjmat_b = bct.weight_conversion(bct.threshold_proportional(m, 1),'binarize')

	tmp_df = pd.DataFrame()

	#initiate
	Left_p_wn = np.zeros(1) 
	Left_p_bn = np.zeros(1) 
	Right_p_wn = np.zeros(1) 
	Right_p_bn = np.zeros(1) 
	Left_p_tn = np.zeros(1) 
	Right_p_tn = np.zeros(1) 

	Left_p_ww = np.zeros(1) 
	Left_p_bw = np.zeros(1) 
	Right_p_ww = np.zeros(1) 
	Right_p_bw = np.zeros(1) 
	Left_p_tw = np.zeros(1) 
	Right_p_tw = np.zeros(1) 

	
	for CI in range(1,10):
		#total weight, sum
		Right_p_tw += np.nan_to_num(np.nansum(Patient_adjmat[(Cortical_CI==CI) & Right_ROIs_pos,:][:, Right_ROIs_pos]))
		Left_p_tw += np.nan_to_num(np.nansum(Patient_adjmat[(Cortical_CI==CI) & Left_ROIs_pos,:][:, Left_ROIs_pos]))

		#between weight, sum
		#Right_p_bw += np.nan_to_num(np.nansum(Patient_adjmat[(Cortical_CI==CI) & Right_ROIs_pos,:][:,(Cortical_CI!=CI) & Right_ROIs_pos]))
		#Left_p_bw += np.nan_to_num(np.nansum(Patient_adjmat[(Cortical_CI==CI) & Left_ROIs_pos,:][:,(Cortical_CI!=CI) & Left_ROIs_pos]))
		#within weight, sum
		

	for i, CI in enumerate(range(1,10)):

		#between weight, average
		Right_p_bn = np.nan_to_num(np.nanmean(Patient_adjmat[(Cortical_CI==CI) & Right_ROIs_pos,:][:,(Cortical_CI!=CI) & Right_ROIs_pos]))
		Left_p_bn = np.nan_to_num(np.nanmean(Patient_adjmat[(Cortical_CI==CI) & Left_ROIs_pos,:][:,(Cortical_CI!=CI) & Left_ROIs_pos]))
		#within weight, average
		Right_p_wn = np.nan_to_num(np.nanmean(Patient_adjmat[(Cortical_CI==CI) & Right_ROIs_pos,:][:,(Cortical_CI==CI) & Right_ROIs_pos]))
		Left_p_wn = np.nan_to_num(np.nanmean(Patient_adjmat[(Cortical_CI==CI) & Left_ROIs_pos,:][:,(Cortical_CI==CI) & Left_ROIs_pos]))
		
		#total weight, average
		Right_p_tn = np.nan_to_num(np.nanmean(Patient_adjmat[(Cortical_CI==CI) & Right_ROIs_pos,:][:, Right_ROIs_pos]))
		Left_p_tn = np.nan_to_num(np.nanmean(Patient_adjmat[(Cortical_CI==CI) & Left_ROIs_pos,:][:, Left_ROIs_pos]))

		Right_p_ww = np.nan_to_num(np.nansum(Patient_adjmat[(Cortical_CI==CI) & Right_ROIs_pos,:][:,(Cortical_CI==CI) & Right_ROIs_pos]))
		Left_p_ww = np.nan_to_num(np.nansum(Patient_adjmat[(Cortical_CI==CI) & Left_ROIs_pos,:][:,(Cortical_CI==CI) & Left_ROIs_pos]))


		#get modularity Q
		Right_sum_bw = np.zeros(1) 
		Left_sum_bw = np.zeros(1) 

		for CI2 in range(1,10):
			Right_tmp_bw = np.nan_to_num(np.nansum(Patient_adjmat[(Cortical_CI==CI) & Right_ROIs_pos,:][:,(Cortical_CI==CI2) & Right_ROIs_pos]))
			Right_sum_bw += (Right_tmp_bw / Right_p_tw)**2

			Left_tmp_bw = np.nan_to_num(np.nansum(Patient_adjmat[(Cortical_CI==CI) & Left_ROIs_pos,:][:,(Cortical_CI==CI2) & Left_ROIs_pos]))
			Left_sum_bw += (Left_tmp_bw / Left_p_tw)**2

		Right_Q = (Right_p_ww / Right_p_tw) - Right_sum_bw
		Left_Q = (Left_p_ww / Left_p_tw) - Left_sum_bw


		tmp_df.set_value(i, 'SubjID', patient)
		tmp_df.set_value(i, 'CI', CI)
		tmp_df.set_value(i, 'Left_Between_network_connectivity_weight', Left_p_bn)
		tmp_df.set_value(i, 'Right_Between_network_connectivity_weight', Right_p_bn)
		tmp_df.set_value(i, 'Left_Within_network_connectivity_weight', Left_p_wn)
		tmp_df.set_value(i, 'Right_Within_network_connectivity_weight', Right_p_wn)
		tmp_df.set_value(i, 'Left_Total_network_connectivity_weight', Left_p_tn)
		tmp_df.set_value(i, 'Right_Total_network_connectivity_weight', Right_p_tn)
		tmp_df.set_value(i, 'Left_Q', Left_Q)
		tmp_df.set_value(i, 'Right_Q', Right_Q)
	patient_df = patient_df.append(tmp_df)	

patient_df.to_csv(path_to_data_folder + '/patient_nework_df.csv', index = False)
################################################################
###### old stuff from HBM poster
################################################################

# from __future__ import division, print_function
# import numpy as np
# import scipy.io as sio
# import pickle
# import glob
# import pandas as pd
# import bct
# import os
# from surfer import Brain
# import FuncParcel
# from brainx import weighted_modularity
# import networkx as nx
# from sklearn.metrics import normalized_mutual_info_score
# #from ggplot import *

# # what to do?
# calulate_template_partition = False	
# identify_patient_cortical_targets = False
# calculate_z_scores = True
# visuazlie_template_partition = False
# visuazlie_patient_partition = False
# visualize_patient_cortical_target = False
# visualize_hubs = False
# run_template_partition_across_densities = False
# cal_sub_parition_by_densities = False
# cal_NMI = False

# # vector of cortical ROI index
# Cortical_ROIs = np.loadtxt('Data/Cortical_ROI_index')
# Cortical_CI = np.loadtxt('Data/Cortical_CI')

# # list of subjects
# #filename = '/home/despoB/connectome-thalamus/MGH/usable_sub'
# #Control_Subj = [line.rstrip('\n') for line in open(filename)]
# #Control_Subj.remove('Sub0094_Ses1')
# Control_Subj = ['1103', '1220', '1306', '1223', '1314', '1311', '1318', '1313', '1326', '1325', '1328', '1329', '1333', '1331', '1335', '1338', '1336', '1339', '1337', '1344']
# #Control_Subj = ['114', '116', '117', '118', '119', '201', '203', '204', '205', '206', '207', '208', '209', '210', '211', '212', '213', '214', '215', '216', '217', '219', '220']
# thalamic_patients = ['128', '162', '163', '168', '176']
# striatal_patients = ['b117', 'b122', 'b138', 'b143', 'b153']
# patients = thalamic_patients + striatal_patients
# Group = ['Control'] * len(Control_Subj) + ['Thalamic_patients'] * len(thalamic_patients) + ['Striatal_patietns'] * len(striatal_patients)

# # get template partion and nodal properties from MGH data
# if calulate_template_partition:

# 	# get partition
# 	AveMat = np.loadtxt('Data/CorticalAveMat')
# 	#W = bct.binarize(bct.threshold_proportional(AveMat, 0.05))
# 	graph = nx.from_numpy_matrix(bct.binarize(bct.threshold_proportional(AveMat, 0.075)))

# 	template_q = 0
# 	for i in xrange(0,20):
# 		print(i)
# 		louvain = weighted_modularity.LouvainCommunityDetection(graph)
# 		weighted_partitions = louvain.run()
# 		if weighted_partitions[0].modularity() > template_q:
# 			template_q = weighted_partitions[0].modularity()
# 			weighted_partition = weighted_partitions[0]
# 			template_ci = FuncParcel.convert_partition_dict_to_array(FuncParcel.convert_partition_to_dict(weighted_partition.communities), 297)

	
# 	#template_ci, template_q = bct.modularity_und(bct.binarize(bct.threshold_proportional(AveMat, 0.08))) #threshold at 0.05 cost
# 	#template_ci = np.loadtxt('Data/Cortical_CI') #use previously generated CI
# 	template_ci = template_ci.astype(int)
# 	np.savetxt('Data/MGH_CI', template_ci)
# 	Cortical_ROI_Coordinate = np.loadtxt('Data/Cortical_ROI_Coordinate')

# 	# get pc and wmd
# 	#template_pc = bct.participation_coef(bct.binarize(bct.threshold_proportional(AveMat, 0.08)), template_ci)
# 	#template_wmd = bct.module_degree_zscore(bct.binarize(bct.threshold_proportional(AveMat, 0.08)), template_ci)
# 	template_pc = FuncParcel.convert_graph_metric_dict_to_array(FuncParcel.participation_coefficient(weighted_partition), 297)
# 	template_wmd = FuncParcel.convert_graph_metric_dict_to_array(FuncParcel.within_community_degree(weighted_partition), 297)
# 	template_pc[np.isnan(template_pc)] = np.ma.masked
# 	template_wmd[np.isnan(template_wmd)] = np.ma.masked
# 	np.savetxt('Data/Cortical_PC', template_pc)
# 	np.savetxt('Data/Cortical_WMD', template_wmd)
# 	#outputdata
	
# 	template_nodal_data = pd.DataFrame()
# 	template_nodal_data['ROI'] = Cortical_ROIs
# 	template_nodal_data['PC'] = template_pc
# 	template_nodal_data['WMD'] = template_wmd
# 	template_nodal_data['Ci'] = template_ci
# 	#template_nodal_data['Coordinate'] = Cortical_ROI_Coordinate
# 	template_nodal_data.to_csv('Data/template_nodal_data.csv')

# 	#write out hubs
# 	connector_hubs =  template_nodal_data.ROI[np.argsort(template_pc)[::-1][0:30]].values
# 	connector_hubs = connector_hubs.astype(int)
# 	np.savetxt('Data/connector_hubs', connector_hubs, fmt='%3.d')

# 	provincial_hubs =  template_nodal_data.ROI[np.argsort(template_wmd)[::-1][0:30]].values
# 	provincial_hubs = provincial_hubs.astype(int)
# 	np.savetxt('Data/provincial_hubs', provincial_hubs, fmt='%3.d')

# 	both_hubs = np.intersect1d(provincial_hubs,connector_hubs)
# 	both_hubs = both_hubs.astype(int)
# 	np.savetxt('Data/both_hubs', both_hubs)

# if identify_patient_cortical_targets:


# 	#striatal patients
# 	#file_path = '/home/despoB/kaihwang/Rest/AdjMatrices/*Ses1_FIX_striatalcortical_corrmat'
# 	#AdjMat = FuncParcel.average_corrmat(file_path)
# 	#np.savetxt('/home/despoB/kaihwang/Rest/Striatum_parcel/StriatalCorticalAveMat', AdjMat)

# 	path_to_adjmat = '/home/despoB/kaihwang/Rest/Striatum_parcel/StriatalCorticalAveMat'
# 	path_to_list_of_subcorticalcortical_ROIs = '/home/despoB/kaihwang/bin/FuncParcel/Data/Striatalcortical_ROIs_index'
# 	path_to_list_of_subcortical_voxels = '/home/despoB/kaihwang/bin/FuncParcel/Data/striatal_voxel_index'
# 	path_to_list_of_cortical_ROIs = '/home/despoB/kaihwang/bin/FuncParcel/Data/Cortical_ROI_index'
	
# 	Subcorticalcortical_Targets = FuncParcel.parcel_subcortical_network(path_to_adjmat, path_to_list_of_subcorticalcortical_ROIs, path_to_list_of_subcortical_voxels, path_to_list_of_cortical_ROIs)
# 	Striatal_Patients_Cortical_Targets, Striatal_Patients_Cortical_NonTargets = FuncParcel.subcortical_patients_cortical_target(Subcorticalcortical_Targets, striatal_patients, 1)


# 	#thalamic patients
# 	#file_path = '/home/despoB/kaihwang/Rest/AdjMatrices/*Ses1_FIX_thalamocortical_corrmat'
# 	#AdjMat = FuncParcel.average_corrmat(file_path)
# 	#np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/ThalamoCorticalAveMat', AdjMat)

# 	path_to_adjmat = '/home/despoB/kaihwang/Rest/Thalamic_parcel/ThalamoCorticalAveMat'
# 	path_to_list_of_subcorticalcortical_ROIs = '/home/despoB/kaihwang/bin/FuncParcel/Data/Thalamocortical_ROIs_index'
# 	path_to_list_of_subcortical_voxels = '/home/despoB/kaihwang/bin/FuncParcel/Data/Thalamic_voxel_index'
# 	path_to_list_of_cortical_ROIs = '/home/despoB/kaihwang/bin/FuncParcel/Data/Cortical_ROI_index'
	
# 	Subcorticalcortical_Targets = FuncParcel.parcel_subcortical_network(path_to_adjmat, path_to_list_of_subcorticalcortical_ROIs, path_to_list_of_subcortical_voxels, path_to_list_of_cortical_ROIs)
# 	Thalamic_Patients_Cortical_Targets, Thalamic_Patients_Cortical_NonTargets = FuncParcel.subcortical_patients_cortical_target(Subcorticalcortical_Targets, thalamic_patients, 1)

# 	#combine dictionaries
# 	Patients_Cortical_Targets = Thalamic_Patients_Cortical_Targets.copy()
# 	Patients_Cortical_Targets.update(Striatal_Patients_Cortical_Targets)
# 	Patients_Cortical_NonTargets = Thalamic_Patients_Cortical_NonTargets.copy()
# 	Patients_Cortical_NonTargets.update(Striatal_Patients_Cortical_NonTargets)

# 	#save dict
# 	output = open('Data/Patients_Cortical_Targets.pkl', 'wb')
# 	pickle.dump(Patients_Cortical_Targets, output)
# 	output.close()

# 	output = open('Data/Patients_Cortical_NonTargets.pkl', 'wb')
# 	pickle.dump(Patients_Cortical_NonTargets, output)
# 	output.close()

# if calculate_z_scores:
# 	# create control's dataframe
# 	OlderControlGlobalData = pd.DataFrame()
# 	OlderControlNodalData = pd.DataFrame()
# 	# load old controls
# 	for s in Control_Subj:
# 		fn = '/home/despoB/kaihwang/Rest/Graph/gsetCI_%s.mat' %s
# 		GlobalData, NodalData = FuncParcel.convert_matlab_graph_str(fn, s, Cortical_ROIs)
# 		OlderControlGlobalData = OlderControlGlobalData.append(GlobalData)
# 		OlderControlNodalData = OlderControlNodalData.append(NodalData)
# 	OlderControlNodalData['Group'] = 'Control'
# 	OlderControlGlobalData['Group'] = 'Control'

# 	# load patients' graph metrics, convert global metrics into z_score 
# 	PatientsGlobalData = pd.DataFrame()
# 	PatientsNodalData = pd.DataFrame()
# 	for p in thalamic_patients:
# 		fn = '/home/despoB/kaihwang/Rest/Graph/gsetCI_%s.mat' %p
# 		GlobalData, NodalData = FuncParcel.convert_matlab_graph_str(fn, p, Cortical_ROIs)
# 		GlobalData['Q_zscore'] = FuncParcel.cal_graph_z_score(GlobalData, OlderControlGlobalData, np.array([]), 'Q')
# 		GlobalData['CC_zscore'] = FuncParcel.cal_graph_z_score(GlobalData, OlderControlGlobalData, np.array([]), 'CC')
# 		GlobalData['Group'] = 'Thalamic_patients'
# 		NodalData['Group'] = 'Thalamic_patients'
# 		PatientsGlobalData = PatientsGlobalData.append(GlobalData)
# 		PatientsNodalData = PatientsNodalData.append(NodalData)

# 	for p in striatal_patients:
# 		fn = '/home/despoB/kaihwang/Rest/Graph/gsetCI_%s.mat' %p
# 		GlobalData, NodalData = FuncParcel.convert_matlab_graph_str(fn, p, Cortical_ROIs)
# 		GlobalData['Q_zscore'] = FuncParcel.cal_graph_z_score(GlobalData, OlderControlGlobalData, np.array([]), 'Q')
# 		GlobalData['CC_zscore'] = FuncParcel.cal_graph_z_score(GlobalData, OlderControlGlobalData, np.array([]), 'CC')
# 		GlobalData['Group'] = 'Striatal_patients'
# 		#NodalData['PC_zscore'] = FuncParcel.cal_graph_z_score(NodalData, OlderControlNodalData, True,'PC')
# 		#NodalData['localE_zscore'] = FuncParcel.cal_graph_z_score(NodalData, OlderControlNodalData, True,'localE')
# 		NodalData['Group'] = 'Striatal_patients'
# 		PatientsGlobalData = PatientsGlobalData.append(GlobalData)
# 		PatientsNodalData = PatientsNodalData.append(NodalData)

	
# 	#organize nodal properties
# 	connector_hubs = np.loadtxt('Data/connector_hubs')
# 	connector_hubs = connector_hubs.astype(int)
# 	provincial_hubs = np.loadtxt('Data/provincial_hubs')
# 	provincial_hubs = provincial_hubs.astype(int)

# 	PatientsNodalData['node_type'] = 'none_hub'
# 	PatientsNodalData['node_type'].loc[PatientsNodalData['ROI'].isin(connector_hubs)] = 'connector_hub'
# 	PatientsNodalData['node_type'].loc[PatientsNodalData['ROI'].isin(provincial_hubs)] = 'provincial_hub'

# 	OlderControlNodalData['node_type'] = 'none_hub'
# 	OlderControlNodalData['node_type'].loc[OlderControlNodalData['ROI'].isin(connector_hubs)] = 'connector_hub'
# 	OlderControlNodalData['node_type'].loc[OlderControlNodalData['ROI'].isin(provincial_hubs)] = 'provincial_hub'
# 	#ControlDf = OlderControlNodalData.groupby(['Subject','ROI', 'Density']).aggregate(np.nanmean).reset_index()

# 	PatientsNodalData['target'] = False
# 	PatientsNodalData['non_target'] = False

# 	#load pickle saving target and nontarget
# 	pkl_file = open('Data/Patients_Cortical_Targets.pkl', 'rb')
# 	Patients_Cortical_Targets = pickle.load(pkl_file)
# 	pkl_file.close()

# 	pkl_file = open('Data/Patients_Cortical_NonTargets.pkl', 'rb')
# 	Patients_Cortical_NonTargets = pickle.load(pkl_file)
# 	pkl_file.close()

# 	for p in patients:
# 		#fn = 'Data/%s_cortical_target' %p 
# 		patient_target = Patients_Cortical_Targets[p]  #np.loadtxt(fn)
# 		#fn = 'Data/%s_cortical_nontarget' %p
# 		patient_nontarget = Patients_Cortical_NonTargets[p] #np.loadtxt(fn) 
# 		patient_target = patient_target.astype(int)
# 		patient_nontarget = patient_nontarget.astype(int)
# 		PatientsNodalData['target'].loc[(PatientsNodalData['Subject']==p) & (PatientsNodalData['ROI'].isin(patient_target))] = True
# 		PatientsNodalData['non_target'].loc[(PatientsNodalData['Subject']==p) & (PatientsNodalData['ROI'].isin(patient_nontarget))] = True

		
# 	#convert nodal graph metrics into z_score
# 	PatientsNodalZscoreData = pd.DataFrame()
# 	for p in patients:
# 		#fn = 'Data/%s_cortical_target' %p 
# 		patient_target = Patients_Cortical_Targets[p] 
# 		#fn = 'Data/%s_cortical_nontarget' %p
# 		patient_nontarget = Patients_Cortical_NonTargets[p]
# 		patient_target = patient_target.astype(int)
# 		patient_nontarget = patient_nontarget.astype(int)
# 		graph_metrics = ['Between_Module_Weight', 'Within_Module_Weight', 'Between_Module_Degree', 'Within_Module_Degree', 'localE', 'PC', 'WMD']
# 		patient_zDF= pd.DataFrame()

# 		#connector_hub 
# 		tmp_dict = {}
# 		target_selector = np.intersect1d(patient_target, connector_hubs)
# 		nontarget_selector = np.intersect1d(patient_nontarget, connector_hubs)
# 		if target_selector.any() & nontarget_selector.any():
# 			for metric in graph_metrics:
# 				tmp_dict['Target_' + metric] = FuncParcel.cal_graph_z_score(PatientsNodalData[PatientsNodalData['Subject']==p], OlderControlNodalData, target_selector , metric)
# 				tmp_dict['nonTarget_' + metric] = FuncParcel.cal_graph_z_score(PatientsNodalData[PatientsNodalData['Subject']==p], OlderControlNodalData, nontarget_selector, metric)
# 				tmp_dict['Density'] = np.linspace(0.01, 0.25, num = 49)
# 				tmp_dict['Subject'] = [p] * 49
# 				tmp_dict['node_type'] = ['connector_hub'] * 49

# 			tmpdf = pd.DataFrame(tmp_dict, columns=['Subject', 'Density', 'node_type', \
# 				'Target_Between_Module_Degree', 'Target_Within_Module_Degree', \
# 				'nonTarget_Between_Module_Degree', 'nonTarget_Within_Module_Degree', \
# 				'Target_Between_Module_Weight', 'Target_Within_Module_Weight', 'Target_PC', 'Target_WMD', 'Target_localE',\
# 				'nonTarget_Between_Module_Weight', 'nonTarget_Within_Module_Weight', 'nonTarget_PC', 'nonTarget_WMD', 'nonTarget_localE'])
# 			patient_zDF = patient_zDF.append(tmpdf)

# 		#provincial_hub
# 		tmp_dict = {}
# 		target_selector = np.intersect1d(patient_target, provincial_hubs)
# 		nontarget_selector = np.intersect1d(patient_nontarget, provincial_hubs)
# 		if target_selector.any() & nontarget_selector.any(): 
# 			for metric in graph_metrics:
# 				tmp_dict['Target_' + metric] = FuncParcel.cal_graph_z_score(PatientsNodalData[PatientsNodalData['Subject']==p], OlderControlNodalData, target_selector, metric)
# 				tmp_dict['nonTarget_' + metric] = FuncParcel.cal_graph_z_score(PatientsNodalData[PatientsNodalData['Subject']==p], OlderControlNodalData, nontarget_selector, metric)
# 				tmp_dict['Density'] = np.linspace(0.01, 0.25, num = 49)
# 				tmp_dict['Subject'] = [p] * 49
# 				tmp_dict['node_type'] = ['provincial_hub'] * 49

# 			tmpdf = pd.DataFrame(tmp_dict, columns=['Subject', 'Density', 'node_type', \
# 				'Target_Between_Module_Degree', 'Target_Within_Module_Degree', \
# 				'nonTarget_Between_Module_Degree', 'nonTarget_Within_Module_Degree', \
# 				'Target_Between_Module_Weight', 'Target_Within_Module_Weight', 'Target_PC', 'Target_WMD', 'Target_localE',\
# 				'nonTarget_Between_Module_Weight', 'nonTarget_Within_Module_Weight', 'nonTarget_PC', 'nonTarget_WMD', 'nonTarget_localE'])
# 			patient_zDF = patient_zDF.append(tmpdf)	
			
# 		#non_hub
# 		target_selector = np.array(list(set(patient_target) - set(np.intersect1d(patient_target, connector_hubs)) - set(np.intersect1d(patient_target, provincial_hubs))))
# 		nontarget_selector = np.array(list(set(patient_nontarget) - set(np.intersect1d(patient_nontarget, connector_hubs)) - set(np.intersect1d(patient_nontarget, provincial_hubs))))
# 		tmp_dict = {}
# 		if target_selector.any() & nontarget_selector.any():
# 			for metric in graph_metrics:
# 				tmp_dict['Target_' + metric] = FuncParcel.cal_graph_z_score(PatientsNodalData[PatientsNodalData['Subject']==p], OlderControlNodalData, target_selector , metric)
# 				tmp_dict['nonTarget_' + metric] = FuncParcel.cal_graph_z_score(PatientsNodalData[PatientsNodalData['Subject']==p], OlderControlNodalData, nontarget_selector, metric)
# 				tmp_dict['Density'] = np.linspace(0.01, 0.25, num = 49)
# 				tmp_dict['Subject'] = [p] * 49
# 				tmp_dict['node_type'] = ['non_hub'] * 49

# 			tmpdf = pd.DataFrame(tmp_dict, columns=['Subject', 'Density', 'node_type', \
# 				'Target_Between_Module_Degree', 'Target_Within_Module_Degree', \
# 				'nonTarget_Between_Module_Degree', 'nonTarget_Within_Module_Degree', \
# 				'Target_Between_Module_Weight', 'Target_Within_Module_Weight', 'Target_PC', 'Target_WMD', 'Target_localE',\
# 				'nonTarget_Between_Module_Weight', 'nonTarget_Within_Module_Weight', 'nonTarget_PC', 'nonTarget_WMD', 'nonTarget_localE'])
# 			patient_zDF = patient_zDF.append(tmpdf)
			
# 		#all
# 		tmp_dict = {}
# 		for metric in graph_metrics:
# 			tmp_dict['Target_' + metric] = FuncParcel.cal_graph_z_score(PatientsNodalData[PatientsNodalData['Subject']==p], OlderControlNodalData, patient_target , metric)
# 			tmp_dict['nonTarget_' + metric] = FuncParcel.cal_graph_z_score(PatientsNodalData[PatientsNodalData['Subject']==p], OlderControlNodalData, patient_nontarget , metric)
# 			tmp_dict['Density'] = np.linspace(0.01, 0.25, num = 49)
# 			tmp_dict['Subject'] = [p] * 49
# 			tmp_dict['node_type'] = ['all'] * 49
	
# 		tmpdf = pd.DataFrame(tmp_dict, columns=['Subject', 'Density', 'node_type', \
# 			'Target_Between_Module_Degree', 'Target_Within_Module_Degree', \
# 			'nonTarget_Between_Module_Degree', 'nonTarget_Within_Module_Degree', \
# 			'Target_Between_Module_Weight', 'Target_Within_Module_Weight', 'Target_PC', 'Target_WMD',  'Target_localE',\
# 			'nonTarget_Between_Module_Weight', 'nonTarget_Within_Module_Weight', 'nonTarget_PC', 'nonTarget_WMD',  'nonTarget_localE'])
# 		patient_zDF = patient_zDF.append(tmpdf)
		
# 		PatientsNodalZscoreData = PatientsNodalZscoreData.append(patient_zDF)

# 	PatientsNodalZscoreData['Group']='na'
# 	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='128'] = 'Thalamic_patients'
# 	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='162'] = 'Thalamic_patients'
# 	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='163'] = 'Thalamic_patients'
# 	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='168'] = 'Thalamic_patients'
# 	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='176'] = 'Thalamic_patients'
# 	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='b117'] = 'Striatal_patients'
# 	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='b122'] = 'Striatal_patients'
# 	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='b138'] = 'Striatal_patients'
# 	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='b143'] = 'Striatal_patients'
# 	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='b153'] = 'Striatal_patients'
# 	PatientsNodalZscoreData.to_csv('Data/PatientsNodalZscoreData.csv')		

# 	#combine dataframes
# 	GraphGlobalData = pd.DataFrame()
# 	GraphNodalData = pd.DataFrame()

# 	GraphGlobalData = PatientsGlobalData.append(OlderControlGlobalData)
# 	GraphNodalData = PatientsNodalData.append(OlderControlNodalData)

# 	GraphGlobalData.to_csv('Data/GraphGlobalData.csv')
# 	GraphNodalData.to_csv('Data/GraphNodalData.csv')

# 	# calculate hemispheric difference
# 	GraphGlobalData['RightvLeft_Q'] = GraphGlobalData['right_Q'] - GraphGlobalData['left_Q']
# 	GraphGlobalData['LeftvRight_Q'] = GraphGlobalData['left_Q'] - GraphGlobalData['right_Q']
# 	GraphGlobalData.to_csv('Data/GraphGlobalData.csv')



# #run partiion across threshold on the MGH's avemat across a range of densities
# if run_template_partition_across_densities:
# 	GraphNodalData = pd.DataFrame.from_csv('Data/GraphNodalData.csv')
# 	GraphNodalData['Subject'] = GraphNodalData['Subject'].astype('string').values
# 	AveMat = np.loadtxt('Data/CorticalAveMat')
# 	MGH_template_partition = pd.DataFrame()
# 	row_count = 0
# 	for d in np.unique(GraphNodalData.Density): #extracting density vectors from previosly created density values to avoid floating error....
# 		ci, q = bct.modularity_louvain_und(bct.binarize(bct.threshold_proportional(AveMat, d)))
# 		for roi in np.arange(0, len(ci)):
# 			MGH_template_partition.loc[row_count, 'Ci'] = ci[roi]
# 			MGH_template_partition.loc[row_count, 'ROI'] = Cortical_ROIs[roi]
# 			MGH_template_partition.loc[row_count, 'Density'] = d
# 			row_count = row_count +1
# 	MGH_template_partition.to_csv('Data/MGH_partition.csv')


# # #run parition on subject data
# if cal_sub_parition_by_densities:
# 	Subject_partition = pd.DataFrame()
# 	row_count = 0
# 	GraphNodalData = pd.DataFrame.from_csv('Data/GraphNodalData.csv')
# 	GraphNodalData['Subject'] = GraphNodalData['Subject'].astype('string').values
# 	for d in np.unique(GraphNodalData.Density): #extracting density vectors from previosly created density values to avoid floating error....
# 		for s in patients+Control_Subj:
# 			fn = '/home/despoB/kaihwang/Rest/AdjMatrices/t%s_Full_WashU333_corrmat' %s
# 			Mat = np.loadtxt(fn) 
# 			ci, q = bct.modularity_louvain_und(bct.binarize(bct.threshold_proportional(Mat, d)))
# 			for roi in np.arange(0, len(ci)):
# 				Subject_partition.loc[row_count, 'Ci'] = ci[roi]
# 				Subject_partition.loc[row_count, 'ROI'] = Cortical_ROIs[roi]
# 				Subject_partition.loc[row_count, 'Density'] = d
# 				Subject_partition.loc[row_count, 'Subject'] = s
# 				Subject_partition.loc[row_count, 'Group'] = GraphNodalData[GraphNodalData['Subject']==s]['Group'].values[0]
# 				row_count = row_count +1
# 	Subject_partition.to_csv('Data/Subject_partition.csv')

# # # calculate mutual information
# if cal_NMI:

# 	NMI_dataframe = pd.DataFrame()
# 	row_count = 0
# 	for s in Control_Subj + patients:
# 		fn = '/home/despoB/kaihwang/bin/FuncParcel/Data/Subject_Partition/%s_ci' %s
# 		subject_ci = np.genfromtxt(fn)
# 		template_ci = np.loadtxt('Data/MGH_CI') #MGH_template_partition['Ci'].values
# 		subject_ci = subject_ci.astype(int)
# 		template_ci = template_ci.astype(int)

# 		#take out single partitions
# 		for i in np.unique(subject_ci):
# 			if np.count_nonzero(subject_ci==i) ==1:
# 				subject_ci[subject_ci==i] = np.ma.masked


# 		NMI_dataframe.loc[row_count,'NMI'] = normalized_mutual_info_score(subject_ci, template_ci)
# 		NMI_dataframe.loc[row_count,'Subject'] = s
# 		NMI_dataframe.loc[row_count,'Group'] = Group[row_count]
# 		row_count = row_count+1
# 	NMI_dataframe.to_csv('Data/NMI_dataframe.csv')


# #try to visulize template graph partition
# if visuazlie_template_partition:
# 	template_ci = np.loadtxt('Data/MGH_CI')
# 	Cortical_ROI_Coordinate = np.loadtxt('Data/Cortical_ROI_Coordinate')
# 	subjects_dir = os.environ["SUBJECTS_DIR"]
# 	subject_id, surface = 'fsaverage', 'inflated'
# 	hemi = 'rh'
# 	brain = Brain(subject_id, hemi, surface, views=['med'], config_opts=dict(background="white"))

# 	#bmap = brewer2mpl.get_map('Paired', 'Qualitative', 12)
# 	colors = ['blue', 'cyan', 'purple', 'yellow', '#ffc966', 'grey', 'black', 'red']
# 	c_i = 0
# 	for i in xrange(0, int(np.max(template_ci)+1)):
# 		coor = Cortical_ROI_Coordinate[template_ci==i]
# 		print(coor)
# 		if len(coor)>3:
# 			#brain.add_foci(coor[coor[:,0]<0], map_surface="white", color=colors[c_i], hemi="lh" )
# 			brain.add_foci(coor[coor[:,0]>0], map_surface="white", color=colors[c_i], hemi="rh" )
# 			c_i = c_i+1
# 	brain.save_image('Data/template_parition_med.png')
# 	#brain.close()	

# #try to visulize patients' graph partition
# if visuazlie_patient_partition:
# 	Cortical_ROI_Coordinate = np.loadtxt('Data/Cortical_ROI_Coordinate')
# 	subject_id, surface = 'fsaverage', 'inflated'
# 	hemi = 'rh'
# 	subjects_dir = os.environ["SUBJECTS_DIR"]
# 	colors = ['#00ffff', '#000000', '#0000ff', '#ff00ff', '#008000', '#808080', '#00ff00', '#800000', '#000080', '#808000', '#800080', '#ff0000', '#c0c0c0', '#008080', '#ffffff', '#ffff00']
# 	for s in ['176']:
# 		fn = 'Data/Subject_Partition/%s_ci' %s
# 		patient_ci = np.loadtxt(fn)
# 		brain = Brain(subject_id, hemi, surface, views=['lat'], config_opts=dict(background="white"))
# 		c_i = 0
# 		for i in xrange(0, int(np.max(patient_ci)+1)):
# 			coor = Cortical_ROI_Coordinate[patient_ci==i]
# 			print(coor)
# 			if len(coor)>3:
# 				#brain.add_foci(coor[coor[:,0]<0], map_surface="white", color=colors[c_i], hemi="lh" )
# 				brain.add_foci(coor[coor[:,0]>0], map_surface="white", color=colors[c_i], hemi="rh" )
# 				c_i = c_i+1
# 		fn = 'Data/Subject_Partition/%s_lat_ci.png' %s		
# 		brain.save_image(fn)			
# 		#brain.close()		
			


# # visualize hubs
# if visualize_hubs:
# 	subjects_dir = os.environ["SUBJECTS_DIR"]
# 	subject_id, surface = 'fsaverage', 'inflated'
# 	hemi = 'split'
# 	Cortical_ROI_Coordinate = np.loadtxt('Data/Cortical_ROI_Coordinate')
# 	connector_hubs = np.loadtxt('Data/connector_hubs')
# 	provincial_hubs = np.loadtxt('Data/provincial_hubs')
# 	both_hubs = np.loadtxt('Data/both_hubs')
# 	brain2 = Brain(subject_id, hemi, surface, views=['lat', 'med'], config_opts=dict(background="white"))
# 	coor = Cortical_ROI_Coordinate[np.in1d(Cortical_ROIs,connector_hubs)]
# 	brain2.add_foci(coor[coor[:,0]<0], scale_factor=2,  map_surface="white", color='red', hemi="lh" )
# 	brain2.add_foci(coor[coor[:,0]>0], scale_factor=2,map_surface="white", color='red', hemi="rh" )
# 	coor = Cortical_ROI_Coordinate[np.in1d(Cortical_ROIs,provincial_hubs)]
# 	brain2.add_foci(coor[coor[:,0]<0], scale_factor=2,map_surface="white", color='blue', hemi="lh" )
# 	brain2.add_foci(coor[coor[:,0]>0], scale_factor=2,map_surface="white", color='blue', hemi="rh" )
# 	coor = Cortical_ROI_Coordinate[np.in1d(Cortical_ROIs,both_hubs)]
# 	brain2.add_foci(coor[coor[:,0]<0], scale_factor=2,map_surface="white", color='black', hemi="lh" )
# 	brain2.add_foci(coor[coor[:,0]>0], scale_factor=2,map_surface="white", color='black', hemi="rh" )
# 	brain2.save_image('hubs.png')	

# #visualize patients cortical target
# if visualize_patient_cortical_target:
# 	subjects_dir = os.environ["SUBJECTS_DIR"]
# 	subject_id, surface = 'fsaverage', 'inflated'
# 	hemi = 'split'
# 	patient_target = np.loadtxt('Data/128_cortical_target')
# 	patient_target = patient_target.astype(int)
# 	brain3 = Brain(subject_id, hemi, surface,
# 	              config_opts=dict(width=800, height=400, background="white"))
# 	coor = Cortical_ROI_Coordinate[template_nodal_data['ROI'].isin(patient_target).values]
# 	brain3.add_foci(coor[coor[:,0]<0], map_surface="white", color='red', hemi="lh" )
# 	brain3.add_foci(coor[coor[:,0]>0], map_surface="white", color='red', hemi="rh" )



# # plot patient nodal data
# #plotData = pd.DataFrame.from_csv('PatientsNodalZscoreData.csv')
# #plotData = GraphNodalData.loc[GraphNodalData['target']==True].groupby(['Subject', 'Density','node_type']).aggregate(np.nanmean).reset_index()
# #ggplot(aes(x = 'Density', y = 'Target_localE', color = 'Subject'),data = plotData) \
# #+ geom_line() + xlim(0.05, 0.15) + facet_wrap('node_type', ncol = 4)

# #cleanup
# #reset_selective GlobalData
# #reset_selective NodalData
# #reset_selective OlderControlGlobalData
# #reset_selective OlderControlNodalData


