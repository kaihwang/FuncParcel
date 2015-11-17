#analyze patient data

from brain_graphs import *
from FuncParcel import *
import matplotlib.pylab as plt
from scipy.stats.mstats import zscore as zscore

################################################################
###### Setup
################################################################

path_to_ROIs = '/home/despoB/connectome-thalamus/ROIs'
path_to_data_folder = '/home/despoB/kaihwang/bin/FuncParcel/Data'
Parcel_path = '/home/despoB/connectome-thalamus/Thalamic_parcel'

Cortical_ROIs = np.loadtxt(path_to_ROIs+'/Gordon_333', dtype = int)
Thalamus_voxels = np.loadtxt(path_to_ROIs+'/thalamus_voxel_indices', dtype = int)
ROIs = np.append(Cortical_ROIs, Thalamus_voxels)

Cortical_CI = np.loadtxt(path_to_ROIs + '/Gordon_consensus_CI')
Thalamus_CIs = np.loadtxt(Parcel_path + '/Thalamus_clusters_corticalCI')

Cortical_plus_thalamus_CI = np.append(Cortical_CI, Thalamus_CIs)
Cortical_ROIs_positions = np.arange(0,len(Cortical_CI),1)
Thalamus_voxel_positions = np.arange(len(Cortical_CI),len(Cortical_plus_thalamus_CI),1)


#Cortical_targets= pickle.load(open(path_to_data_folder +'/Cortical_targets', "rb"))
#Cortical_nontargets= pickle.load(open(path_to_data_folder +'/Cortical_nontargets', "rb"))

path_to_lesion_masks = '/home/despoB/connectome-thalamus/Lesion_Masks/'
path_to_adjmat = '/home/despoB/connectome-thalamus/NotBackedUp/AdjMatrices/'

thalamic_patients = ['128', '162', '163', '168', '176']

Partition_CIs = np.unique(Cortical_CI[Cortical_CI!=0])
Network_names = ['DF', 'CO', 'SM', 'FP', 'OP', 'V', 'RT', 'ST', 'Attn', 'T']

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
	cd 

################################################################
###### Load adj mats for all subjects, in 3d array ROIxROIxsubject
################################################################

AdjMat_Files = glob.glob(path_to_adjmat + 'MGH*Gordon_333_cortical_corrmat')
Control_AdjMats = np.loadtxt(AdjMat_Files[0])
#Control_AdjMats = zscore(Control_AdjMats, axis = None)
for f in AdjMat_Files[1:]:

	M = np.loadtxt(f)
	#M = zscore(M, axis = None)
	Control_AdjMats = np.dstack((Control_AdjMats,M))

save_object(Control_AdjMats, path_to_data_folder +'/Control_AdjMats_Gordon')


AdjMat_Files = glob.glob(path_to_adjmat + 'Tha*Gordon_333_cortical_corrmat')
Patient_AdjMats = np.loadtxt(AdjMat_Files[0])
Patient_AdjMats = np.nan_to_num(Patient_AdjMats)
#Patient_AdjMats = zscore(Patient_AdjMats, axis = None)
for f in AdjMat_Files[1:]:

	M = np.loadtxt(f)
	M = np.nan_to_num(M)
	#M = zscore(M, axis = None)
	Patient_AdjMats = np.dstack((Patient_AdjMats,M))

save_object(Control_AdjMats, path_to_data_folder +'/Patient_AdjMats_Gordon')
#np.dstack((a,a)).shape

################################################################
###### Look at changes in overal network modular structure
################################################################

control_df = pd.DataFrame()
for s in range(0, np.shape(Control_AdjMats)[2]):
	right_M = Control_AdjMats[Gordon_right_ROI_positions,:,:][:,Gordon_right_ROI_positions,:][:,:,s].copy()
	left_M = Control_AdjMats[Gordon_left_ROI_positions,:,:][:,Gordon_left_ROI_positions,:][:,:,s].copy()
	right_M = np.nan_to_num(right_M)
	left_M = np.nan_to_num(left_M)
	M = np.nan_to_num(Control_AdjMats[:,:,s].copy())

	tmp_df = pd.DataFrame()
	for i, p in enumerate(np.arange(0.02, 0.16, 0.01)):
		right_Q = cal_modularity_w_imposed_community(bct.threshold_proportional(right_M, p), Cortical_CI[Gordon_right_ROI_positions])
		left_Q = cal_modularity_w_imposed_community(bct.threshold_proportional(left_M, p), Cortical_CI[Gordon_left_ROI_positions])

		tmp_df.set_value(i, 'SubjID', s)
		tmp_df.set_value(i, 'Density', p)
		tmp_df.set_value(i, 'left_Q', left_Q)
		tmp_df.set_value(i, 'right_Q', right_Q)
		tmp_df.set_value(i, 'DF_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 1))
		tmp_df.set_value(i, 'CO_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 2))
		tmp_df.set_value(i, 'SM_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 3))
		tmp_df.set_value(i, 'FP_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 4))
		tmp_df.set_value(i, 'OP_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 5))
		tmp_df.set_value(i, 'V_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 6))
		tmp_df.set_value(i, 'RS_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 7))
		tmp_df.set_value(i, 'ST_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 8))
		tmp_df.set_value(i, 'ATTN_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 9))
		tmp_df.set_value(i, 'T_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 12))

	control_df = control_df.append(tmp_df)

Patient_df = pd.DataFrame()
for s in range(0, np.shape(Patient_AdjMats)[2]):
	right_M = Patient_AdjMats[Gordon_right_ROI_positions,:,:][:,Gordon_right_ROI_positions,:][:,:,s].copy()
	left_M = Patient_AdjMats[Gordon_left_ROI_positions,:,:][:,Gordon_left_ROI_positions,:][:,:,s].copy()
	right_M = np.nan_to_num(right_M)
	left_M = np.nan_to_num(left_M)
	M = np.nan_to_num(Patient_AdjMats[:,:,s].copy())

	tmp_df = pd.DataFrame()
	for i, p in enumerate(np.arange(0.02, 0.16, 0.01)):
		right_Q = cal_modularity_w_imposed_community(bct.threshold_proportional(right_M, p), Cortical_CI[Gordon_right_ROI_positions])
		left_Q = cal_modularity_w_imposed_community(bct.threshold_proportional(left_M, p), Cortical_CI[Gordon_left_ROI_positions])

		tmp_df.set_value(i, 'SubjID', s)
		tmp_df.set_value(i, 'Density', p)
		tmp_df.set_value(i, 'left_Q', left_Q)
		tmp_df.set_value(i, 'right_Q', right_Q)
		tmp_df.set_value(i, 'DF_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 1))
		tmp_df.set_value(i, 'CO_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 2))
		tmp_df.set_value(i, 'SM_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 3))
		tmp_df.set_value(i, 'FP_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 4))
		tmp_df.set_value(i, 'OP_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 5))
		tmp_df.set_value(i, 'V_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 6))
		tmp_df.set_value(i, 'RS_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 7))
		tmp_df.set_value(i, 'ST_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 8))
		tmp_df.set_value(i, 'ATTN_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 9))
		tmp_df.set_value(i, 'T_Q', cal_modularity_community(bct.threshold_proportional(M, p), Cortical_CI, 12))

	Patient_df = Patient_df.append(tmp_df)

Patient_df['L-R_Q'] = Patient_df['left_Q'] - Patient_df['right_Q']
Patient_df['R-L_Q'] = Patient_df['right_Q'] - Patient_df['left_Q']

control_df['R-L_Q'] = control_df['right_Q'] - control_df['left_Q']
control_df['L-R_Q'] = control_df['left_Q'] - control_df['right_Q']



def cal_modularity_w_imposed_community(M, CI):
	Total_weight = M.sum()
	Q=0.0
	for i in np.unique(CI):
		Within_weight = np.sum(M[CI==i,:][:,CI==i])
		Within_weight_ratio = Within_weight / Total_weight
		Between_weight = 0.0
		for j in np.unique(CI):
			if i !=j:
				Between_weight += (np.sum(M[CI==i,:][:,CI==j]) / Total_weight)
		Between_weight_ratio = (Between_weight)**2
		Q += (Within_weight_ratio - Between_weight_ratio)
	return Q

def cal_modularity_community(M, CI, targetCI):
	Total_weight = M.sum()
	Q=0.0
	for i in np.unique(targetCI):
		Within_weight = np.sum(M[CI==i,:][:,CI==i])
		Within_weight_ratio = Within_weight / Total_weight
		Between_weight = 0.0
		for j in np.unique(CI):
			if i !=j:
				Between_weight += (np.sum(M[CI==i,:][:,CI==j]) / Total_weight)
		Between_weight_ratio = (Between_weight)**2
		Q += (Within_weight_ratio - Between_weight_ratio)
	return Q

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








