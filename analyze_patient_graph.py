#analyze patient data
from __future__ import division, print_function
#from brain_graphs import *
from FuncParcel import *
import matplotlib.pylab as plt
from scipy.stats.mstats import zscore as zscore
from brain_graphs import *

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
Thalamus_CIs = np.loadtxt(Parcel_path + '/MGH_Thalamus_WTA_CIs')
Thalamus_Morel = np.loadtxt(Parcel_path + '/Morel_parcel') #morel atlas
Thalamus_FSL = np.loadtxt(Parcel_path + '/fsl_thalamus_ana_parcel') #fsl atls
Morel_mask = np.loadtxt('/home/despoB/connectome-thalamus/Thalamic_parcel/morel_mask')
mask_value = Morel_mask==0
Thalamus_CIs[mask_value] = 0
Thalamus_Morel[mask_value] = 0
Thalamus_FSL[mask_value] = 0

Cortical_plus_thalamus_CI = np.append(Cortical_CI, Thalamus_CIs)
Cortical_ROIs_positions = np.arange(0,len(Cortical_CI),1)
Thalamus_voxel_positions = np.arange(len(Cortical_CI),len(Cortical_plus_thalamus_CI),1)


#Cortical_targets= pickle.load(open(path_to_data_folder +'/Cortical_targets', "rb"))
#Cortical_nontargets= pickle.load(open(path_to_data_folder +'/Cortical_nontargets', "rb"))

path_to_lesion_masks = '/home/despoB/connectome-thalamus/Lesion_Masks/'
path_to_adjmat = '/home/despoB/connectome-thalamus/NotBackedUp/AdjMatrices/'

thalamic_patients = ['176', '128', '168', '163']  #S1-4

Partition_CIs = np.unique(Cortical_CI[Cortical_CI!=0])
Network_names = ['DF', 'CO', 'SM', 'FP', 'latO', 'mO', 'mT', 'T', 'sFP', 'T']

################################################################
###### look at each subject's lesioned voxel distribution
################################################################

#load lesioned voxels indices
Lesioned_voxels = {}
for patient in thalamic_patients:
	Lesioned_voxels[patient] = np.loadtxt(path_to_lesion_masks+ patient+'_lesioned_voxels', dtype='int')

PC = pickle.load(open('/home/despoB/kaihwang/Rest/Graph/MGH_avemat_tha_nodal_pcorr_PCs','rb')) 

Lesioned_PC = {}
for patient in thalamic_patients:
	Lesioned_PC[patient] = np.sum(PC[:,333:][:,np.in1d(Thalamus_voxels, Lesioned_voxels[patient])],axis=1)[::-1]

Lesioned_CIs = {}
Lesioned_Morel = {}
Lesioned_FSL = {}
for patient in thalamic_patients:
	Lesioned_CIs[patient] = Thalamus_CIs[np.in1d(Thalamus_voxels, Lesioned_voxels[patient])]
	Lesioned_Morel[patient] = Thalamus_Morel[np.in1d(Thalamus_voxels, Lesioned_voxels[patient])]
	Lesioned_FSL[patient] = Thalamus_FSL[np.in1d(Thalamus_voxels, Lesioned_voxels[patient])]



################################################################
###### Organize adj mats for all subjects, in 3d array ROIxROIxsubject
################################################################

def pool_control_adjmats():
	AdjMat_Files = glob.glob(path_to_adjmat + 'MGH*Gordon_333_cortical_corrmat')
	Control_AdjMats = np.loadtxt(AdjMat_Files[0])
	#Control_AdjMats = zscore(Control_AdjMats, axis = None)
	for f in AdjMat_Files[1:]:

		M = np.loadtxt(f)
		#M = zscore(M, axis = None)
		Control_AdjMats = np.dstack((Control_AdjMats,M))
	save_object(Control_AdjMats, path_to_data_folder +'/Control_AdjMats_Gordon')
	return Control_AdjMats

def pool_patient_adjmats():
	M_S1 = M = np.loadtxt('/home/despoB/kaihwang/Rest/NotBackedUp/ParMatrices/Tha_176_Gordon_333_cortical_corrmat')
	M_S2 = M = np.loadtxt('/home/despoB/kaihwang/Rest/NotBackedUp/ParMatrices/Tha_128_Gordon_333_cortical_corrmat')
	M_S3 = M = np.loadtxt('/home/despoB/kaihwang/Rest/NotBackedUp/ParMatrices/Tha_168_Gordon_333_cortical_corrmat')
	M_S4 = M = np.loadtxt('/home/despoB/kaihwang/Rest/NotBackedUp/ParMatrices/Tha_163_Gordon_333_cortical_corrmat')
	Patient_AdjMats = np.dstack((M_S1, M_S2, M_S3, M_S4))
	return Patient_AdjMats

#Control_AdjMats = save_control_adjmats()
#Patient_AdjMats = pool_patient_adjmats()

################################################################
###### Analyze patient graph, Look at changes in overal network modular structure
################################################################
Gordon_right_ROI_positions = np.loadtxt(path_to_ROIs + '/Gordon_right_ROIs_positions')
Gordon_right_ROI_positions = Gordon_right_ROI_positions.astype(bool)
Gordon_left_ROI_positions = np.loadtxt(path_to_ROIs + '/Gordon_left_ROIs_positions')
Gordon_left_ROI_positions  = Gordon_left_ROI_positions.astype(bool)

def cal_modularity_w_indiv_infomap(AdjMats):
	'''calculate Q separately for whole brain, left/right hemispheres, using subject-level infomap parition.
	'''
	output_df = pd.DataFrame()
	for s in range(0, np.shape(AdjMats)[2]):
		#load matrices
		right_M = AdjMats[Gordon_right_ROI_positions,:,:][:,Gordon_right_ROI_positions,:][:,:,s].copy()
		left_M = AdjMats[Gordon_left_ROI_positions,:,:][:,Gordon_left_ROI_positions,:][:,:,s].copy()
		right_M = np.nan_to_num(right_M)
		left_M = np.nan_to_num(left_M)
		M = np.nan_to_num(AdjMats[:,:,s].copy())

		#partitions
		g= recursive_network_partition(matrix=M.copy(), min_cost=.01,max_cost=0.15, min_community_size=5 ,min_weight=2)
		CI = np.array(g.community.membership)

		r_g= recursive_network_partition(matrix=right_M.copy(), min_cost=.01,max_cost=0.15, min_community_size=5 ,min_weight=2)
		r_CI = np.array(r_g.community.membership)

		l_g= recursive_network_partition(matrix=left_M.copy(), min_cost=.01,max_cost=0.15, min_community_size=5 ,min_weight=2)
		l_CI = np.array(l_g.community.membership)
		
		tmp_df = pd.DataFrame()
		for i, p in enumerate(np.arange(0.01, 0.16, 0.01)):
			#cal Q using individual level whole-brain partition
			Q_indwbci = cal_modularity_w_imposed_community(bct.threshold_proportional(M, p, copy=True), CI)
			right_Q_indwbci = cal_modularity_w_imposed_community(bct.threshold_proportional(right_M, p, copy=True), CI[Gordon_right_ROI_positions])
			left_Q_indwbci = cal_modularity_w_imposed_community(bct.threshold_proportional(left_M, p, copy=True), CI[Gordon_left_ROI_positions])

			#cal Q using group level partition
			Q_groupwbci = cal_modularity_w_imposed_community(bct.threshold_proportional(M, p, copy=True), Cortical_CI)
			right_Q_groupwbci = cal_modularity_w_imposed_community(bct.threshold_proportional(right_M, p, copy=True), Cortical_CI[Gordon_right_ROI_positions])
			left_Q_groupwbci = cal_modularity_w_imposed_community(bct.threshold_proportional(left_M, p, copy=True), Cortical_CI[Gordon_left_ROI_positions])

			#cal LH and RH Q using parition restricted to a singel hemisphere
			right_Q_indhemici = cal_modularity_w_imposed_community(bct.threshold_proportional(right_M, p, copy=True), r_CI)
			left_Q_indhemici = cal_modularity_w_imposed_community(bct.threshold_proportional(left_M, p, copy=True), l_CI)
			
			#organize output
			tmp_df.set_value(i, 'SubjID', s)
			tmp_df.set_value(i, 'Density', p)
			tmp_df.set_value(i, 'indwbci_Q', Q_indwbci)
			tmp_df.set_value(i, 'left_indwbci_Q', left_Q_indwbci)
			tmp_df.set_value(i, 'right_indwbci_Q', right_Q_indwbci)
			tmp_df.set_value(i, 'groupwbci_Q', Q_groupwbci)
			tmp_df.set_value(i, 'left_groupwbci_Q', left_Q_groupwbci)
			tmp_df.set_value(i, 'right_groupwbci_Q', right_Q_groupwbci)
			tmp_df.set_value(i, 'left_indhemici_Q', left_Q_indhemici)
			tmp_df.set_value(i, 'right_indhemici_Q', right_Q_indhemici)			
		
		output_df = output_df.append(tmp_df)
	return output_df	

Patient_AdjMats = pool_patient_adjmats()
Control_AdjMats = pickle.load(open(path_to_data_folder +'/Control_AdjMats_Gordon', 'rb'))

Control_df = cal_modularity_w_indiv_infomap(Control_AdjMats)
Patient_df = cal_modularity_w_indiv_infomap(Patient_AdjMats)

Patient_df.to_csv(path_to_data_folder + '/patient_df.csv', index = False)
Control_df.to_csv(path_to_data_folder + '/control_df.csv', index = False)

################################################################
##### cal Z-score
################################################################
def cal_z (pdf, cdf, metric):
	z_score = (pdf[metric] -cdf.groupby(['Density']).aggregate(np.nanmean).reset_index()
			[metric]) / cdf.groupby(['Density']).aggregate(np.nanstd).reset_index()[metric]
	return z_score


#### write out patient dataframe

#thalamic_patients = ['176', '128', '168', '163'] 

Q_df = pd.DataFrame()
for s in np.unique(Patient_df['SubjID']):
	tmp_df = pd.DataFrame()
	Q = cal_z(Patient_df[Patient_df['SubjID']==int(s)], Control_df, 'groupwbci_Q')
	L_Q = cal_z(Patient_df[Patient_df['SubjID']==int(s)], Control_df, 'left_groupwbci_Q')
	R_Q = cal_z(Patient_df[Patient_df['SubjID']==int(s)], Control_df, 'right_groupwbci_Q')
	for i,d in enumerate(np.arange(0.01, 0.16, 0.01)):
		tmp_df.set_value(i, 'SubjID', 'S'+str(int(s)+1))
		tmp_df.set_value(i, 'Density', d)
		tmp_df.set_value(i, 'Whole Brain Q', Q[i])
		if int(s) == 0: #176 bilateral lesion
			tmp_df.set_value(i, 'Lesioned Hemisphere Q', np.nan)
			tmp_df.set_value(i, 'Intact Hemisphere Q', np.nan)
		elif int(s) == 2: #168 lesion on the right
			tmp_df.set_value(i, 'Lesioned Hemisphere Q', R_Q[i])
			tmp_df.set_value(i, 'Intact Hemisphere Q', L_Q[i])
		else:
			tmp_df.set_value(i, 'Lesioned Hemisphere Q', L_Q[i])
			tmp_df.set_value(i, 'Intact Hemisphere Q', R_Q[i])
		tmp_df.set_value(i, 'PC Damage Score', Lesioned_PC[thalamic_patients[int(s)]][i])
	Q_df = Q_df.append(tmp_df)	

Q_df.to_csv(path_to_data_folder + '/Q_df.csv', index = False)


################################################################
###### compare partitions
################################################################
from sklearn.metrics.cluster import adjusted_mutual_info_score as aNMI
Cortical_CI = np.loadtxt(path_to_ROIs + '/Gordon_consensus_CI')

def compare_partitions():
	for p in thalamic_patients:
		fn = '/home/despoB/kaihwang/Rest/NotBackedUp/ParMatrices/Tha_%s_Gordon_333_cortical_corrmat' %p
		M = np.loadtxt(fn)
		g = recursive_network_partition(matrix=M, min_cost=.01,max_cost=0.15, min_community_size=5 ,min_weight=2)
		print(str(p) +": "+str(aNMI(Cortical_CI, gS1_163.community.membership))

	nmi = []
	for s in np.arange(0,62):
		fn = '/home/despoB/connectome-thalamus/Graph/NKI_1400_%s_graph.output' %s
		g = pickle.load(open(fn, "rb"))
		nmi += [aNMI(Cortical_CI, g.community.membership)]
	return nmi	

#compare_partitions()	









