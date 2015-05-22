from __future__ import division, print_function
import numpy as np
import scipy as sp
import scipy.io as sio
import matplotlib.pyplot as plt
import pickle
import csv
import glob
import pandas as pd
import bct
import os
#from ggplot import *

# import functions
import FuncParcel

# vector of cortical ROI index
Cortical_ROIs = np.loadtxt('Data/Cortical_ROI_index')
Cortical_CI = np.loadtxt('Data/Cortical_CI')

# list of subjects
Control_Subj = ['1103', '1220', '1306', '1223', '1314', '1311', '1318', '1313', '1326', '1325', '1328', '1329', '1333', '1331', '1335', '1338', '1336', '1339', '1337', '1344', '1340']
#Control_Subj = ['114', '116', '117', '118', '119', '201', '203', '204', '205', '206', '207', '208', '209', '210', '211', '212', '213', '214', '215', '216', '217', '219', '220']
thalamic_patients = ['128', '162', '163', '168', '176']
striatal_patients = ['b117', 'b122', 'b138', 'b143', 'b153']
patients = thalamic_patients + striatal_patients

if True:
	# create control's dataframe
	OlderControlGlobalData = pd.DataFrame()
	OlderControlNodalData = pd.DataFrame()
	# load old controls
	for s in Control_Subj:
		fn = '/home/despoB/kaihwang/Rest/Graph/gsetCI_%s.mat' %s
		GlobalData, NodalData = FuncParcel.convert_matlab_graph_str(fn, s, Cortical_ROIs)
		OlderControlGlobalData = OlderControlGlobalData.append(GlobalData)
		OlderControlNodalData = OlderControlNodalData.append(NodalData)
	OlderControlNodalData['Group'] = 'Control'
	OlderControlGlobalData['Group'] = 'Control'

	#load young controls

	# convert patients' global graph metrics into z-score
	PatientsGlobalData = pd.DataFrame()
	PatientsNodalData = pd.DataFrame()
	for p in thalamic_patients:
		fn = '/home/despoB/kaihwang/Rest/Graph/gsetCI_%s.mat' %p
		GlobalData, NodalData = FuncParcel.convert_matlab_graph_str(fn, p, Cortical_ROIs)
		GlobalData['Q_zscore'] = FuncParcel.cal_graph_z_score(GlobalData, OlderControlGlobalData, 'Q')
		GlobalData['CC_zscore'] = FuncParcel.cal_graph_z_score(GlobalData, OlderControlGlobalData, 'CC')
		GlobalData['Group'] = 'Thalamic_patients'
		NodalData['PC_zscore'] = FuncParcel.cal_graph_z_score(NodalData, OlderControlNodalData, 'PC')
		NodalData['localE_zscore'] = FuncParcel.cal_graph_z_score(NodalData, OlderControlNodalData, 'localE')
		NodalData['Between_Module_Weight_zscore'] = FuncParcel.cal_graph_z_score(NodalData, OlderControlNodalData, 'Between_Module_Weight')
		NodalData['Within_Module_Weight_zscore'] = FuncParcel.cal_graph_z_score(NodalData, OlderControlNodalData, 'Within_Module_Weight')
		NodalData['Between_Module_Degree_zscore'] = FuncParcel.cal_graph_z_score(NodalData, OlderControlNodalData, 'Between_Module_Degree')
		NodalData['Within_Module_Degree_zscore'] = FuncParcel.cal_graph_z_score(NodalData, OlderControlNodalData, 'Within_Module_Degree')
		NodalData['Group'] = 'Thalamic_patients'
		PatientsGlobalData = PatientsGlobalData.append(GlobalData)
		PatientsNodalData = PatientsNodalData.append(NodalData)

	for p in striatal_patients:
		fn = '/home/despoB/kaihwang/Rest/Graph/gsetCI_%s.mat' %p
		GlobalData, NodalData = FuncParcel.convert_matlab_graph_str(fn, p, Cortical_ROIs)
		GlobalData['Q_zscore'] = FuncParcel.cal_graph_z_score(GlobalData, OlderControlGlobalData, 'Q')
		GlobalData['CC_zscore'] = FuncParcel.cal_graph_z_score(GlobalData, OlderControlGlobalData, 'CC')
		GlobalData['Group'] = 'Striatal_patients'
		NodalData['PC_zscore'] = FuncParcel.cal_graph_z_score(NodalData, OlderControlNodalData, 'PC')
		NodalData['localE_zscore'] = FuncParcel.cal_graph_z_score(NodalData, OlderControlNodalData, 'localE')
		NodalData['Between_Module_Weight_zscore'] = FuncParcel.cal_graph_z_score(NodalData, OlderControlNodalData, 'Between_Module_Weight')
		NodalData['Within_Module_Weight_zscore'] = FuncParcel.cal_graph_z_score(NodalData, OlderControlNodalData, 'Within_Module_Weight')
		NodalData['Between_Module_Degree_zscore'] = FuncParcel.cal_graph_z_score(NodalData, OlderControlNodalData, 'Between_Module_Degree')
		NodalData['Within_Module_Degree_zscore'] = FuncParcel.cal_graph_z_score(NodalData, OlderControlNodalData, 'Within_Module_Degree')
		NodalData['Group'] = 'Striatal_patients'
		PatientsGlobalData = PatientsGlobalData.append(GlobalData)
		PatientsNodalData = PatientsNodalData.append(NodalData)

	#combine dataframes
	GraphGlobalData = pd.DataFrame()
	GraphNodalData = pd.DataFrame()

	GraphGlobalData = PatientsGlobalData.append(OlderControlGlobalData)
	GraphNodalData = PatientsNodalData.append(OlderControlNodalData)

	GraphGlobalData.to_csv('Data/GraphGlobalData.csv')
	GraphNodalData.to_csv('Data/GraphNodalData.csv')



#GraphGlobalData = pd.DataFrame.from_csv('Data/GraphGlobalData.csv')

# calculate hemispheric difference
GraphGlobalData['RightvLeft_Q'] = GraphGlobalData['right_Q'] - GraphGlobalData['left_Q']
GraphGlobalData['LeftvRight_Q'] = GraphGlobalData['left_Q'] - GraphGlobalData['right_Q']
GraphGlobalData.to_csv('Data/GraphGlobalData.csv')


# use BrainX to get MGH subjects partition and nodal roles...
file_path = '/home/despoB/kaihwang/Rest/AdjMatrices/*Ses1_Full_WashU333_corrmat'
AveMat = FuncParcel.average_corrmat(file_path)
np.savetxt('Data/CorticalAveMat', AveMat)

# thresholding
#m = bct.binarize(bct.threshold_proportional(AveMat, .05))

# run partiion across threshold on the MGH's avemat
MGH_template_partition = pd.DataFrame()
row_count = 0
for d in np.unique(GraphNodalData.Density): #extracting density vectors from previosly created density values to avoid floating error....
	ci, q = bct.modularity_und(bct.binarize(bct.threshold_proportional(AveMat, d)))
	for roi in np.arange(0, len(ci)):
		MGH_template_partition.loc[row_count, 'Ci'] = ci[roi]
		MGH_template_partition.loc[row_count, 'ROI'] = Cortical_ROIs[roi]
		MGH_template_partition.loc[row_count, 'Density'] = d
		row_count = row_count +1
MGH_template_partition.to_csv('Data/MGH_partition.csv')

#run parition on subject data
Subject_partition = pd.DataFrame()
row_count = 0
for d in np.unique(GraphNodalData.Density): #extracting density vectors from previosly created density values to avoid floating error....
	for s in patients+Control_Subj:
		fn = '/home/despoB/kaihwang/Rest/AdjMatrices/t%s_Full_WashU333_corrmat' %s
		Mat = np.loadtxt(fn) 
		ci, q = bct.modularity_und(bct.binarize(bct.threshold_proportional(Mat, d)))
		for roi in np.arange(0, len(ci)):
			Subject_partition.loc[row_count, 'Ci'] = ci[roi]
			Subject_partition.loc[row_count, 'ROI'] = Cortical_ROIs[roi]
			Subject_partition.loc[row_count, 'Density'] = d
			Subject_partition.loc[row_count, 'Subject'] = s
			row_count = row_count +1
Subject_partition.to_csv('Data/Subject_partition.csv')

# calculate mutual information
from sklearn.metrics import normalized_mutual_info_score

NMI_dataframe = pd.DataFrame()
row_count = 0
for d in np.unique(GraphNodalData.Density):
	for s in patients+Control_Subj:
		tmp_df = Subject_partition[Subject_partition.Density==d][['Subject','ROI','Ci']]
		subject_ci = tmp_df[tmp_df.Subject==s]['Ci'].values
		template_ci = MGH_template_partition[MGH_template_partition.Density == d]['Ci'].values
		NMI_dataframe.loc[row_count,'NMI'] = normalized_mutual_info_score(subject_ci, template_ci)
		NMI_dataframe.loc[row_count,'Subject'] = s
		NMI_dataframe.loc[row_count,'Density'] = d
		row_count = row_count+1
NMI_dataframe.to_csv('Data/NMI_dataframe.csv')


#nmi = normalized_mutual_info_score(ci, ci)



# explore ggplot
#ggplot(aes(x = 'Density', y = 'Q_zscore', color = 'Subject'),data = GraphGlobalData[GraphGlobalData.Group == 'Thalamic_patients']) \
#+ geom_line() + xlim(0.05, 0.15)

#cleanup
#reset_selective GlobalData
#reset_selective NodalData
#reset_selective OlderControlGlobalData
#reset_selective OlderControlNodalData


