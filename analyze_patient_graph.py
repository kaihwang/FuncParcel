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
from ggplot import *

# import functions
import FuncParcel

# vector of cortical ROI index
Cortical_ROIs = np.loadtxt('Data/Cortical_ROI_index')

# list of subjects
Control_Subj = ['1103', '1220', '1306', '1223', '1314', '1311', '1318', '1313', '1326', '1325', '1328', '1329', '1333', '1331', '1335', '1338', '1336', '1339', '1337', '1344', '1340']
thalamic_patients = ['128', '162', '163', '168', '176']
striatal_patients = ['b117', 'b122', 'b138', 'b143', 'b153']
patients = thalamic_patients + striatal_patients

if os.path.isfile('data/GraphGlobalData.csv') & os.path.isfile('data/GraphNodalData.csv') is False:
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



GraphGlobalData = pd.DataFrame.from_csv('data/GraphGlobalData.csv')

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

# run partiion across threshold on the ave mat
MGH_template_partition = pd.DataFrame()
row_count = 0
for d in np.arange(0.01, 0.255, 0.005):
	ci, q = bct.modularity_louvain_und(bct.binarize(bct.threshold_proportional(AveMat, d)))
	for roi in np.arange(0, len(ci)):
		MGH_template_partition.loc[row_count, 'Ci'] = ci[roi]
		MGH_template_partition.loc[row_count, 'ROI'] = Cortical_ROIs[roi]
		MGH_template_partition.loc[row_count, 'Density'] = d
		row_count = row_count +1

# calculate mutual information
#from sklearn.metrics import normalized_mutual_info_score
#nmi = normalized_mutual_info_score(ci, ci)



# explore ggplot
#ggplot(aes(x = 'Density', y = 'Q_zscore', color = 'Subject'),data = GraphGlobalData[GraphGlobalData.Group == 'Thalamic_patients']) \
#+ geom_line() + xlim(0.05, 0.15)

#cleanup
#reset_selective GlobalData
#reset_selective NodalData
#reset_selective OlderControlGlobalData
#reset_selective OlderControlNodalData


