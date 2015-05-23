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
import surfer
from surfer import Brain
import bct
import brewer2mpl
import FuncParcel

# what to do?
calculate_z_scores = True
calulate_template_partition = False
visuazlie_template_partition = False
visualize_patient_cortical_target = False
visualize_hubs = False

# vector of cortical ROI index
Cortical_ROIs = np.loadtxt('Data/Cortical_ROI_index')
Cortical_CI = np.loadtxt('Data/Cortical_CI')

# list of subjects
Control_Subj = ['1103', '1220', '1306', '1223', '1314', '1311', '1318', '1313', '1326', '1325', '1328', '1329', '1333', '1331', '1335', '1338', '1336', '1339', '1337', '1344', '1340']
#Control_Subj = ['114', '116', '117', '118', '119', '201', '203', '204', '205', '206', '207', '208', '209', '210', '211', '212', '213', '214', '215', '216', '217', '219', '220']
thalamic_patients = ['128', '162', '163', '168', '176']
striatal_patients = ['b117', 'b122', 'b138', 'b143', 'b153']
patients = thalamic_patients + striatal_patients

if calculate_z_scores:
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

	# calculate hemispheric difference
	GraphGlobalData['RightvLeft_Q'] = GraphGlobalData['right_Q'] - GraphGlobalData['left_Q']
	GraphGlobalData['LeftvRight_Q'] = GraphGlobalData['left_Q'] - GraphGlobalData['right_Q']
	GraphGlobalData.to_csv('Data/GraphGlobalData.csv')


# thresholding
#m = bct.binarize(bct.threshold_proportional(AveMat, .05))

# run partiion across threshold on the MGH's avemat
# MGH_template_partition = pd.DataFrame()
# row_count = 0
# for d in np.unique(GraphNodalData.Density): #extracting density vectors from previosly created density values to avoid floating error....
# 	ci, q = bct.modularity_und(bct.binarize(bct.threshold_proportional(AveMat, d)))
# 	for roi in np.arange(0, len(ci)):
# 		MGH_template_partition.loc[row_count, 'Ci'] = ci[roi]
# 		MGH_template_partition.loc[row_count, 'ROI'] = Cortical_ROIs[roi]
# 		MGH_template_partition.loc[row_count, 'Density'] = d
# 		row_count = row_count +1
# MGH_template_partition.to_csv('Data/MGH_partition.csv')


# #run parition on subject data
# Subject_partition = pd.DataFrame()
# row_count = 0
# for d in np.unique(GraphNodalData.Density): #extracting density vectors from previosly created density values to avoid floating error....
# 	for s in patients+Control_Subj:
# 		fn = '/home/despoB/kaihwang/Rest/AdjMatrices/t%s_Full_WashU333_corrmat' %s
# 		Mat = np.loadtxt(fn) 
# 		ci, q = bct.modularity_und(bct.binarize(bct.threshold_proportional(Mat, d)))
# 		for roi in np.arange(0, len(ci)):
# 			Subject_partition.loc[row_count, 'Ci'] = ci[roi]
# 			Subject_partition.loc[row_count, 'ROI'] = Cortical_ROIs[roi]
# 			Subject_partition.loc[row_count, 'Density'] = d
# 			Subject_partition.loc[row_count, 'Subject'] = s
# 			row_count = row_count +1
# Subject_partition.to_csv('Data/Subject_partition.csv')

# # calculate mutual information
# from sklearn.metrics import normalized_mutual_info_score

# NMI_dataframe = pd.DataFrame()
# row_count = 0
# for d in np.unique(GraphNodalData.Density):
# 	for s in patients+Control_Subj:
# 		tmp_df = Subject_partition[Subject_partition.Density==d][['Subject','ROI','Ci']]
# 		subject_ci = tmp_df[tmp_df.Subject==s]['Ci'].values
# 		template_ci = MGH_template_partition[MGH_template_partition.Density == d]['Ci'].values
# 		NMI_dataframe.loc[row_count,'NMI'] = normalized_mutual_info_score(subject_ci, template_ci)
# 		NMI_dataframe.loc[row_count,'Subject'] = s
# 		NMI_dataframe.loc[row_count,'Density'] = d
# 		row_count = row_count+1
# NMI_dataframe.to_csv('Data/NMI_dataframe.csv')


# gent template partion and nodal properties from MGH data
if calulate_template_partition:
	AveMat = np.loadtxt('Data/CorticalAveMat')

	template_ci, template_q = bct.modularity_und(bct.binarize(bct.threshold_proportional(AveMat, 0.052))) #threshold at 0.075 cost
	template_pc = bct.participation_coef(bct.binarize(bct.threshold_proportional(AveMat, 0.052)), template_ci)
	template_wmd = bct.module_degree_zscore(bct.binarize(bct.threshold_proportional(AveMat, 0.052)), template_ci)

	Cortical_ROI_Coordinate = np.loadtxt('Data/Cortical_ROI_Coordinate')
	template_nodal_data = pd.DataFrame()
	template_nodal_data['ROI'] = Cortical_ROIs
	template_nodal_data['PC'] = template_pc
	template_nodal_data['WMD'] = template_wmd
	template_nodal_data['Ci'] = template_ci
	#template_nodal_data['Coordinate'] = Cortical_ROI_Coordinate
	template_nodal_data.to_csv('Data/template_nodal_data.csv')

	connector_hubs =  template_nodal_data.ROI[template_nodal_data.PC>0.65].values
	connector_hubs = connector_hubs.astype(int)
	np.savetxt('Data/connector_hubs', connector_hubs)

	provincial_hubs =  template_nodal_data.ROI[template_nodal_data.WMD>1.2].values
	provincial_hubs = provincial_hubs.astype(int)
	np.savetxt('Data/provincial_hubs', provincial_hubs)

	both_hubs = np.intersect1d(provincial_hubs,connector_hubs)
	np.savetxt('Data/both_hubs', both_hubs)


#try to visulize template graph partition
if visuazlie_template_partition:
	subjects_dir = os.environ["SUBJECTS_DIR"]
	subject_id, surface = 'fsaverage', 'inflated'
	hemi = 'split'
	brain = Brain(subject_id, hemi, surface,
	              config_opts=dict(width=800, height=400, background="white"))

	#bmap = brewer2mpl.get_map('Paired', 'Qualitative', 12)
	colors = ['#00ffff', '#000000', '#0000ff', '#ff00ff', '#008000', '#808080', '#00ff00', '#800000', '#000080', '#808000', '#800080', '#ff0000', '#c0c0c0', '#008080', '#ffffff', '#ffff00']
	c_i = 0
	for i in range(1, int(np.max(template_ci)+1)):
		coor = Cortical_ROI_Coordinate[template_ci==i]
		print(coor)
		if len(coor)>1:
			brain.add_foci(coor[coor[:,0]<0], map_surface="white", color=colors[c_i], hemi="lh" )
			brain.add_foci(coor[coor[:,0]>0], map_surface="white", color=colors[c_i], hemi="rh" )
		c_i = c_i+1



# visualize hubs
if visualize_hubs:
	subjects_dir = os.environ["SUBJECTS_DIR"]
	subject_id, surface = 'fsaverage', 'inflated'
	hemi = 'split'
	brain2 = Brain(subject_id, hemi, surface,
	              config_opts=dict(width=800, height=400, background="white"))
	coor = Cortical_ROI_Coordinate[template_nodal_data['ROI'].isin(connector_hubs).values]
	brain2.add_foci(coor[coor[:,0]<0], map_surface="white", color='red', hemi="lh" )
	brain2.add_foci(coor[coor[:,0]>0], map_surface="white", color='red', hemi="rh" )
	coor = Cortical_ROI_Coordinate[template_nodal_data['ROI'].isin(provincial_hubs).values]
	brain2.add_foci(coor[coor[:,0]<0], map_surface="white", color='blue', hemi="lh" )
	brain2.add_foci(coor[coor[:,0]>0], map_surface="white", color='blue', hemi="rh" )
	coor = Cortical_ROI_Coordinate[template_nodal_data['ROI'].isin(both_hubs).values]
	brain.add_foci(coor[coor[:,0]<0], map_surface="white", color='black', hemi="lh" )
	brain.add_foci(coor[coor[:,0]>0], map_surface="white", color='black', hemi="rh" )

#visualize patients cortical target
if visualize_patient_cortical_target:
	subjects_dir = os.environ["SUBJECTS_DIR"]
	subject_id, surface = 'fsaverage', 'inflated'
	hemi = 'split'
	patient_target = np.loadtxt('Data/128_cortical_target')
	patient_target = patient_target.astype(int)
	brain3 = Brain(subject_id, hemi, surface,
	              config_opts=dict(width=800, height=400, background="white"))
	coor = Cortical_ROI_Coordinate[template_nodal_data['ROI'].isin(patient_target).values]
	brain3.add_foci(coor[coor[:,0]<0], map_surface="white", color='red', hemi="lh" )
	brain3.add_foci(coor[coor[:,0]>0], map_surface="white", color='red', hemi="rh" )

# enter patietn's cortical target nodal data and node type
GraphNodalData = pd.DataFrame.from_csv('Data/GraphNodalData.csv')
GraphNodalData['Subject'] = GraphNodalData['Subject'].astype('string')
connector_hubs = np.loadtxt('Data/connector_hubs')
provincial_hubs = np.loadtxt('Data/provincial_hubs')
GraphNodalData['node_type'] = 'none_hub'
GraphNodalData['node_type'].loc[GraphNodalData['ROI'].isin(connector_hubs)] = 'connector_hub'
GraphNodalData['node_type'].loc[GraphNodalData['ROI'].isin(provincial_hubs)] = 'provincial_hub'

GraphNodalData['target'] = False
GraphNodalData['non_target'] = False
for p in patients:
	fn = 'Data/%s_cortical_target' %p 
	patient_target = np.loadtxt(fn)
	fn = 'Data/%s_cortical_nontarget' %p
	patient_nontarget = np.loadtxt(fn) 
	patient_target = patient_target.astype(int)
	patient_nontarget = patient_nontarget.astype(int)
	GraphNodalData['target'].loc[(GraphNodalData['Subject']==p) & (GraphNodalData['ROI'].isin(patient_target))] = True
	GraphNodalData['non_target'].loc[(GraphNodalData['Subject']==p) & (GraphNodalData['ROI'].isin(patient_nontarget))] = True
GraphNodalData.to_csv('Data/GraphNodalData.csv')

# explore ggplot
#ggplot(aes(x = 'Density', y = 'Q_zscore', color = 'Subject'),data = GraphGlobalData[GraphGlobalData.Group == 'Thalamic_patients']) \
#+ geom_line() + xlim(0.05, 0.15)

#cleanup
#reset_selective GlobalData
#reset_selective NodalData
#reset_selective OlderControlGlobalData
#reset_selective OlderControlNodalData


