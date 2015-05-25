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
#import brewer2mpl
import FuncParcel
#from ggplot import *

# what to do?
calculate_z_scores = True
calulate_template_partition = True
visuazlie_template_partition = False
visualize_patient_cortical_target = False
visualize_hubs = False
organize_patient_target_nontarget_data = False
run_template_partition_across_densities = False
cal_sub_parition_by_densities = False
cal_NMI = False

# vector of cortical ROI index
Cortical_ROIs = np.loadtxt('Data/Cortical_ROI_index')
Cortical_CI = np.loadtxt('Data/Cortical_CI')

# list of subjects
Control_Subj = ['1103', '1220', '1306', '1223', '1314', '1311', '1318', '1313', '1326', '1325', '1328', '1329', '1333', '1331', '1335', '1338', '1336', '1339', '1337', '1344', '1340']
#Control_Subj = ['114', '116', '117', '118', '119', '201', '203', '204', '205', '206', '207', '208', '209', '210', '211', '212', '213', '214', '215', '216', '217', '219', '220']
thalamic_patients = ['128', '162', '163', '168', '176']
striatal_patients = ['b117', 'b122', 'b138', 'b143', 'b153']
patients = thalamic_patients + striatal_patients


# get template partion and nodal properties from MGH data
if calulate_template_partition:
	AveMat = np.loadtxt('Data/CorticalAveMat')

	#template_ci, template_q = bct.modularity_und(bct.binarize(bct.threshold_proportional(AveMat, 0.05))) #threshold at 0.05 cost
	template_ci = np.loadtxt('Data/MGH_CI') #use previously generated CI at .05 cost
	template_ci = template_ci.astype(int)
	template_pc = bct.participation_coef(bct.binarize(bct.threshold_proportional(AveMat, 0.05)), template_ci)
	template_wmd = bct.module_degree_zscore(bct.binarize(bct.threshold_proportional(AveMat, 0.05)), template_ci)

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

	# load patients' graph metrics, convert global metrics into z_score 
	PatientsGlobalData = pd.DataFrame()
	PatientsNodalData = pd.DataFrame()
	for p in thalamic_patients:
		fn = '/home/despoB/kaihwang/Rest/Graph/gsetCI_%s.mat' %p
		GlobalData, NodalData = FuncParcel.convert_matlab_graph_str(fn, p, Cortical_ROIs)
		GlobalData['Q_zscore'] = FuncParcel.cal_graph_z_score(GlobalData, OlderControlGlobalData, np.array([]), 'Q')
		GlobalData['CC_zscore'] = FuncParcel.cal_graph_z_score(GlobalData, OlderControlGlobalData, np.array([]), 'CC')
		GlobalData['Group'] = 'Thalamic_patients'
		NodalData['Group'] = 'Thalamic_patients'
		PatientsGlobalData = PatientsGlobalData.append(GlobalData)
		PatientsNodalData = PatientsNodalData.append(NodalData)

	for p in striatal_patients:
		fn = '/home/despoB/kaihwang/Rest/Graph/gsetCI_%s.mat' %p
		GlobalData, NodalData = FuncParcel.convert_matlab_graph_str(fn, p, Cortical_ROIs)
		GlobalData['Q_zscore'] = FuncParcel.cal_graph_z_score(GlobalData, OlderControlGlobalData, np.array([]), 'Q')
		GlobalData['CC_zscore'] = FuncParcel.cal_graph_z_score(GlobalData, OlderControlGlobalData, np.array([]), 'CC')
		GlobalData['Group'] = 'Striatal_patients'
		#NodalData['PC_zscore'] = FuncParcel.cal_graph_z_score(NodalData, OlderControlNodalData, True,'PC')
		#NodalData['localE_zscore'] = FuncParcel.cal_graph_z_score(NodalData, OlderControlNodalData, True,'localE')
		NodalData['Group'] = 'Striatal_patients'
		PatientsGlobalData = PatientsGlobalData.append(GlobalData)
		PatientsNodalData = PatientsNodalData.append(NodalData)

	
	#organize nodal properties
	connector_hubs = np.loadtxt('Data/connector_hubs')
	connector_hubs = connector_hubs.astype(int)
	provincial_hubs = np.loadtxt('Data/provincial_hubs')
	provincial_hubs = provincial_hubs.astype(int)

	PatientsNodalData['node_type'] = 'none_hub'
	PatientsNodalData['node_type'].loc[PatientsNodalData['ROI'].isin(connector_hubs)] = 'connector_hub'
	PatientsNodalData['node_type'].loc[PatientsNodalData['ROI'].isin(provincial_hubs)] = 'provincial_hub'

	OlderControlNodalData['node_type'] = 'none_hub'
	OlderControlNodalData['node_type'].loc[OlderControlNodalData['ROI'].isin(connector_hubs)] = 'connector_hub'
	OlderControlNodalData['node_type'].loc[OlderControlNodalData['ROI'].isin(provincial_hubs)] = 'provincial_hub'
	#ControlDf = OlderControlNodalData.groupby(['Subject','ROI', 'Density']).aggregate(np.nanmean).reset_index()

	PatientsNodalData['target'] = False
	PatientsNodalData['non_target'] = False
	for p in patients:
		fn = 'Data/%s_cortical_target' %p 
		patient_target = np.loadtxt(fn)
		fn = 'Data/%s_cortical_nontarget' %p
		patient_nontarget = np.loadtxt(fn) 
		patient_target = patient_target.astype(int)
		patient_nontarget = patient_nontarget.astype(int)
		PatientsNodalData['target'].loc[(PatientsNodalData['Subject']==p) & (PatientsNodalData['ROI'].isin(patient_target))] = True
		PatientsNodalData['non_target'].loc[(PatientsNodalData['Subject']==p) & (PatientsNodalData['ROI'].isin(patient_nontarget))] = True

		
	#convert nodal graph metrics into z_score
	PatientsNodalZscoreData = pd.DataFrame()
	for p in patients:
		fn = 'Data/%s_cortical_target' %p 
		patient_target = np.loadtxt(fn)
		fn = 'Data/%s_cortical_nontarget' %p
		patient_nontarget = np.loadtxt(fn) 
		patient_target = patient_target.astype(int)
		patient_nontarget = patient_nontarget.astype(int)
		graph_metrics = ['Between_Module_Weight', 'Within_Module_Weight', 'localE', 'PC']
		patient_zDF= pd.DataFrame()

		#connector_hub 
		tmp_dict = {}
		target_selector = np.intersect1d(patient_target, connector_hubs)
		nontarget_selector = np.intersect1d(patient_nontarget, connector_hubs)
		if target_selector.any() & nontarget_selector.any():
			for metric in graph_metrics:
				tmp_dict['Target_' + metric] = FuncParcel.cal_graph_z_score(PatientsNodalData[PatientsNodalData['Subject']==p], OlderControlNodalData, target_selector , metric)
				tmp_dict['nonTarget_' + metric] = FuncParcel.cal_graph_z_score(PatientsNodalData[PatientsNodalData['Subject']==p], OlderControlNodalData, nontarget_selector, metric)
				tmp_dict['Density'] = np.linspace(0.01, 0.25, num = 49)
				tmp_dict['Subject'] = [p] * 49
				tmp_dict['node_type'] = ['connector_hub'] * 49

			tmpdf = pd.DataFrame(tmp_dict, columns=['Subject', 'Density', 'node_type', \
				'Target_Between_Module_Weight', 'Target_Within_Module_Weight', 'Target_PC', 'Target_localE',\
				'nonTarget_Between_Module_Weight', 'nonTarget_Within_Module_Weight', 'nonTarget_PC', 'nonTarget_localE'])
			patient_zDF = patient_zDF.append(tmpdf)

		#provincial_hub
		tmp_dict = {}
		target_selector = np.intersect1d(patient_target, provincial_hubs)
		nontarget_selector = np.intersect1d(patient_nontarget, provincial_hubs)
		if target_selector.any() & nontarget_selector.any(): 
			for metric in graph_metrics:
				tmp_dict['Target_' + metric] = FuncParcel.cal_graph_z_score(PatientsNodalData[PatientsNodalData['Subject']==p], OlderControlNodalData, target_selector, metric)
				tmp_dict['nonTarget_' + metric] = FuncParcel.cal_graph_z_score(PatientsNodalData[PatientsNodalData['Subject']==p], OlderControlNodalData, nontarget_selector, metric)
				tmp_dict['Density'] = np.linspace(0.01, 0.25, num = 49)
				tmp_dict['Subject'] = [p] * 49
				tmp_dict['node_type'] = ['provincial_hub'] * 49

			tmpdf = pd.DataFrame(tmp_dict, columns=['Subject', 'Density', 'node_type', \
				'Target_Between_Module_Weight', 'Target_Within_Module_Weight', 'Target_PC', 'Target_localE',\
				'nonTarget_Between_Module_Weight', 'nonTarget_Within_Module_Weight', 'nonTarget_PC', 'nonTarget_localE'])
			patient_zDF = patient_zDF.append(tmpdf)	
			
		#non_hub
		target_selector = np.array(list(set(patient_nontarget) - set(np.intersect1d(patient_nontarget, connector_hubs)) - set(np.intersect1d(patient_nontarget, provincial_hubs))))
		nontarget_selector = np.array(list(set(patient_nontarget) - set(np.intersect1d(patient_nontarget, connector_hubs)) - set(np.intersect1d(patient_nontarget, provincial_hubs))))
		tmp_dict = {}
		for metric in graph_metrics:
			tmp_dict['Target_' + metric] = FuncParcel.cal_graph_z_score(PatientsNodalData[PatientsNodalData['Subject']==p], OlderControlNodalData, target_selector , metric)
			tmp_dict['nonTarget_' + metric] = FuncParcel.cal_graph_z_score(PatientsNodalData[PatientsNodalData['Subject']==p], OlderControlNodalData, nontarget_selector, metric)
			tmp_dict['Density'] = np.linspace(0.01, 0.25, num = 49)
			tmp_dict['Subject'] = [p] * 49
			tmp_dict['node_type'] = ['non_hub'] * 49

		tmpdf = pd.DataFrame(tmp_dict, columns=['Subject', 'Density', 'node_type', \
			'Target_Between_Module_Weight', 'Target_Within_Module_Weight', 'Target_PC', 'Target_localE',\
			'nonTarget_Between_Module_Weight', 'nonTarget_Within_Module_Weight', 'nonTarget_PC', 'nonTarget_localE'])
		patient_zDF = patient_zDF.append(tmpdf)
			
		#all
		tmp_dict = {}
		for metric in graph_metrics:
			tmp_dict['Target_' + metric] = FuncParcel.cal_graph_z_score(PatientsNodalData[PatientsNodalData['Subject']==p], OlderControlNodalData, patient_target , metric)
			tmp_dict['nonTarget_' + metric] = FuncParcel.cal_graph_z_score(PatientsNodalData[PatientsNodalData['Subject']==p], OlderControlNodalData, patient_nontarget , metric)
			tmp_dict['Density'] = np.linspace(0.01, 0.25, num = 49)
			tmp_dict['Subject'] = [p] * 49
			tmp_dict['node_type'] = ['all'] * 49
	
		tmpdf = pd.DataFrame(tmp_dict, columns=['Subject', 'Density', 'node_type', \
			'Target_Between_Module_Weight', 'Target_Within_Module_Weight', 'Target_PC', 'Target_localE',\
			'nonTarget_Between_Module_Weight', 'nonTarget_Within_Module_Weight', 'nonTarget_PC', 'nonTarget_localE'])
		patient_zDF = patient_zDF.append(tmpdf)
		
		PatientsNodalZscoreData = PatientsNodalZscoreData.append(patient_zDF)

	PatientsNodalZscoreData['Group']='na'
	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='128'] = 'Thalamic_patients'
	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='162'] = 'Thalamic_patients'
	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='163'] = 'Thalamic_patients'
	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='168'] = 'Thalamic_patients'
	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='176'] = 'Thalamic_patients'
	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='b117'] = 'Striatal_patients'
	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='b122'] = 'Striatal_patients'
	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='b138'] = 'Striatal_patients'
	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='b143'] = 'Striatal_patients'
	PatientsNodalZscoreData['Group'].loc[PatientsNodalZscoreData['Subject']=='b153'] = 'Striatal_patients'
	PatientsNodalZscoreData.to_csv('Data/PatientsNodalZscoreData.csv')		

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



#run partiion across threshold on the MGH's avemat across a range of densities
if run_template_partition_across_densities:
	GraphNodalData = pd.DataFrame.from_csv('Data/GraphNodalData.csv')
	GraphNodalData['Subject'] = GraphNodalData['Subject'].astype('string').values
	AveMat = np.loadtxt('Data/CorticalAveMat')
	MGH_template_partition = pd.DataFrame()
	row_count = 0
	for d in np.unique(GraphNodalData.Density): #extracting density vectors from previosly created density values to avoid floating error....
		ci, q = bct.modularity_louvain_und(bct.binarize(bct.threshold_proportional(AveMat, d)))
		for roi in np.arange(0, len(ci)):
			MGH_template_partition.loc[row_count, 'Ci'] = ci[roi]
			MGH_template_partition.loc[row_count, 'ROI'] = Cortical_ROIs[roi]
			MGH_template_partition.loc[row_count, 'Density'] = d
			row_count = row_count +1
	MGH_template_partition.to_csv('Data/MGH_partition.csv')


# #run parition on subject data
if cal_sub_parition_by_densities:
	Subject_partition = pd.DataFrame()
	row_count = 0
	GraphNodalData = pd.DataFrame.from_csv('Data/GraphNodalData.csv')
	GraphNodalData['Subject'] = GraphNodalData['Subject'].astype('string').values
	for d in np.unique(GraphNodalData.Density): #extracting density vectors from previosly created density values to avoid floating error....
		for s in patients+Control_Subj:
			fn = '/home/despoB/kaihwang/Rest/AdjMatrices/t%s_Full_WashU333_corrmat' %s
			Mat = np.loadtxt(fn) 
			ci, q = bct.modularity_louvain_und(bct.binarize(bct.threshold_proportional(Mat, d)))
			for roi in np.arange(0, len(ci)):
				Subject_partition.loc[row_count, 'Ci'] = ci[roi]
				Subject_partition.loc[row_count, 'ROI'] = Cortical_ROIs[roi]
				Subject_partition.loc[row_count, 'Density'] = d
				Subject_partition.loc[row_count, 'Subject'] = s
				Subject_partition.loc[row_count, 'Group'] = GraphNodalData[GraphNodalData['Subject']==s]['Group'].values[0]
				row_count = row_count +1
	Subject_partition.to_csv('Data/Subject_partition.csv')

# # calculate mutual information
if cal_NMI:
	from sklearn.metrics import normalized_mutual_info_score
	MGH_template_partition = pd.DataFrame.from_csv('Data/MGH_partition.csv')
	GraphNodalData = pd.DataFrame.from_csv('Data/GraphNodalData.csv')
	GraphNodalData['Subject'] = GraphNodalData['Subject'].astype('string').values
	NMI_dataframe = pd.DataFrame()
	row_count = 0
	for d in np.unique(GraphNodalData.Density):
		for s in patients+Control_Subj:
			tmp_df = GraphNodalData[GraphNodalData.Density==d][['Subject','ROI','Ci', 'Group']]
			subject_ci = tmp_df[tmp_df.Subject==s]['Ci'].values
			template_ci = MGH_template_partition[MGH_template_partition.Density == d]['Ci'].values
			subject_ci = subject_ci.astype(int)
			template_ci = template_ci.astype(int)
			NMI_dataframe.loc[row_count,'NMI'] = normalized_mutual_info_score(subject_ci, template_ci)
			NMI_dataframe.loc[row_count,'Subject'] = s
			NMI_dataframe.loc[row_count,'Density'] = d
			NMI_dataframe.loc[row_count,'Group'] = tmp_df[tmp_df['Subject']==s]['Group'].values[0]
			row_count = row_count+1
	NMI_dataframe.to_csv('Data/NMI_dataframe.csv')


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



# plot patient nodal data
#plotData = pd.DataFrame.from_csv('PatientsNodalZscoreData.csv')
#plotData = GraphNodalData.loc[GraphNodalData['target']==True].groupby(['Subject', 'Density','node_type']).aggregate(np.nanmean).reset_index()
#ggplot(aes(x = 'Density', y = 'Target_localE', color = 'Subject'),data = plotData) \
#+ geom_line() + xlim(0.05, 0.15) + facet_wrap('node_type', ncol = 4)

#cleanup
#reset_selective GlobalData
#reset_selective NodalData
#reset_selective OlderControlGlobalData
#reset_selective OlderControlNodalData


