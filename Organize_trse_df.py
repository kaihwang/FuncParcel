from FuncParcel import *
import pandas as pd
#from ggplot import *
import pickle


path_to_data_folder = '/home/despoB/kaihwang/bin/FuncParcel/Data/'
Data_path='/home/despoB/kaihwang/TRSE/TRSEPPI/Graph/'
Conditions = ['categ_face', 'categ_scene', 'relev_face', 'relev_scene', 'irrel_face', 'irrel_scene']
ROIs = ['Gordon_plus_Morel', 'Gordon_plus_thalamus_WTA']

Morel_ind = {'VM':348, 'VA':347, 'VP':346, 'LGN': 345, 'PuA':344, 'PuL':343, 'PuI':342, 'PuM':341, 'Po':340, 'MGN':339, 'IL':338, 'LP':337, 'MD':336,'VL':335,'AN':334}
Func_ind ={'sFP':341, 'T':340, 'mT':339, 'mO':338, 'latO':337, 'FP':336, 'SM':335, 'CO':334, 'DM':333}

INDEX = [Morel_ind, Func_ind]

# load behav data from Courtney
RT_df = pd.read_csv(path_to_data_folder+'bseries_trse_RT.csv')
Accu_df = pd.read_csv(path_to_data_folder+'bseries_trse_Accu.csv')


for condition in Conditions:
	for roi, ind in zip(ROIs,INDEX):
		
		for key in ind:
			RT_df[condition+'_'+key] = 0
			Accu_df[condition+'_'+key] = 0
		
		for s in RT_df['sub']:
			pc = np.loadtxt(Data_path+str(s)+'_'+roi+'_'+condition+'_meanPC')
			
			for key in ind:		
				RT_df.loc[RT_df['sub']==s, condition+'_'+key] = pc[ind[key]] 
				Accu_df.loc[Accu_df['sub']==s, condition+'_'+key] = pc[ind[key]] 

RT_df.to_csv(path_to_data_folder+'bseries_trse_RT.csv')
Accu_df.to_csv(path_to_data_folder+'bseries_trse_Accu.csv')