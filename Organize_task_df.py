#organized data for TRSE, TDSig, HCP
from FuncParcel import *
import pandas as pd
import pickle


#common variables
path_to_data_folder = '/home/despoB/kaihwang/bin/FuncParcel/Data/'
Data_path='/home/despoB/kaihwang/Rest/Graph/'
ROIs = ['Gordon_plus_Morel_LPI', 'Gordon_plus_thalamus_WTA_LPI']
Morel_ind = {'VM':348, 'VA':347, 'VP':346, 'LGN': 345, 'PuA':344, 'PuL':343, 'PuI':342, 'PuM':341, 'Po':340, 'MGN':339, 'IL':338, 'LP':337, 'MD':336,'VL':335,'AN':334}
Func_ind ={'sFP':341, 'T':340, 'mT':339, 'mO':338, 'latO':337, 'FP':336, 'SM':335, 'CO':334, 'DM':333}
INDEX = [Morel_ind, Func_ind]

# load behav data from Courtney
def organize_TRSE():
	RT_df = pd.read_csv(path_to_data_folder+'bseries_trse_RT.csv')
	Accu_df = pd.read_csv(path_to_data_folder+'bseries_trse_Accu.csv')
	TRSE_Conditions = ['Cat', 'FH', 'HF', 'Both']
	for condition in TRSE_Conditions:
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


def organize_TDSigEI():	
	df = pd.read_csv(path_to_data_folder+'TDSigEI_behav.csv')
	TD_Conditions = ['FH', 'HF', 'Fp', 'Hp']
	for condition in TD_Conditions:
		for roi, ind in zip(ROIs,INDEX):
			
			for key in ind:
				df[condition+'_'+key] = 0
			
			for s in df['sub']:
				pc = np.loadtxt(Data_path+str(s)+'_'+roi+'_'+condition+'_meanPC')
				
				for key in ind:		
					df.loc[df['sub']==s, condition+'_'+key] = pc[ind[key]] 
	df.to_csv(path_to_data_folder+'TDSigEI_behav.csv')			


def get_HCP_task_performance(subject,task):
	#this is from Maxwell
	bdf = pd.read_csv('/home/despoB/mb3152/dynamic_mod/os_behavior_data.csv')
	all_performance = []

	try:
		files = glob.glob('/home/despoB/mb3152/scanner_performance_data/%s_tfMRI_*%s*_Stats.csv' %(subject,task))
		performance = []
		for f in files:
			df = pd.read_csv(f)
			if task == 'WM':
				t_performance = np.mean(df['Value'][[24,27,30,33]])
			if task == 'RELATIONAL':
				t_performance = np.mean([df['Value'][0],df['Value'][1]])
			if task == 'LANGUAGE':
				t_performance = np.mean([df['Value'][2],df['Value'][5]])
				s1 = bdf['ReadEng_AgeAdj'][bdf.Subject == int(subject)] 
				s2 = bdf['PicVocab_AgeAdj'][bdf.Subject == int(subject)]
				t_performance = np.nanmean([t_performance,s1,s2])
			if task == 'SOCIAL':
				t_performance = np.mean([df['Value'][0],df['Value'][5]])
			performance.append(t_performance)
		all_performance.append(np.mean(performance))
	except:
		all_performance.append(np.nan)
	return np.array(all_performance)

def organize_HCP():
	subj = np.loadtxt('/home/despoB/kaihwang/bin/FC_Scripts/hcp_subj')
	HCP_Conditions = ['LANGUAGE', 'WM', 'RELATIONAL']
	ROIs = ['Morel', 'Thalamus_WTA']

	df = pd.DataFrame()
	df['sub'] = subj.astype('int')
	
	for condition in HCP_Conditions:
		df[condition] = np.nan
		for s in df['sub']:
			try:
				df.loc[df['sub']==s, condition] = get_HCP_task_performance(str(s),condition)
			except:
				df.loc[df['sub']==s, condition] = np.nan	
				
		for roi, ind in zip(ROIs,INDEX):

			for key in ind:
				df[condition+'_'+key] = 0

			for s in df['sub']:
				try:
					pc = np.loadtxt(Data_path+str(s)+'_'+roi+'_'+condition+'_meanPC')
				except:
					pc = [np.nan]*500

				for key in ind:		
					df.loc[df['sub']==s, condition+'_'+key] = pc[ind[key]] 
	
	df.to_csv(path_to_data_folder+'HCP_behav.csv')



if __name__ == "__main__":

	#organize_TDSigEI()
	organize_TRSE()
	#organize_HCP()






