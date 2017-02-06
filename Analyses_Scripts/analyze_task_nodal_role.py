# get graph metrics for TRSE, TDSigEI, hcp datasets
from FuncParcel import *

# function
def cal_parcel_pcorr_graph(Data_path, subject, pipeline, condition, roi, Cortical_ROIs, Thalamus_parcel_positions, CI):
	fn = Data_path + str(subject) +'_' + pipeline +'_' + condition + '_' + roi + '_pcorr_mat' 
	adj = np.loadtxt(fn)
	adj[np.isnan(adj)]=0

	_, _, cost_thresholds = map_subcortical_cortical_targets(adj, Cortical_ROIs, Thalamus_parcel_positions)

	PCs = []
	for c in cost_thresholds:
		Par_adj = adj.copy()
		Par_adj[Par_adj<c]=0
		PCs += [bct.participation_coef(Par_adj, CI)]

	mean_PC = np.sum(PCs,axis=0)/13.5

	return mean_PC, PCs

def cal_parcel_fullcorr_graph(Data_path, subject, pipeline, condition, roi, Cortical_ROIs, Thalamus_parcel_positions, CI):
	fn = Data_path + str(subject) +'_' + pipeline +'_' + condition + '_' + roi + '_corrmat.txt'  
	adj = np.loadtxt(fn)
	adj[np.isnan(adj)]=0

	PCs = []
	for c in np.arange(0.01,0.16, 0.01):
		M = bct.threshold_proportional(adj, c, copy=True)
		PCs += [bct.participation_coef(M, CI)]

	mean_PC = np.sum(PCs,axis=0)/13.5

	return mean_PC, PCs	

# general variables	
ROIs = ['Gordon_plus_Morel_LPI', 'Gordon_plus_thalamus_WTA_LPI']
path_to_ROIs = '/home/despoB/connectome-thalamus/ROIs'
Cortical_CI = np.loadtxt(path_to_ROIs + '/Gordon_consensus_CI')
Cortical_ROIs = np.loadtxt(path_to_ROIs+'/Gordon_333', dtype = int)
Output_path = '/home/despoB/kaihwang/Rest/Graph/'
Data_path='/home/despoB/kaihwang/Rest/NotBackedUp/ParMatrices/'

#organiz for func parcel
Thalamus_Parcels = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
func_CI = np.append(Cortical_CI, Thalamus_Parcels)
#np.append(Cortical_CI, Thalamus_Parcels)
func_Thalamus_parcel_positions = np.arange(len(Cortical_CI),len(np.append(Cortical_CI, Thalamus_Parcels)),1)

#organize for 
Thalamus_Parcels = np.array([1]*16) #an extra empty "400" thalamic parcel was included... 
morel_CI = np.append(Cortical_CI, Thalamus_Parcels)
Morel_Thalamus_parcel_positions = np.arange(len(Cortical_CI),len(np.append(Cortical_CI, Thalamus_Parcels)),1)

# put list together
CIs = [morel_CI, func_CI]
Parcel_positions = [Morel_Thalamus_parcel_positions, func_Thalamus_parcel_positions]


def run_TRSE():
	## run TRSE
	subj = np.loadtxt('/home/despoB/kaihwang/bin/FC_Scripts/trse_subjects')
	Conditions = ['Cat', 'FH', 'HF', 'Both']
	pipeline = 'FIR'

	for subject in subj:
		subject = subject.astype("int")
		for condition in Conditions:
			for roi, CI, pos in zip (ROIs, CIs, Parcel_positions) :
				meanPC, _= cal_parcel_fullcorr_graph(Data_path, subject, pipeline, condition, roi, Cortical_ROIs, pos, CI)
				fn = Output_path + str(subject) + '_' + roi + '_' + condition + '_' +'meanPC'
				np.savetxt(fn, meanPC)

def run_TDSigEI():
	## run TDSigEI
	subj = np.loadtxt('/home/despoB/kaihwang/bin/FC_Scripts/TDSigEI_subj')
	Conditions = ['FH', 'HF', 'Fp', 'Hp']
	pipeline = 'nusiance'
	for subject in subj:
		subject = subject.astype("int")
		for condition in Conditions:
			for roi, CI, pos in zip (ROIs, CIs, Parcel_positions) :
				meanPC, _= cal_parcel_fullcorr_graph(Data_path,  subject, pipeline, condition, roi, Cortical_ROIs, pos, CI)
				fn = Output_path + str(subject) + '_' + roi + '_' + condition + '_' +'meanPC'
				np.savetxt(fn, meanPC)	

def run_HCP():
	## run HCP
	subj = np.loadtxt('/home/despoB/kaihwang/bin/FC_Scripts/hcp_subj')
	ROIs = ['Morel', 'Thalamus_WTA']
	
	for subject in subj:
		
		subject = subject.astype("int")

		Conditions = ['LANGUAGE', 'WM', 'RELATIONAL']
		phase = ['LR', 'RL']	

		#average left and right encoding matrices
		for condition in Conditions:
			try:
				for roi in ROIs:
					Left_plus_right_adj=[]
					for p in phase:		
						fn = Data_path + str(subject) +'_' + condition + '_' + p + '_' + roi + '_pcorr_mat'  #for TRSE and TDSig have to add FIR
						adj = np.loadtxt(fn)
						adj[np.isnan(adj)]=0
						Left_plus_right_adj += [adj]
					ave_adj = np.mean(Left_plus_right_adj, axis = 0)
					fn = Data_path + str(subject) +'_' + condition + '_' + roi + '_pcorr_mat'	
					np.savetxt(fn, ave_adj)
			except:
				print("file can't be found for %s" %subject)		

		for condition in Conditions:	
			try:		
				for roi, CI, pos in zip (ROIs, CIs, Parcel_positions) :
					meanPC, _= cal_parcel_graph(Data_path,  subject, condition, roi, Cortical_ROIs, pos, CI)	
					fn = Output_path + str(subject) + '_' + roi + '_' + condition + '_' +'meanPC'
					np.savetxt(fn, meanPC)	
			except:
				meanPC = [np.nan]*400
				fn = Output_path + str(subject) + '_' + roi + '_' + condition + '_' +'meanPC'
				np.savetxt(fn, meanPC)	
		print("subject %s done" %subject)		

if __name__ == "__main__":

	#run_TDSigEI()
	run_TRSE()
	#run_HCP()


