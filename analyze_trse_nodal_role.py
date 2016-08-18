# get graph metrics for TRSE, TDSigEI, hcp datasets
from FuncParcel import *

# function
def cal_parcel_graph(Data_path, subject,condition, roi, Cortical_ROIs, Thalamus_parcel_positions, CI):
	fn = Data_path + str(subject) +'_' + condition + '_' + roi + '_pcorr_mat'  #for TRSE and TDSig have to add FIR
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

# general variables	
ROIs = ['Morel_LPI', 'Thalamus_WTA_LPI']
path_to_ROIs = '/home/despoB/connectome-thalamus/ROIs'
Cortical_CI = np.loadtxt(path_to_ROIs + '/Gordon_consensus_CI')
Cortical_ROIs = np.loadtxt(path_to_ROIs+'/Gordon_333', dtype = int)

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
	TRSE_Data_path='/home/despoB/kaihwang/Rest/NotBackedUp/ParMatrices/'
	subj = np.loadtxt('/home/despoB/kaihwang/bin/FC_Scripts/trse_subjects')
	Conditions = ['Cat', 'FH', 'HF', 'Both']


	for subject in subj:
		subject = subject.astype("int")
		for condition in Conditions:
			for roi, CI, pos in zip (ROIs, CIs, Parcel_positions) :
				meanPC, _= cal_parcel_graph(TRSE_Data_path, subject, condition, roi, Cortical_ROIs, pos, CI)
				fn = TRSE_Data_path + str(subject) + '_' + roi + '_' + condition + '_' +'meanPC'
				np.savetxt(fn, meanPC)

def run_TDSigEI():
	## run TRSE
	TD_Data_path='/home/despoB/kaihwang/Rest/NotBackedUp/ParMatrices/'
	subj = np.loadtxt('/home/despoB/kaihwang/bin/FC_Scripts/TDSigEI_subj')
	Conditions = ['FH', 'HF', 'Fp', 'Hp']

	for subject in subj:
		subject = subject.astype("int")
		for condition in Conditions:
			for roi, CI, pos in zip (ROIs, CIs, Parcel_positions) :
				meanPC, _= cal_parcel_graph(TD_Data_path,  subject, condition, roi, Cortical_ROIs, pos, CI)
				fn = TD_Data_path + str(subject) + '_' + roi + '_' + condition + '_' +'meanPC'
				np.savetxt(fn, meanPC)	

def run_HCP():
	## run TRSE
	TD_Data_path='/home/despoB/kaihwang/Rest/NotBackedUp/ParMatrices/'
	subj = np.loadtxt('/home/despoB/kaihwang/bin/FC_Scripts/hcp_subj')
	Conditions = ['LANGUAGE_LR', 'LANGUAGE_RL', 'WM_LR', 'WM_RL', 'RELATIONAL_LR', 'RELATIONAL_RL']

	for subject in subj:
		subject = subject.astype("int")
		for condition in Conditions:
			for roi, CI, pos in zip (ROIs, CIs, Parcel_positions) :
				meanPC, _= cal_parcel_graph(TD_Data_path,  subject, condition, roi, Cortical_ROIs, pos, CI)
				fn = TD_Data_path + str(subject) + '_' + roi + '_' + condition + '_' +'meanPC'
				np.savetxt(fn, meanPC)								

if __name__ == "__main__":

	#run_TDSigEI()
	#run_TRSE()


