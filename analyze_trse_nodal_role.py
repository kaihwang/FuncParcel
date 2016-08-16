# get graph metrics from TRSE datasets
from FuncParcel import *


def cal_parcel_graph(subject,condition,roi,CI):

	Data_path='/home/despoB/kaihwang/TRSE/TRSEPPI/Graph/'
	fn = Data_path + str(subject) +'_' + condition + '_' + roi + '_LPI_bcorrmat.txt' 
	adj = np.loadtxt(fn)
	adj[np.isnan(adj)]=0

	PCs = []
	for ix, c in enumerate(np.arange(0.01,0.16, 0.01)):
		M = bct.threshold_proportional(adj, c, copy=True)
		PCs += [bct.participation_coef(M, CI)]

	mean_PC = np.sum(PCs,axis=0)/13.5

	return mean_PC, PCs

path_to_ROIs = '/home/despoB/connectome-thalamus/ROIs'
Cortical_CI = np.loadtxt(path_to_ROIs + '/Gordon_consensus_CI')
Thalamus_Parcels = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
func_CI = np.append(Cortical_CI, Thalamus_Parcels)
Thalamus_Parcels = np.array([1]*16) #an extra empty "400" thalamic parcel was included... 
morel_CI = np.append(Cortical_CI, Thalamus_Parcels)
CIs = [morel_CI, func_CI]

subj = np.loadtxt('/home/despoB/kaihwang/bin/FC_Scripts/trse_subjects')
Conditions = ['categ_face', 'categ_scene', 'relev_face', 'relev_scene', 'irrel_face', 'irrel_scene']
ROIs = ['Gordon_plus_Morel', 'Gordon_plus_thalamus_WTA']
Data_path='/home/despoB/kaihwang/TRSE/TRSEPPI/Graph/'

for subject in subj:
	subject = subject.astype("int")
	for condition in Conditions:
		for roi, CI in zip (ROIs, CIs) :
			meanPC, _= cal_parcel_graph(subject,condition,roi,CI)
			fn = Data_path + str(subject) + '_' + roi + '_' + condition + '_' +'meanPC'
			np.savetxt(fn, meanPC)


