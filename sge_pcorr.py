import sys
from FuncParcel import *

# script to do partial corr

subject = sys.stdin.read().strip('\n')

ts_path = '/home/despoB/connectome-thalamus/NotBackedUp/TS/'

pcorr_path = '/home/despoB/connectome-thalamus/Partial_CorrMats/'

fn = ts_path + subject + 'Thalamus_indices_TS_000.netts'
subcortical_ts = np.loadtxt(fn)

ROIs = ['Craddock_300_cortical'] #'Craddock_300_cortical' 'Cortical_CI', 'Cortical_ROIs' 'Craddock_300_cortical'
for roi in ROIs:
	fn = ts_path + subject + roi + '_TS_000.netts' 
	cortical_ts = np.loadtxt(fn)

	pcorr_mat = pcorr_subcortico_cortical_connectivity(subcortical_ts, cortical_ts)

	fn = pcorr_path + subject + roi + '_pcorr_mat'
	#np.savetxt(fn, pcorr_mat, fmt='%.4f')
	save_object(pcorr_mat, fn)

