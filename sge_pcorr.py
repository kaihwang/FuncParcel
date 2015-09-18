import sys
from FuncParcel import *
from functools import partial
from multiprocessing import Pool
from itertools import product

pool = Pool(16)

# script to do partial corr

subject = sys.stdin.read().strip('\n')

ts_path = '/home/despoB/connectome-thalamus/NotBackedUp/TS/'

pcorr_path = '/home/despoB/connectome-thalamus/Partial_CorrMats/'

st = ts_path + subject + 'Thalamus_indices_TS_000.netts'
#st = np.loadtxt(fn)

ROIs = ['Craddock_300_cortical'] #'Craddock_300_cortical' 'Cortical_CI', 'Cortical_ROIs' 'Craddock_300_cortical'
for roi in ROIs:
	ct = ts_path + subject + roi + '_TS_000.netts' 
	#ct = np.loadtxt(fn)

	#create output
	pcorr_mat = np.zeros((np.loadtxt(st).shape[0],np.loadtxt(ct).shape[0]),dtype=np.float)

	#start parfor
	pf = partial(par_pcorr_subcortico_cortical_connectivity, subcortical_ts=st, cortical_ts=ct)
	chunksize = 1
	for ind, res in enumerate(pool.imap(pf, product(range(np.loadtxt(st).shape[0]), range(np.loadtxt(ct).shape[0])), chunksize=160)):
		pcorr_mat.flat[ind] = res

	#pcorr_mat = pcorr_subcortico_cortical_connectivity(subcortical_ts, cortical_ts)

	fn = pcorr_path + subject + roi + '_pcorr_mat'
	#np.savetxt(fn, pcorr_mat, fmt='%.4f')
	save_object(pcorr_mat, fn)









