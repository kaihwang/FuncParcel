import FuncParcel
import numpy as np

#thalamo Cortical
Old_ThalamoCorticalAveMat = FuncParcel.average_corrmat('/home/despoB/connectome-thalamus/Partial_CorrMats/MGH*Craddock*')
np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/MGH_Craddock_300_cortical_plus_thalamus_parcorrmatavg', Old_ThalamoCorticalAveMat)

Old_ThalamoCorticalAveMat = FuncParcel.average_corrmat('/home/despoB/connectome-thalamus/Partial_CorrMats/NKI*645*Craddock*')
np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/NKI_mx_645_Craddock_300_cortical_plus_thalamus_parcorrmatavg', Old_ThalamoCorticalAveMat)

Old_ThalamoCorticalAveMat = FuncParcel.average_corrmat('/home/despoB/connectome-thalamus/Partial_CorrMats/NKI*1400*Craddock*')
np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/NKI_mx_1400_Craddock_300_cortical_plus_thalamus_parcorrmatavg', Old_ThalamoCorticalAveMat)


#thalamo cortex
# Old_ThalamoCorticalAveMat = FuncParcel.average_corrmat('/home/despoB/connectome-thalamus/Partial_CorrMats/MGH*ROI*')
# np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/MGH_cortex_plus_thalamus_parcorrmatavg', Old_ThalamoCorticalAveMat)

# Old_ThalamoCorticalAveMat = FuncParcel.average_corrmat('/home/despoB/connectome-thalamus/Partial_CorrMats/NKI*645*ROI*')
# np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/NKI_mx_645_cortex_plus_thalamus_parcorrmatavg', Old_ThalamoCorticalAveMat)

# Old_ThalamoCorticalAveMat = FuncParcel.average_corrmat('/home/despoB/connectome-thalamus/Partial_CorrMats/NKI*1400*ROI*')
# np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/NKI_mx_1400_cortex_plus_thalamus_parcorrmatavg', Old_ThalamoCorticalAveMat)


#Old_ThalamoCorticalAveMat = FuncParcel.average_corrmat('//home/despoB/connectome-thalamus/NotBackedUp/ParMatrices/*NKI*CI*2500*')
#np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/NKI_std_2500_cortical_network_plus_thalamus_corrmatavg', Old_ThalamoCorticalAveMat)

#Old_CorticalAveMat = FuncParcel.average_corrmat('/home/despoB/kaihwang/Rest/AdjMatrx`ices/t1*FIX*')
#striatal
# Old_StriatalCorticalAveMat = FuncParcel.average_corrmat('/home/despoB/kaihwang/Rest/AdjMatrices/t1*FIX_striatalcortical*')
# np.savetxt('/home/despoB/kaihwang/Rest/Striatum_parcel/Old_StriatalCorticalAveMat', Old_StriatalCorticalAveMat)

# Control_Subj = ['1103', '1220', '1306', '1223', '1314', '1311', '1318', '1313', '1326', '1325', '1328', '1329', '1333', '1331', '1335', '1338', '1336', '1339', '1337', '1344']

# M_Sum =[];
# for s in Control_Subj:
# 	fn = '/home/despoB/kaihwang/Rest/AdjMatrices/t%s_Full_WashU333_corrmat' %s
# 	M = np.loadtxt(f)
# 	M_Sum += [M]

# Old_CorticalAveMat = sum(M_Sum)/len(Control_Subj)
# np.savetxt('/home/despoB/kaihwang/bin/FuncParcel/Data/Old_CorticalAveMat', Old_ThalamoCorticalAveMat)

