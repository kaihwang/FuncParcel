import FuncParcel
import numpy as np

#thalamo
Old_ThalamoCorticalAveMat = FuncParcel.average_corrmat('/home/despoB/kaihwang/Rest/AdjMatrices/t1*FIX_thalamocortical*')
np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/Old_ThalamoCorticalAveMat', Old_ThalamoCorticalAveMat)
#Old_CorticalAveMat = FuncParcel.average_corrmat('/home/despoB/kaihwang/Rest/AdjMatrices/t1*FIX*')
#striatal
Old_StriatalCorticalAveMat = FuncParcel.average_corrmat('/home/despoB/kaihwang/Rest/AdjMatrices/t1*FIX_striatalcortical*')
np.savetxt('/home/despoB/kaihwang/Rest/Striatum_parcel/Old_StriatalCorticalAveMat', Old_StriatalCorticalAveMat)

Control_Subj = ['1103', '1220', '1306', '1223', '1314', '1311', '1318', '1313', '1326', '1325', '1328', '1329', '1333', '1331', '1335', '1338', '1336', '1339', '1337', '1344']

M_Sum =[];
for s in Control_Subj:
	fn = '/home/despoB/kaihwang/Rest/AdjMatrices/t%s_Full_WashU333_corrmat' %s
	M = np.loadtxt(f)
	M_Sum += [M]

Old_CorticalAveMat = sum(M_Sum)/len(Control_Subj)
np.savetxt('/home/despoB/kaihwang/bin/FuncParcel/Data/Old_CorticalAveMat', Old_ThalamoCorticalAveMat)

