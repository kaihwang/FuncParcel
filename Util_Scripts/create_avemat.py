import FuncParcel
import numpy as np
#Script for generating ave matrices

#thalamocortical
#MGH_Gordon_333_thalamocortical_pcorr_avemat = FuncParcel.average_corrmat('/home/despoB/connectome-thalamus/Partial_CorrMats/MGH*Gordon*cortical*')
#np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/MGH_Gordon_333_thalamocortical_pcorr_avemat', MGH_Gordon_333_thalamocortical_pcorr_avemat)

#NKI_645_Gordon_333_thalamocortical_pcorr_avemat = FuncParcel.average_corrmat('/home/despoB/connectome-thalamus/Partial_CorrMats/NKI*645*Gordon*cortical*')
#np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/NKI_645_Gordon_333_thalamocortical_pcorr_avemat', NKI_645_Gordon_333_thalamocortical_pcorr_avemat)

#NKI_1400_Gordon_333_thalamocortical_pcorr_avemat = FuncParcel.average_corrmat('/home/despoB/connectome-thalamus/Partial_CorrMats/NKI*1400*Gordon*cortical*')
#np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/NKI_1400_Gordon_333_thalamocortical_pcorr_avemat', NKI_1400_Gordon_333_thalamocortical_pcorr_avemat)

#corticocortical
#MGH_Gordon_333_cortical_avemat = FuncParcel.average_corrmat('/home/despoB/kaihwang/Rest/NotBackedUp/AdjMatrices/MGH*Gordon*cortical*', np_txt=True, pickle_object=False)
#np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/MGH_Gordon_333_cortical_avemat', MGH_Gordon_333_cortical_avemat)

#NKI_645_Gordon_333_cortical_avemat = FuncParcel.average_corrmat('/home/despoB/kaihwang/Rest/NotBackedUp/AdjMatrices/NKI*Gordon*cortical*645*', np_txt=True, pickle_object=False)
#np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/NKI_645_Gordon_333_cortical_avemat', NKI_645_Gordon_333_cortical_avemat)

#NKI_1400_Gordon_333_cortical_avemat = FuncParcel.average_corrmat('/home/despoB/kaihwang/Rest/NotBackedUp/AdjMatrices/NKI*Gordon*cortical*1400*', np_txt=True, pickle_object=False)
#np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/NKI_1400_Gordon_333_cortical_avemat', NKI_1400_Gordon_333_cortical_avemat)

# MGH_Gordon_ConsensusCI_thalamocortical_avemat = FuncParcel.average_corrmat('/home/despoB/connectome-thalamus/Partial_CorrMats/MGH*Gordon*consensus_CI*')
# np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/MGH_Gordon_ConsesnsusCI_thalamocortical_pcorr_avemat', MGH_Gordon_ConsensusCI_thalamocortical_avemat)

# NKI_645_Gordon_ConsensusCI_thalamocortical_avemat = FuncParcel.average_corrmat('/home/despoB/connectome-thalamus/Partial_CorrMats/NKI*645*Gordon*consensus_CI*')
# np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/NKI_645_Gordon_ConsesnsusCI_thalamocortical_pcorr_avemat', NKI_645_Gordon_ConsensusCI_thalamocortical_avemat)

# NKI_1400_Gordon_ConsensusCI_thalamocortical_avemat = FuncParcel.average_corrmat('/home/despoB/connectome-thalamus/Partial_CorrMats/NKI*1400*Gordon*consensus_CI*')
# np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/NKI_1400_Gordon_ConsesnsusCI_thalamocortical_pcorr_avemat', NKI_1400_Gordon_ConsensusCI_thalamocortical_avemat)


#thalamocortical and corticocortical using Gordon and WTA parcels
#MGH_Gordon_ThaWTA_avemat = FuncParcel.average_corrmat('/home/despoB/kaihwang/Rest/NotBackedUp/ParMatrices/MGH*_Gordon_333_cortical_ThaWTA_pcorr_mat')
#np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/MGH_Gordon_ThaWTA_avemat', MGH_Gordon_ThaWTA_avemat)

MGH_Gordon_333_cortical_avemat = FuncParcel.average_corrmat('/home/despoB/kaihwang/Rest/NotBackedUp/ParMatrices/MGH*_Gordon_333_cortical_corrmat', np_txt=True, pickle_object=False)
np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/MGH_Gordon_333_cortical_avemat', MGH_Gordon_333_cortical_avemat)

MGH_Gordon_333_thalamocortical_pcorr_avemat = FuncParcel.average_corrmat('/home/despoB/kaihwang/Rest/NotBackedUp/ParMatrices/Partial_CorrMats/MGH*Gordon*cortical*')
np.savetxt('/home/despoB/kaihwang/Rest/Thalamic_parcel/MGH_Gordon_333_thalamocortical_pcorr_avemat', MGH_Gordon_333_thalamocortical_pcorr_avemat)
