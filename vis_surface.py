## script to visualizee volume results in surcace, using pysurfer
import os
from surfer import Brain, project_volume_data
import numpy as np

#load data to visualize
#mask_file = "Cortical_CI.nii.gz"

### to visualize ROI's CI assignment
#load coordinates
subject_id = "fsaverage"
subjects_dir = os.environ["SUBJECTS_DIR"]

brain = Brain("fsaverage", "rh", "inflated", views=['lat', 'med'], background="white")
#brain = Brain("fsaverage", "rh", "inflated", background="white")
coords = np.loadtxt('/Volumes/neuro/Rest/ROIs/Gordon_coordinates')
CIs = np.loadtxt('/Volumes/neuro/Rest/ROIs/Gordon_consensus_CI')

colors = {1:'red', 2:'purple', 3:'green', 4:'yellow', 5:'cyan', 6:'blue', 7:'brown', 8:'pink', 9:'teal', 12:'pink'}

for i, coord in enumerate(coords):
	if CIs[i] != 0:
		if coord[0] > 0:
			brain.add_foci(coord, map_surface="white", color=colors[CIs[i]])

brain2 = Brain("fsaverage", "lh", "inflated", views=['lat', 'med'], background="white")
colors = {1:'red', 2:'purple', 3:'green', 4:'yellow', 5:'cyan', 6:'blue', 7:'brown', 8:'pink', 9:'teal', 12:'pink'}

for i, coord in enumerate(coords):
	if CIs[i] != 0:
		if coord[0] < 0:
			brain2.add_foci(coord, map_surface="white", color=colors[CIs[i]])



### to visualize cortical ROI PC values
subject_id = "fsaverage"
subjects_dir = os.environ["SUBJECTS_DIR"]

brain = Brain("fsaverage", "rh", "inflated", views=['lat', 'med'], background="white")
#brain = Brain("fsaverage", "rh", "inflated", background="white")
coords = np.loadtxt('/Volumes/neuro/Rest/ROIs/Gordon_coordinates')
PC = np.loadtxt('//Volumes/neuro/Rest/Thalamic_parcel/roipc')

from matplotlib import cm
for i, coord in enumerate(coords):
	if coord[0] > 0:
		brain.add_foci(coord, map_surface="white", color=cm.spring(int(PC[i]/100*255))[0:3])


###### to visualize cortical ROI wmd values
subject_id = "fsaverage"
subjects_dir = os.environ["SUBJECTS_DIR"]

brain = Brain("fsaverage", "rh", "inflated", views=['lat', 'med'], background="white")
#brain = Brain("fsaverage", "rh", "inflated", background="white")
coords = np.loadtxt('/Volumes/neuro/Rest/ROIs/Gordon_coordinates')
wmd = np.loadtxt('//Volumes/neuro/Rest/Thalamic_parcel/wmdroi')

from matplotlib import cm
for i, coord in enumerate(coords):
	if coord[0] > 0:
		brain.add_foci(coord, map_surface="white", color=cm.spring(int((wmd[i]/100+2.1)/4*255))[0:3])



#### to visualize AN, MD, PuM, Intra ROI target

subject_id = "fsaverage"
subjects_dir = os.environ["SUBJECTS_DIR"]

brain = Brain("fsaverage", "rh", "inflated", views=['lat', 'med'], background="white")
#brain = Brain("fsaverage", "rh", "inflated", background="white")
coords = np.loadtxt('/Volumes/neuro/Rest/ROIs/Gordon_coordinates')
CIs = np.loadtxt('/Volumes/neuro/Rest/ROIs/Gordon_consensus_CI')
targets = np.loadtxt('/Volumes/neuro/Rest/ROIs/Morel_PuM_targets')

colors = {1:'red', 2:'purple', 3:'green', 4:'yellow', 5:'cyan', 6:'blue', 7:'brown', 8:'pink', 9:'teal', 12:'pink'}

for i, coord in enumerate(coords):
	if CIs[i] != 0:
		if np.in1d(i+1, targets):
			if coord[0] > 0:
				brain.add_foci(coord, map_surface="white", color=colors[CIs[i]])

brain2 = Brain("fsaverage", "lh", "inflated", views=['lat', 'med'], background="white")
colors = {1:'red', 2:'purple', 3:'green', 4:'yellow', 5:'cyan', 6:'blue', 7:'brown', 8:'pink', 9:'teal', 12:'pink'}

for i, coord in enumerate(coords):
	if CIs[i] != 0:
		if np.in1d(i+1, targets):
			if coord[0] < 0:
				brain2.add_foci(coord, map_surface="white", color=colors[CIs[i]])


# 0 
# 1 DF Red
# 2 CO Purple
# 3 SM Green
# 4 FP Yellow
# 5 Occipital-Parietal Cyan
# 6 Visual Blue
# 7 Retronspinal brown
# 8 Superiro temporal Pink
# 9 Attn Teal
# 10 Inferior temporal (delete)
# 11 OFC (delete)
# 12 middle temporal Pink

### to visualize volume results
# mask_file = "MGH_cortical_nodal_role_wmd.nii.gz"

# reg_file = os.path.join(os.environ["FREESURFER_HOME"], "average/mni152.register.dat")


# #cmap = ['red', 'blue', 'cyan', 'yellow', 'teal', 'purple', 'pink', 'green', 'white', 'white', 'magenta']
# #cmap = ['Navy','Crimson','Khaki','lightblue','Lime']
# cmap = 'spring'
# #brain = Brain("fsaverage", "rh", "inflated", cortex='bone')
# subject_id = "fsaverage"
# hemi = "rh"
# surface = "pial"

# brain = Brain(subject_id, hemi, surface)
# #mask = project_volume_data(mask_file, "rh", reg_file, projmeth ='frac', smooth_fwhm=0, projarg=[0,1,0.1], projsum="max")
# mask = project_volume_data(mask_file, "rh", reg_file)

# #mask = project_volume_data(mask_file, "rh", reg_file, smooth_fwhm=0, projsum="max").astype(int)
# #brain.add_data(mask, min=-100, max=300, thresh=-200, colormap=cmap, alpha=1, colorbar=True, smoothing_steps = [], remove_existing = True)
