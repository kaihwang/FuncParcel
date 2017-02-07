## script to visualizee volume results in surcace, using pysurfer
import os
from surfer import Brain, project_volume_data
import numpy as np

### script to visualize cognitive components
subject_id = "fsaverage"
subjects_dir = os.environ["SUBJECTS_DIR"]

brain = Brain("fsaverage", "lh", "inflated", views=['lat'], background="white")
overlay_file = "/Volumes/neuro/Rest/ROIs/Yeo_12.nii.gz"
reg_file = os.path.join(os.environ["FREESURFER_HOME"], "average/mni152.register.dat")
zstat = project_volume_data(overlay_file, "lh", reg_file)
zstat = project_volume_data(overlay_file, "lh", subject_id="fsaverage", smooth_fwhm=4)
brain.add_data(zstat, min=1, max=5, thresh=1, colormap="hot", colorbar=False)

brain.show_view("medial")
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

brain = Brain("fsaverage", "lh", "inflated", views=['lat', 'med'], background="white")
CIs = np.loadtxt('/Volumes/neuro/Rest/ROIs/Gordon_consensus_CI')
coords = np.loadtxt('/Volumes/neuro/Rest/ROIs/Gordon_coordinates')
PC = np.loadtxt('//Volumes/neuro/Rest/Thalamic_parcel/roiPC')
PC = PC/100
PC[PC>0.8] = 1
#PC[PC<0.3] = 0
PC = PC/0.9

from matplotlib import cm
for i, coord in enumerate(coords):
	if CIs[i] != 0:
		if coord[0] < 0:
			brain.add_foci(coord, map_surface="white", color=cm.inferno(int(PC[i]*255))[0:3])

# save 1 image
brain.save_image("pc_lat.tiff")

###### to visualize cortical ROI wmd values
subject_id = "fsaverage"
subjects_dir = os.environ["SUBJECTS_DIR"]

brain = Brain("fsaverage", "lh", "inflated", views=['lat', 'med'], background="white")
CIs = np.loadtxt('/Volumes/neuro/Rest/ROIs/Gordon_consensus_CI')
coords = np.loadtxt('/Volumes/neuro/Rest/ROIs/Gordon_coordinates')
WMD = np.loadtxt('//Volumes/neuro/Rest/Thalamic_parcel/roiWMD')
WMD = WMD/100
WMD = WMD+1.5
WMD[WMD<0] =0
WMD = WMD/3
WMD[WMD>1] = 1

from matplotlib import cm
for i, coord in enumerate(coords):
	if CIs[i] != 0:
		if coord[0] < 0:
			#if wmd[i] > 80:
			brain.add_foci(coord, map_surface="white", color=cm.bwr(int(WMD[i]*255))[0:3])
brain.save_image("WMD_lat.tiff")


### to visualize hubs
subject_id = "fsaverage"
subjects_dir = os.environ["SUBJECTS_DIR"]

brain = Brain("fsaverage", "rh", "inflated", views=['lat', 'med'], background="white")
coords = np.loadtxt('/Volumes/neuro/Rest/ROIs/Gordon_coordinates')
PC = np.loadtxt('//Volumes/neuro/Rest/Thalamic_parcel/roipc')
wmd = np.loadtxt('//Volumes/neuro/Rest/Thalamic_parcel/wmdroi')

for i, coord in enumerate(coords):
	if coord[0] > 0:
		if wmd[i] > 80:
			brain.add_foci(coord, map_surface="white", color='indigo')
		if 	PC[i] > 61:
			brain.add_foci(coord, map_surface="white", color='darkgreen')
		if (wmd[i] > 61) & (PC[i] > 80):
			brain.add_foci(coord, map_surface="white", color='gold')	


#### to visualize AN, MD, PuM, Intra ROI target

subject_id = "fsaverage"
subjects_dir = os.environ["SUBJECTS_DIR"]

brain = Brain("fsaverage", "rh", "inflated", views=['lat', 'med'], background="white")
#brain = Brain("fsaverage", "rh", "inflated", background="white")
coords = np.loadtxt('/Volumes/neuro/Rest/ROIs/Gordon_coordinates')
CIs = np.loadtxt('/Volumes/neuro/Rest/ROIs/Gordon_consensus_CI')
targets = np.loadtxt('/Volumes/neuro/Rest/ROIs/Morel_VL_targets')

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



