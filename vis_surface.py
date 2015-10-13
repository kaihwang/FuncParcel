## script to visualizee volume results in surcace, using pysurfer


import os
from surfer import Brain, project_volume_data
import numpy as np

#load data to visualize
#mask_file = "Cortical_CI.nii.gz"


mask_file = "MGH_cortical_nodal_role_wmd.nii.gz"

reg_file = os.path.join(os.environ["FREESURFER_HOME"], "average/mni152.register.dat")


#cmap = ['red', 'blue', 'cyan', 'yellow', 'teal', 'purple', 'pink', 'green', 'white', 'white', 'magenta']
#cmap = ['Navy','Crimson','Khaki','lightblue','Lime']
cmap = 'spring'
#brain = Brain("fsaverage", "rh", "inflated", cortex='bone')
subject_id = "fsaverage"
hemi = "rh"
surface = "pial"

brain = Brain(subject_id, hemi, surface)
#mask = project_volume_data(mask_file, "rh", reg_file, projmeth ='frac', smooth_fwhm=0, projarg=[0,1,0.1], projsum="max")
mask = project_volume_data(mask_file, "rh", reg_file)

#mask = project_volume_data(mask_file, "rh", reg_file, smooth_fwhm=0, projsum="max").astype(int)
#brain.add_data(mask, min=-100, max=300, thresh=-200, colormap=cmap, alpha=1, colorbar=True, smoothing_steps = [], remove_existing = True)
