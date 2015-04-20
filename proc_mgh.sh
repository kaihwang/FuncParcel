#!/bin/bash

#get MGH data.

Source='/home/despoB/harvard_gsp/preprocessed/pipeline_three_pipelines__scrub'

WD='/home/despoB/connectome-thalamus/MGH'



for s in $(ls ${Source}); do

	if [ ! -d ${WD}/${s} ]; then
		mkdir ${WD}/${s}
	fi

	if [ ! -d ${WD}/${s}/run1 ]; then
		mkdir ${WD}/${s}/run1
	fi

	if [ ! -d ${WD}/${s}/run2 ]; then
		mkdir ${WD}/${s}/run2
	fi

	
	mkdir ${WD}/${s}/run1/compcor_pipeline/
	mkdir ${WD}/${s}/run2/compcor_pipeline/
	mkdir ${WD}/${s}/run1/gbreg_pipeline/
	mkdir ${WD}/${s}/run2/gbreg_pipeline/

	ln -s /home/despoB/harvard_gsp/preprocessed/pipeline_three_pipelines__scrub/${s}/functional_mni/_scan_Scan_02_BOLD1/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_threshold_0.2/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global0.motion1.quadratic1.gm0.compcor1.csf0/_bandpass_freqs_0.009.0.08/scrubbed_preprocessed_antswarp.nii.gz ${WD}/${s}/run1/compcor_pipeline/scrubbed_preprocessed_antswarp.nii.gz
	ln -s /home/despoB/harvard_gsp/preprocessed/pipeline_three_pipelines__scrub/${s}/functional_mni/_scan_Scan_03_BOLD2/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_threshold_0.2/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global0.motion1.quadratic1.gm0.compcor1.csf0/_bandpass_freqs_0.009.0.08/scrubbed_preprocessed_antswarp.nii.gz ${WD}/${s}/run2/compcor_pipeline/scrubbed_preprocessed_antswarp.nii.gz
	ln -s /home/despoB/harvard_gsp/preprocessed/pipeline_three_pipelines__scrub/${s}/functional_mni/_scan_Scan_02_BOLD1/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_threshold_0.2/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global1.motion1.quadratic1.gm0.compcor0.csf1/_bandpass_freqs_0.009.0.08/scrubbed_preprocessed_antswarp.nii.gz ${WD}/${s}/run1/gbreg_pipeline/scrubbed_preprocessed_antswarp.nii.gz
	ln -s /home/despoB/harvard_gsp/preprocessed/pipeline_three_pipelines__scrub/${s}/functional_mni/_scan_Scan_03_BOLD2/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_threshold_0.2/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global1.motion1.quadratic1.gm0.compcor0.csf1/_bandpass_freqs_0.009.0.08/scrubbed_preprocessed_antswarp.nii.gz ${WD}/${s}/run2/gbreg_pipeline/scrubbed_preprocessed_antswarp.nii.gz


done