#!/bin/bash

#get MGH data.

Source='/home/despoB/harvard_gsp/preprocessed/pipeline_three_pipelines__scrub'

WD='/home/despoB/connectome-thalamus/MGH'


cd ${Source}
for s in $(ls -d Sub*); do

	if [ ! -d ${WD}/${s} ]; then
		mkdir ${WD}/${s}
	fi

	if [ ! -d ${WD}/${s}/run1 ]; then
		mkdir ${WD}/${s}/run1
	fi

	if [ ! -d ${WD}/${s}/run2 ]; then
		mkdir ${WD}/${s}/run2
	fi

	if [ ! -d ${WD}/${s}/run1/compcor_pipeline/ ]; then
	mkdir ${WD}/${s}/run1/compcor_pipeline/
	fi

	if [ ! -d ${WD}/${s}/run2/compcor_pipeline/ ]; then
	mkdir ${WD}/${s}/run2/compcor_pipeline/
	fi


	if [ ! -d ${WD}/${s}/run1/gbreg_pipeline/ ]; then
	mkdir ${WD}/${s}/run1/gbreg_pipeline/
	fi


	if [ ! -d ${WD}/${s}/run2/gbreg_pipeline/ ]; then
	mkdir ${WD}/${s}/run2/gbreg_pipeline/
	fi


	if [ ! -d ${WD}/${s}/run1/reg_pipeline/ ]; then
	mkdir ${WD}/${s}/run1/reg_pipeline/
	fi


	if [ ! -d ${WD}/${s}/run2/reg_pipeline/ ]; then
	mkdir ${WD}/${s}/run2/reg_pipeline/
	fi

	if [ ! -e ${WD}/${s}/run1/compcor_pipeline/scrubbed_preprocessed_antswarp.nii.gz ]; then
	ln -s /home/despoB/harvard_gsp/preprocessed/pipeline_three_pipelines__scrub/${s}/functional_mni/_scan_Scan_02_BOLD1/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_threshold_0.2/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global0.motion1.quadratic1.gm0.compcor1.csf0/_bandpass_freqs_0.009.0.08/scrubbed_preprocessed_antswarp.nii.gz ${WD}/${s}/run1/compcor_pipeline/scrubbed_preprocessed_antswarp.nii.gz
	fi

	if [ ! -e ${WD}/${s}/run2/compcor_pipeline/scrubbed_preprocessed_antswarp.nii.gz ]; then
	ln -s /home/despoB/harvard_gsp/preprocessed/pipeline_three_pipelines__scrub/${s}/functional_mni/_scan_Scan_03_BOLD2/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_threshold_0.2/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global0.motion1.quadratic1.gm0.compcor1.csf0/_bandpass_freqs_0.009.0.08/scrubbed_preprocessed_antswarp.nii.gz ${WD}/${s}/run2/compcor_pipeline/scrubbed_preprocessed_antswarp.nii.gz
	fi

	if [ ! -e ${WD}/${s}/run1/gbreg_pipeline/scrubbed_preprocessed_antswarp.nii.gz ]; then
	ln -s /home/despoB/harvard_gsp/preprocessed/pipeline_three_pipelines__scrub/${s}/functional_mni/_scan_Scan_02_BOLD1/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_threshold_0.2/_compcor_ncomponents_5_selector_pc10.linear1.wm1.global1.motion1.quadratic1.gm0.compcor0.csf1/_bandpass_freqs_0.009.0.08/scrubbed_preprocessed_antswarp.nii.gz ${WD}/${s}/run1/gbreg_pipeline/scrubbed_preprocessed_antswarp.nii.gz
	fi

	if [ ! -e ${WD}/${s}/run2/gbreg_pipeline/scrubbed_preprocessed_antswarp.nii.gz ]; then
	ln -s /home/despoB/harvard_gsp/preprocessed/pipeline_three_pipelines__scrub/${s}/functional_mni/_scan_Scan_03_BOLD2/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_threshold_0.2/_compcor_ncomponents_5_selector_pc10.linear1.wm1.global1.motion1.quadratic1.gm0.compcor0.csf1/_bandpass_freqs_0.009.0.08/scrubbed_preprocessed_antswarp.nii.gz ${WD}/${s}/run2/gbreg_pipeline/scrubbed_preprocessed_antswarp.nii.gz
	fi

	if [ ! -e ${WD}/${s}/run1/reg_pipeline/scrubbed_preprocessed_antswarp.nii.gz ]; then
	ln -s /home/despoB/harvard_gsp/preprocessed/pipeline_three_pipelines__scrub/${s}/functional_mni/_scan_Scan_02_BOLD1/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_threshold_0.2/_compcor_ncomponents_5_selector_pc10.linear1.wm1.global0.motion1.quadratic1.gm0.compcor0.csf1/_bandpass_freqs_0.009.0.08/scrubbed_preprocessed_antswarp.nii.gz ${WD}/${s}/run1/reg_pipeline/scrubbed_preprocessed_antswarp.nii.gz
	fi

	if [ ! -e ${WD}/${s}/run2/reg_pipeline/scrubbed_preprocessed_antswarp.nii.gz ]; then
	ln -s /home/despoB/harvard_gsp/preprocessed/pipeline_three_pipelines__scrub/${s}/functional_mni/_scan_Scan_03_BOLD2/_csf_threshold_0.96/_gm_threshold_0.7/_wm_threshold_0.96/_threshold_0.2/_compcor_ncomponents_5_selector_pc10.linear1.wm1.global0.motion1.quadratic1.gm0.compcor0.csf1/_bandpass_freqs_0.009.0.08/scrubbed_preprocessed_antswarp.nii.gz ${WD}/${s}/run2/reg_pipeline/scrubbed_preprocessed_antswarp.nii.gz
	fi

done

cd $WD
find -L */*/*/*.nii.gz -type l -delete