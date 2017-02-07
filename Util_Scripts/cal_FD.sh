#!/bin/bash
#script to calculate FD



cd /home/despoB/kaihwang/Rest/MGH
Controls=$(/bin/ls -d *)

for s in $Controls; do
	
	cd /home/despoB/kaihwang/Rest/MGH/${s}/MNINonLinear
	cat rfMRI_REST1_motion_24.1D rfMRI_REST2_motion_24.1D > motion.1D
	1d_tool.py -infile motion.1D -derivative -write motion.diff.1D
	1deval -a motion.diff.1D[0] -expr '2*pi*50*(a/360)' > rot_x.1D
	1deval -a motion.diff.1D[1] -expr '2*pi*50*(a/360)' > rot_y.1D
	1deval -a motion.diff.1D[2] -expr '2*pi*50*(a/360)' > rot_z.1D
	1dcat rot_x.1D rot_y.1D rot_z.1D motion.diff.1D[3..5] > ${s}.motion.diff.despike.1D
	1deval -a ${s}.motion.diff.despike.1D[0] \
	-b ${s}.motion.diff.despike.1D[1] \
	-c ${s}.motion.diff.despike.1D[2] \
	-d ${s}.motion.diff.despike.1D[3] \
	-e ${s}.motion.diff.despike.1D[4] \
	-f ${s}.motion.diff.despike.1D[5] \
	-expr 'abs(a)+abs(b)+abs(c)+abs(d)+abs(e)+abs(f)' > ${s}_FD.1D

	rm motion.diff.1D
	rm rot_x.1D
	rm rot_y.1D
	rm rot_z.1D
	
	1deval -a ${s}_FD.1D -expr 'step(isnegative(a-0.5))' | grep -n 1 | cut -f1 -d: > tmp.1D
	1deval -a tmp.1D -expr 'a-1' > ${s}_scrub.1D
	wc -l ${s}_scrub.1D > scrub_vol

	rm tmp.1D
done


for s in $(cat /home/despoB/kaihwang/Rest/MGH_subjs); do cat /home/despoB/kaihwang/Rest/MGH/${s}/MNINonLinear/scrub_vol >> ~/scrub.csv; done