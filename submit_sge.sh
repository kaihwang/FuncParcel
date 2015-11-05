#!/bin/bash

# bash script to use sge to submit jobs iterating modularity detection

# echo 'python /home/despoB/kaihwang/bin/FuncParcel/create_avemat.py' > ~/tmp/avg_partition.sh
# qsub -V -M kaihwang -m e -e ~/tmp -o ~/tmp ~/tmp/avg_partition.sh

# for Subject in 1103 1220 1306 1223 1314 1311 1318 1313 1326 1325 1328 1329 1333 1331 1335 1338 1336 1339 1337 1344 1340 128 162 163 168 176 b116 b117 b120 b121 b122 b138 b143 b153; do

# 	sed "s/128/${Subject}/g" < python_brainx_mod.sh > ~/tmp/python_brainx_${Subject}.sh
# 	qsub -V -M kaihwang -m e -e ~/tmp -o ~/tmp ~/tmp/python_brainx_${Subject}.sh

# done


WD='/home/despoB/connectome-thalamus/NKI'
SCRIPT='/home/despoB/kaihwang/bin/FuncParcel'
cd ${WD}

for s in $(ls -d 0*); do

	for seqq in __mx_1400 __mx_645; do #__mx_1400 __mx_645
		sed "s/128/NKI_${s}${seqq}_/g" < ${SCRIPT}/python_brainx_mod.sh > ~/tmp/pcorr${s}${seqq}.sh
		qsub -V -M kaihwang -m e -e ~/tmp -o ~/tmp ~/tmp/pcorr${s}${seqq}.sh

	done

done


WD='/home/despoB/connectome-thalamus/MGH'
SCRIPT='/home/despoB/kaihwang/bin/FuncParcel'
cd ${WD}

for s in $(ls -d Sub*); do
	if [ -e ${WD}/${s}/MNINonLinear/rfMRI_REST_ncsreg.nii.gz ]; then
		sed "s/128/MGH_${s}_/g" < ${SCRIPT}/python_brainx_mod.sh > ~/tmp/pcorr${s}.sh
		qsub -V -M kaihwang -m e -e ~/tmp -o ~/tmp ~/tmp/pcorr${s}.sh
	fi
done
