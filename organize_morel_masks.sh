
cd /home/despoB/kaihwang/Rest/ROIs/MorelAtlasMNI152

cd right-vols-1mm/

rm temp*
3dcalc -a AD.nii.gz -b AM.nii.gz \
-expr "(a*1 + b*2) * equals((a/a)+(b/b),1) + b*2*(astep(a+b,1))" \
-prefix temp2.nii.gz 

i=2
for nu in AV LD MDmc MDpc MV CL CeM CM Pf PuM PuI PuL PuA LP MGN SG Li Po LGNmc LGNpc VPLa VPLp VPM VPI VLa VLpd VLpv VAmc VApc VM; do
	ii=$((i+1))
	3dcalc -a temp${i}.nii.gz -b ${nu}.nii.gz \
	-expr "(a*1 + b*${ii}) * equals((a/a)+(b/b),1) + b*${ii}*(astep((a/a)+(b/b),1))" \
	-prefix temp${ii}.nii.gz
	i=$((i+1))

done

cd ../left-vols-1mm/

rm temp*
3dcalc -a AD.nii.gz -b AM.nii.gz \
-expr "(a*1 + b*2) * equals((a/a)+(b/b),1) + b*2*(astep(a+b,1))" \
-prefix temp2.nii.gz 

i=2
for nu in AV LD MDmc MDpc MV CL CeM CM Pf PuM PuI PuL PuA LP MGN SG Li Po LGNmc LGNpc VPLa VPLp VPM VPI VLa VLpd VLpv VAmc VApc VM; do
	ii=$((i+1))
	3dcalc -a temp${i}.nii.gz -b ${nu}.nii.gz \
	-expr "(a*1 + b*${ii}) * equals((a/a)+(b/b),1) + b*${ii}*(astep((a/a)+(b/b),1))" \
	-prefix temp${ii}.nii.gz
	i=$((i+1))

done


cd /home/despoB/kaihwang/Rest/ROIs/MorelAtlasMNI152
3dcalc -a /home/despoB/kaihwang/Rest/ROIs/MorelAtlasMNI152/right-vols-1mm/temp32.nii.gz \
-b /home/despoB/kaihwang/Rest/ROIs/MorelAtlasMNI152/left-vols-1mm/temp32.nii.gz \
-expr 'a+b' -prefix morel_overlap_mask.nii.gz