# to consolidate Morel atlas partitions

cd /home/despoB/kaihwang/Rest/ROIs/MorelAtlasMNI152

cd right-vols-1mm/

#anterior nucleus
3dcalc -prefix morel_rh_a.nii.gz \
-a AD.nii.gz \
-b AM.nii.gz \
-c AV.nii.gz \
-d LD.nii.gz \
-expr 'a*10+b*11+c*12+d*13' 

#medial group
3dcalc -prefix morel_rh_m.nii.gz \
-a MDmc.nii.gz \
-b MDpc.nii.gz \
-c MV.nii.gz \
-d CL.nii.gz \
-e CeM.nii.gz \
-f CM.nii.gz \
-g Pf.nii.gz \
-expr 'a*20+b*21+c*22+d*23+e*24+f*25+g*26'

#posterior group
3dcalc -prefix morel_rh_p.nii.gz \
-a PuM.nii.gz \
-b PuI.nii.gz \
-c PuL.nii.gz \
-d PuA.nii.gz \
-e LP.nii.gz \
-f MGN.nii.gz \
-g SG.nii.gz \
-h Li.nii.gz \
-i Po.nii.gz \
-j LGNmc.nii.gz \
-k LGNpc.nii.gz \
-expr 'a*30+b*31+c*32+d*33+e*34+f*35+g*36+h*37+i*38+j*39+k*40'

#lateral group
3dcalc -prefix morel_rh_l.nii.gz \
-a VPLa.nii.gz \
-b VPLp.nii.gz \
-c VPM.nii.gz \
-d VPI.nii.gz \
-e VLa.nii.gz \
-f VLpd.nii.gz \
-g VLpv.nii.gz \
-h VAmc.nii.gz \
-i VApc.nii.gz \
-j VM.nii.gz \
-expr 'a*41+b*42+c*43+d*44+e*45+f*46+g*47+h*48+i*49+j*50'

3dcalc -prefix morel_rh.nii.gz \
-a morel_rh_a.nii.gz \
-b morel_rh_m.nii.gz \
-c morel_rh_p.nii.gz \
-d morel_rh_l.nii.gz \
-expr 'a+b+c+d'

#
cd ../left-vols-1mm/

#anterior nucleus
3dcalc -prefix morel_lh_a.nii.gz \
-a AD.nii.gz \
-b AM.nii.gz \
-c AV.nii.gz \
-d LD.nii.gz \
-expr 'a*10+b*11+c*12+d*13' 

#medial group
3dcalc -prefix morel_lh_m.nii.gz \
-a MDmc.nii.gz \
-b MDpc.nii.gz \
-c MV.nii.gz \
-d CL.nii.gz \
-e CeM.nii.gz \
-f CM.nii.gz \
-g Pf.nii.gz \
-expr 'a*20+b*21+c*22+d*23+e*24+f*25+g*26'

#posterior group
3dcalc -prefix morel_lh_p.nii.gz \
-a PuM.nii.gz \
-b PuI.nii.gz \
-c PuL.nii.gz \
-d PuA.nii.gz \
-e LP.nii.gz \
-f MGN.nii.gz \
-g SG.nii.gz \
-h Li.nii.gz \
-i Po.nii.gz \
-j LGNmc.nii.gz \
-k LGNpc.nii.gz \
-expr 'a*30+b*31+c*32+d*33+e*34+f*35+g*36+h*37+i*38+j*39+k*40'

#lateral group
3dcalc -prefix morel_lh_l.nii.gz \
-a VPLa.nii.gz \
-b VPLp.nii.gz \
-c VPM.nii.gz \
-d VPI.nii.gz \
-e VLa.nii.gz \
-f VLpd.nii.gz \
-g VLpv.nii.gz \
-h VAmc.nii.gz \
-i VApc.nii.gz \
-j VM.nii.gz \
-expr 'a*41+b*42+c*43+d*44+e*45+f*46+g*47+h*48+i*49+j*50'

3dcalc -prefix morel_lh.nii.gz \
-a morel_lh_a.nii.gz \
-b morel_lh_m.nii.gz \
-c morel_lh_p.nii.gz \
-d morel_lh_l.nii.gz \
-expr 'a+b+c+d'

cd /home/despoB/kaihwang/Rest/ROIs/
3dcalc -prefix Morel.nii.gz \
-a /home/despoB/kaihwang/Rest/ROIs/MorelAtlasMNI152/left-vols-1mm/morel_lh.nii.gz \
-b /home/despoB/kaihwang/Rest/ROIs/MorelAtlasMNI152/right-vols-1mm/morel_rh.nii.gz \
-expr 'a+b'


### get the borders
cd /home/despoB/kaihwang/Rest/ROIs/MorelAtlasMNI152/right-vols-1mm/

3dcalc -prefix morel_rh_a_mask.nii.gz \
-a AD.nii.gz \
-b AM.nii.gz \
-c AV.nii.gz \
-d LD.nii.gz \
-expr 'a+b+c+d' 

#medial group
3dcalc -prefix morel_rh_m_mask.nii.gz \
-a MDmc.nii.gz \
-b MDpc.nii.gz \
-c MV.nii.gz \
-d CL.nii.gz \
-e CeM.nii.gz \
-f CM.nii.gz \
-g Pf.nii.gz \
-expr 'a+b+c+d+e+f+g'

#posterior group
3dcalc -prefix morel_rh_p_mask.nii.gz \
-a PuM.nii.gz \
-b PuI.nii.gz \
-c PuL.nii.gz \
-d PuA.nii.gz \
-e LP.nii.gz \
-f MGN.nii.gz \
-g SG.nii.gz \
-h Li.nii.gz \
-i Po.nii.gz \
-j LGNmc.nii.gz \
-k LGNpc.nii.gz \
-expr 'a+b+c+d+e+f+g+h+i+j+k'

#lateral group
3dcalc -prefix morel_rh_l_mask.nii.gz \
-a VPLa.nii.gz \
-b VPLp.nii.gz \
-c VPM.nii.gz \
-d VPI.nii.gz \
-e VLa.nii.gz \
-f VLpd.nii.gz \
-g VLpv.nii.gz \
-h VAmc.nii.gz \
-i VApc.nii.gz \
-j VM.nii.gz \
-expr 'a+b+c+d+e+f+g+h+i+j'

3dcalc -prefix morel_rh_mask.nii.gz \
-a morel_rh_a_mask.nii.gz \
-b morel_rh_m_mask.nii.gz \
-c morel_rh_p_mask.nii.gz \
-d morel_rh_l_mask.nii.gz \
-expr 'a+b+c+d'

#left
cd /home/despoB/kaihwang/Rest/ROIs/MorelAtlasMNI152/left-vols-1mm/

#anterior nucleus
3dcalc -prefix morel_lh_a_mask.nii.gz \
-a AD.nii.gz \
-b AM.nii.gz \
-c AV.nii.gz \
-d LD.nii.gz \
-expr 'a+b+c+d' 

#medial group
3dcalc -prefix morel_lh_m_mask.nii.gz \
-a MDmc.nii.gz \
-b MDpc.nii.gz \
-c MV.nii.gz \
-d CL.nii.gz \
-e CeM.nii.gz \
-f CM.nii.gz \
-g Pf.nii.gz \
-expr 'a+b+c+d+e+f+g'

#posterior group
3dcalc -prefix morel_lh_p_mask.nii.gz \
-a PuM.nii.gz \
-b PuI.nii.gz \
-c PuL.nii.gz \
-d PuA.nii.gz \
-e LP.nii.gz \
-f MGN.nii.gz \
-g SG.nii.gz \
-h Li.nii.gz \
-i Po.nii.gz \
-j LGNmc.nii.gz \
-k LGNpc.nii.gz \
-expr 'a+b+c+d+e+f+g+h+i+j+k'

#lateral group
3dcalc -prefix morel_lh_l_mask.nii.gz \
-a VPLa.nii.gz \
-b VPLp.nii.gz \
-c VPM.nii.gz \
-d VPI.nii.gz \
-e VLa.nii.gz \
-f VLpd.nii.gz \
-g VLpv.nii.gz \
-h VAmc.nii.gz \
-i VApc.nii.gz \
-j VM.nii.gz \
-expr 'a+b+c+d+e+f+g+h+i+j'

3dcalc -prefix morel_lh_mask.nii.gz \
-a morel_lh_a_mask.nii.gz \
-b morel_lh_m_mask.nii.gz \
-c morel_lh_p_mask.nii.gz \
-d morel_lh_l_mask.nii.gz \
-expr 'a+b+c+d'

cd /home/despoB/kaihwang/Rest/ROIs/
3dcalc -prefix Morel_mask.nii.gz \
-a /home/despoB/kaihwang/Rest/ROIs/MorelAtlasMNI152/left-vols-1mm/morel_lh_mask.nii.gz \
-b /home/despoB/kaihwang/Rest/ROIs/MorelAtlasMNI152/right-vols-1mm/morel_rh_mask.nii.gz \
-expr 'a+b'
