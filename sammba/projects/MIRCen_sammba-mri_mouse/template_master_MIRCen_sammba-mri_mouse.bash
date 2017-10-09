#!/bin/bash

#bash -xef /home/nadkarni/git/sammba-mri/sammba/projects/MIRCen_sammba-mri_mouse/template_master_MIRCen_sammba-mri_mouse.bash 2>&1 | tee /home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MIRCen_sammba-mri_mouse/$(date +%Y%m%d_%H%M%S).log

export AFNI_DECONFLICT=OVERWRITE

PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/common
PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/projects/MIRCen_sammba-mri_mouse
export PATH

PYTHONPATH=$PYTHONPATH:$PATH
export PYTHONPATH

projectdir=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MIRCen_sammba-mri_mouse
savedir=$projectdir/20170925_template
mkdir $savedir

dicomdirlist=$projectdir/dicomdirs.txt

finres=0.1
conv=0.005
twoblur=1
brainvol=400
#-rbt values might be improveable
Urad=18.3
b=70
t=80

bash template_convert_MIRCen_sammba-mri_mouse.bash $dicomdirlist $savedir 0.1
bash MRIT2_extrcen.bash $savedir $brainvol $Urad $b $t
bash MRIT3_shr.bash $savedir $savedir/UnBmBeCC_mean.nii.gz 1 $conv $twoblur
bash MRIT3_shr.bash $savedir $savedir/shr1_mean.nii.gz 2 $conv $twoblur
bash MRIT4_aff.bash $savedir $savedir/shr2_meanhead.nii.gz 2 3 $conv $twoblur $savedir/shr2_count.nii.gz

3dmask_tool -frac 0.5 -inputs $savedir/aff3_video.nii.gz -prefix $savedir/aff3_frac05mask.nii.gz
3dcalc -a $savedir/aff3_meanhead.nii.gz -expr 'step(0.5-(x+3.95)*(x+3.95)-(y-3.75)*(y-3.75)-(z+0.75)*(z+0.75))' -prefix $savedir/rightcolliculus.nii.gz
3dcalc -a $savedir/aff3_meanhead.nii.gz -expr 'step(0.5-(x-3.75)*(x-3.75)-(y-3.75)*(y-3.75)-(z+0.85)*(z+0.85))' -prefix $savedir/leftcolliculus.nii.gz
3dmask_tool -union -inputs $savedir/aff3_frac05mask.nii.gz $savedir/rightcolliculus.nii.gz $savedir/leftcolliculus.nii.gz -prefix $savedir/aff3_unionmask.nii.gz
3dmask_tool -dilate_inputs 5 -inputs $savedir/aff3_unionmask.nii.gz -prefix $savedir/aff3_unionmaskdil5.nii.gz
weight=$savedir/aff3_unionmaskdil5.nii.gz

bash MRIT5_Qw.bash $savedir $savedir/aff3_meanhead.nii.gz $weight UnCC.nii.gz UnBmBeCCAl2UnCCAl3.aff12.1D 0 4 1
bash MRIT5_Qw.bash $savedir $savedir/Qw1_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw1_WARP.nii.gz 5 7 2
bash MRIT6_Qw.bash $savedir $savedir/Qw2_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw2_WARP.nii.gz 8 13 3
bash MRIT6_Qw.bash $savedir $savedir/Qw3_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw3_WARP.nii.gz 11 9 4

MIC40C57=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MRIatlases/MIC40C57/20170223
dorrbrain=$MIC40C57/brain100.nii.gz
dorratlas=$MIC40C57/labels100.nii.gz
3dAllineate -source $dorrbrain -base $savedir/aff3_mean.nii.gz -prefix $savedir/dorrbrainAl.nii.gz -1Dmatrix_save $savedir/dorrbrainAl.aff12.1D
3dQwarp -base $savedir/Qw4_meanhead.nii.gz -source $savedir/dorrbrainAl.nii.gz -prefix $savedir/dorrbrainAlQw.nii.gz -noneg -blur 0 -nmi
3dNwarpCat -prefix $savedir/dorrbrainAlQwcat_WARP.nii.gz $savedir/dorrbrainAlQw_WARP.nii.gz $savedir/dorrbrainAl.aff12.1D
3dNwarpApply -nwarp $savedir/dorrbrainAlQwcat_WARP.nii.gz -source $MIC40C57/labels100.nii.gz -master $savedir/Qw4_meanhead.nii.gz -ainterp NN -prefix $savedir/dorratlasAlQw.nii.gz

3dresample -dxyz 0.2 0.2 0.2 -prefix $savedir/Qw4_meanhead200.nii.gz -input $savedir/Qw4_meanhead.nii.gz
