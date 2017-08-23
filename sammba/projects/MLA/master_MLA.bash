#!/bin/bash

#bash -xef /home/nadkarni/git/sammba-mri/sammba/projects/MLA/master_MLA.bash 2>&1 | tee /home/Pmamobipet/Tvx-Manips-MD_/MD_1701-Microcebe-Creation-Atlas/$(date +%Y%m%d_%H%M%S).log

export AFNI_DECONFLICT=OVERWRITE

PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/common
PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/projects/MLA
export PATH

projectdir=/home/Pmamobipet/Tvx-Manips-MD_/MD_1701-Microcebe-Creation-Atlas
rawdatadir=$projectdir/MRI-Images-Brutes
savedir=$projectdir/processed_20170731
mkdir $savedir

finres=0.115
conv=0.01
twoblur=2
brainvol=1600
#-rbt values might be improveable
Urad=18.3
b=70
t=80

bash convert_MLA.bash $rawdatadir $savedir $finres
bash MRIT2_extrcen.bash $savedir $brainvol $Urad $b $t

3dcopy $projectdir/bonsatlasprocessing/bonstack_dupflipcomb.hdr $savedir/bons.nii.gz
echo 1 0 0 0 0 0 1 0 0 -1 0 0 > $savedir/x90.1d
echo -1 0 0 0 0 1 0 0 0 0 -1 0 > $savedir/y180.1d
echo -1 0 0 0 -0 -1 0 0 0 0 1 0 > $savedir/z180.1d
nifti_tool -disp_hdr -field srow_x -field srow_y -field srow_z -infiles $savedir/bons.nii.gz -quiet > $savedir/sform.txt
fsl5.0-fslorient -setsform $(cat_matvec -ONELINE $savedir/sform.txt $savedir/x90.1d $savedir/y180.1d $savedir/z180.1d) 0 0 0 1 $savedir/bons.nii.gz
fsl5.0-fslorient -copysform2qform $savedir/bons.nii.gz
3dCM -set 0 0 0 $savedir/bons.nii.gz
3dresample -master $savedir/UnCC_mean.nii.gz -rmode Cubic -prefix $savedir/bons.nii.gz -inset $savedir/bons.nii.gz
bons=$savedir/bons.nii.gz

bash MRIT3_shr.bash $savedir $bons 1 $conv $twoblur
bash MRIT3_shr.bash $savedir $savedir/shr1_mean.nii.gz 2 $conv $twoblur
bash MRIT4_aff.bash $savedir $savedir/shr2_meanhead.nii.gz 2 3 $conv $twoblur $savedir/shr2_count.nii.gz

3dcalc -a $savedir/aff3_meanhead.nii.gz -expr 'step(2-(x+6.15)*(x+6.15)-(y+4.08)*(y+4.08)-(z-5)*(z-5))' -prefix $savedir/rightcolliculus.nii.gz
3dmask_tool -union -inputs $savedir/aff3_video.nii.gz $savedir/rightcolliculus.nii.gz -prefix $savedir/aff3_unionmask.nii.gz
3dmask_tool -dilate_inputs 7 -inputs $savedir/aff3_unionmask.nii.gz -prefix $savedir/aff3_unionmaskdil7.nii.gz
weight=$savedir/aff3_unionmaskdil7.nii.gz

bash MRIT5_Qw.bash $savedir $savedir/aff3_meanhead.nii.gz $weight UnCC.nii.gz UnBmBeCCAl2UnCCAl3.aff12.1D 0 4 1
bash MRIT5_Qw.bash $savedir $savedir/Qw1_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw1_WARP.nii.gz 5 7 2
bash MRIT6_Qw.bash $savedir $savedir/Qw2_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw2_WARP.nii.gz 8 13 3
bash MRIT6_Qw.bash $savedir $savedir/Qw3_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw3_WARP.nii.gz 11 9 4
bash MRIT7_origproc.bash $savedir 4
bash MRIT8_TBM.bash $savedir UnCCQw4_WARP.nii.gz

3dMean -stdev -prefix $savedir/bulkstdev.nii.gz $(find $savedir -name 'bulk.nii.gz')
3dZeropad -master $savedir/aff3_unionmask.nii.gz -prefix $savedir/bulkstdevZP.nii.gz $savedir/bulkstdev.nii.gz
3dcalc -a $savedir/aff3_unionmask.nii.gz -b $savedir/bulkstdevZP.nii.gz -expr 'ispositive(a)*b' -prefix $savedir/bulkstdevZP_masked.nii.gz
