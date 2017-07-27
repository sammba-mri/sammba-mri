#!/bin/bash

#bash -xef /home/nadkarni/git/sammba-mri/sammba/projects/DCLK3/master_DCLK3.bash 2>&1 | tee /home/Pplateforme/Plate-forme_RMN/NachiketNadkarni/DCLK3/$(date +%Y%m%d_%H%M%S).log

PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/common
PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/projects/DCLK3
export PATH

PYTHONPATH=$PYTHONPATH:$PATH
export PYTHONPATH

rawdatadir1=/home/Pplateforme/Plate-forme_RMN/JulienF/Plateforme_IRM/Raw_data/raw_data_117T/data_Souris_DCLK3/fevrier_2017
rawdatadir2=/home/Pplateforme/Plate-forme_RMN/JulienF/Plateforme_IRM/Raw_data/raw_data_117T/data_Souris_DCLK3/juillet_2017
savedir=/home/Pplateforme/Plate-forme_RMN/NachiketNadkarni/DCLK3/processed_20170727/MRIsessions
mkdir $savedir

finres=0.1
conv=0.005
twoblur=1
brainvol=400
#-rbt values might be improveable
Urad=18.3
b=70
t=80

bash convert_DCLK3.bash $rawdatadir1 $savedir 0.1
rm -f $savedir/raw_mean.nii.gz
rm -f $savedir/raw_video.nii.gz
bash convert_DCLK3.bash $rawdatadir2 $savedir 0.1

bash MRIT2_extrcen.bash $savedir $brainvol $Urad $b $t

template=$savedir/UnBmBeCC_mean.nii.gz

bash MRIT3_shr.bash $savedir $template 1 $conv $twoblur
bash MRIT3_shr.bash $savedir $savedir/shr1_mean.nii.gz 2 $conv $twoblur
bash MRIT4_aff.bash $savedir $savedir/shr2_meanhead.nii.gz 2 3 $conv $twoblur $savedir/shr2_count.nii.gz

3dmask_tool -inter -inputs $savedir/shr2_video.nii.gz -prefix $savedir/shr2_intermask.nii.gz
3dmask_tool -dilate_inputs 15 -inputs $savedir/shr2_intermask.nii.gz -prefix $savedir/shr2_intermaskdil15.nii.gz
weight=$savedir/shr2_intermaskdil15.nii.gz

bash MRIT5_Qw.bash $savedir $savedir/aff3_meanhead.nii.gz $weight UnCC.nii.gz UnBmBeCCAl2UnCCAl3.aff12.1D 0 4 1
bash MRIT5_Qw.bash $savedir $savedir/Qw1_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw1_WARP.nii.gz 5 7 2
bash MRIT6_Qw.bash $savedir $savedir/Qw2_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw2_WARP.nii.gz 8 13 3
bash MRIT6_Qw.bash $savedir $savedir/Qw3_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw3_WARP.nii.gz 11 9 4

MIC40C57=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MRIatlases/MIC40C57/20170223
dorrbrain=$MIC40C57/brain100.nii.gz
dorratlas=$MIC40C57/labels100.nii.gz
3dAllineate -source $dorrbrain -base $savedir/shr2_mean.nii.gz -prefix $savedir/dorrbrainAl.nii.gz -1Dmatrix_save $savedir/dorrbrainAl.aff12.1D
3dQwarp -base $savedir/Qw4_meanhead.nii.gz -source $savedir/dorrbrainAl.nii.gz -prefix $savedir/dorrbrainAlQw.nii.gz -noneg -blur 0 -nmi
3dNwarpCat -prefix $savedir/dorrbrainAlQwcat_WARP.nii.gz $savedir/dorrbrainAlQw_WARP.nii.gz $savedir/dorrbrainAl.aff12.1D
3dNwarpApply -nwarp $savedir/dorrbrainAlQwcat_WARP.nii.gz -source $MIC40C57/labels100.nii.gz -master $savedir/Qw4_meanhead.nii.gz -ainterp NN -prefix $savedir/dorratlasAlQw.nii.gz

bash MRIT7_origproc.bash $savedir 4 $dorratlas $savedir/dorrbrainAlQwcat_WARP.nii.gz
bash MRIT8_TBM.bash $savedir UnCCQw4_WARP.nii.gz

3dMean -prefix $savedir/bulkWT.nii.gz $(find $savedir -name 'bulk.nii.gz' | grep WT)
3dMean -prefix $savedir/bulkKO.nii.gz $(find $savedir -name 'bulk.nii.gz' | grep KO)
3dcalc -a $savedir/bulkKO.nii.gz -b $savedir/bulkWT.nii.gz -expr 'b-a' -prefix $savedir/bulkKOminusWT.nii.gz
3dZeropad -master $savedir/Qw4_meanhead.nii.gz -prefix $savedir/bulkKOminusWT_Zp.nii.gz $savedir/bulkKOminusWT.nii.gz

3dMean -prefix $savedir/shearWT.nii.gz $(find $savedir -name 'shear.nii.gz' | grep WT)
3dMean -prefix $savedir/shearKO.nii.gz $(find $savedir -name 'shear.nii.gz' | grep KO)
3dcalc -a $savedir/shearKO.nii.gz -b $savedir/shearWT.nii.gz -expr 'b-a' -prefix $savedir/shearKOminusWT.nii.gz
3dZeropad -master $savedir/Qw4_meanhead.nii.gz -prefix $savedir/shearKOminusWT_Zp.nii.gz $savedir/shearKOminusWT.nii.gz

3dMean -prefix $savedir/vorticityWT.nii.gz $(find $savedir -name 'vorticity.nii.gz' | grep WT)
3dMean -prefix $savedir/vorticityKO.nii.gz $(find $savedir -name 'vorticity.nii.gz' | grep KO)
3dcalc -a $savedir/vorticityKO.nii.gz -b $savedir/vorticityWT.nii.gz -expr 'b-a' -prefix $savedir/vorticityKOminusWT.nii.gz
3dZeropad -master $savedir/Qw4_meanhead.nii.gz -prefix $savedir/vorticityKOminusWT_Zp.nii.gz $savedir/vorticityKOminusWT.nii.gz
