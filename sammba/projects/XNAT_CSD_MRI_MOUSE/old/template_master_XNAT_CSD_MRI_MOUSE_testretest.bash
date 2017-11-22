#!/bin/bash

#bash -xef /home/nadkarni/git/sammba-mri/sammba/projects/XNAT_CSD_MRI_MOUSE_testretest/template_master_XNAT_CSD_MRI_MOUSE_testretest.bash 2>&1 | tee /home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/XNAT_CSD_MRI_MOUSE_testretest/$(date +%Y%m%d_%H%M%S).log

export AFNI_DECONFLICT=OVERWRITE

PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/common
PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/projects/XNAT_CSD_MRI_MOUSE_testretest
export PATH

PYTHONPATH=$PYTHONPATH:$PATH
export PYTHONPATH

projectdir=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/XNAT_CSD_MRI_MOUSE_testretest
savedir1=$projectdir/20170925_goodbaselinetemplate; mkdir -p $savedir1
savedir2=$projectdir/20170925_allgoodtemplate; mkdir -p $savedir2
savedir3=$projectdir/20170925_allbaselinetemplate; mkdir -p $savedir3
savedir4=$projectdir/20170925_goodposttemplate; mkdir -p $savedir4

find $projectdir/dl/baseline -maxdepth 1 -mindepth 1 -type d > $projectdir/allbaseline.txt

conv=0.005
twoblur=1
brainvol=400
#-rbt values might be improveable
Urad=18.3
b=70
t=80

template_maker () {

	savedir=$1; anatlist=$2; brainvol=$3; Urad=$4; b=$5; t=$6 conv=$7; twoblur=$8

	echo savedir $savedir
	echo anatlist $anatlist
	echo brainvol $brainvol
	echo Urad $Urad
	echo b $b
	echo t $t
	echo conv $conv
	echo twoblur $twoblur
	
	bash template_convert_XNAT_CSD_MRI_MOUSE_testretest.bash $savedir $anatlist
	bash MRIT2_extrcen.bash $savedir $brainvol $Urad $b $t
	bash MRIT3_shr.bash $savedir $savedir/UnBmBeCC_mean.nii.gz 1 $conv $twoblur
	bash MRIT3_shr.bash $savedir $savedir/shr1_mean.nii.gz 2 $conv $twoblur
	bash MRIT4_aff.bash $savedir $savedir/shr2_meanhead.nii.gz 2 3 $conv $twoblur $savedir/shr2_count.nii.gz

	3dmask_tool -union -inputs $savedir/aff3_video.nii.gz -prefix $savedir/aff3_unionmask.nii.gz
	3dmask_tool -dilate_inputs 4 -inputs $savedir/aff3_unionmask.nii.gz -prefix $savedir/aff3_unionmaskdil4.nii.gz
	weight=$savedir/aff3_unionmaskdil4.nii.gz

	bash MRIT5_Qw.bash $savedir $savedir/aff3_meanhead.nii.gz $weight UnCC.nii.gz UnBmBeCCAl2UnCCAl3.aff12.1D 0 4 1
	bash MRIT5_Qw.bash $savedir $savedir/Qw1_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw1_WARP.nii.gz 5 7 2
	bash MRIT6_Qw.bash $savedir $savedir/Qw2_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw2_WARP.nii.gz 8 13 3
	bash MRIT6_Qw.bash $savedir $savedir/Qw3_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw3_WARP.nii.gz 11 9 4
	
	for result in raw_mean UnCC_mean shr1_meanhead shr2_meanhead aff3_meanhead	 \
		Qw1_meanhead Qw2_meanhead Qw3_meanhead Qw4_meanhead; do
		meanSI=$(3dmaskave -q $savedir/$result.nii.gz)
		3dcalc -a $savedir/$result.nii.gz -expr "a*1000/$meanSI" -prefix $savedir/"meanSIfixed"_$result.nii.gz
	done

	3dTcat -prefix $savedir/progression.nii.gz									 \
				   $savedir/"meanSIfixed"_raw_mean.nii.gz						 \
				   $savedir/"meanSIfixed"_UnCC_mean.nii.gz						 \
				   $savedir/"meanSIfixed"_shr*_meanhead.nii.gz					 \
				   $savedir/"meanSIfixed"_aff3_meanhead.nii.gz					 \
				   $savedir/"meanSIfixed"_Qw*_meanhead.nii.gz

	rm -f $savedir/meanSIfixed_*
	
}

template_maker $savedir1 $projectdir/goodforbaselinetemplate.txt $brainvol $Urad $b $t $conv $twoblur
template_maker $savedir2 $projectdir/goodforalltemplate.txt $brainvol $Urad $b $t $conv $twoblur
template_maker $savedir3 $projectdir/allbaseline.txt $brainvol $Urad $b $t $conv $twoblur
template_maker $savedir4 $projectdir/goodforposttemplate.txt $brainvol $Urad $b $t $conv $twoblur

MIC40C57=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MRIatlases/MIC40C57/20170223
dorrbrain=$MIC40C57/brain100.nii.gz
dorratlas=$MIC40C57/labels100.nii.gz
3dAllineate -source $dorrbrain -base $savedir4/aff3_mean.nii.gz -prefix $savedir4/dorrbrainAl.nii.gz -1Dmatrix_save $savedir4/dorrbrainAl.aff12.1D
3dQwarp -base $savedir4/Qw4_meanhead.nii.gz -source $savedir4/dorrbrainAl.nii.gz -prefix $savedir4/dorrbrainAlQw.nii.gz -noneg -blur 0 -nmi
3dNwarpCat -prefix $savedir4/dorrbrainAlQwcat_WARP.nii.gz $savedir4/dorrbrainAlQw_WARP.nii.gz $savedir4/dorrbrainAl.aff12.1D
3dNwarpApply -nwarp $savedir4/dorrbrainAlQwcat_WARP.nii.gz -source $MIC40C57/labels100.nii.gz -master $savedir4/Qw4_meanhead.nii.gz -ainterp NN -prefix $savedir4/dorratlasAlQw.nii.gz

3dresample -dxyz 0.2 0.2 0.2 -prefix $savedir4/Qw4_meanhead200.nii.gz -input $savedir4/Qw4_meanhead.nii.gz
