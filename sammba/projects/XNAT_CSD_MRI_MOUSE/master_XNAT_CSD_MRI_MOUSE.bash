#!/bin/bash

#bash -xef /home/nadkarni/git/sammba-mri/sammba/projects/XNAT_CSD_MRI_MOUSE/master_XNAT_CSD_MRI_MOUSE.bash 2>&1 | tee /home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/XNAT_CSD_MRI_MOUSE/$(date +%Y%m%d_%H%M%S).log

export AFNI_DECONFLICT=OVERWRITE

PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/common
PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/projects/XNAT_CSD_MRI_MOUSE
export PATH

PYTHONPATH=$PYTHONPATH:$PATH
export PYTHONPATH

projectdir=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/XNAT_CSD_MRI_MOUSE
savedir=$projectdir/analysis20171018
mkdir $savedir
tsavedir1=$savedir/goodbaselinetemplate; mkdir -p $tsavedir1
tsavedir2=$savedir/allgoodtemplate; mkdir -p $tsavedir2
tsavedir3=$savedir/allbaselinetemplate; mkdir -p $tsavedir3
tsavedir4=$savedir/goodposttemplate; mkdir -p $tsavedir4

find $projectdir/dl/baseline -maxdepth 1 -mindepth 1 -type d > $projectdir/allbaseline.txt

conv=0.005
twoblur=1
brainvol=400
#-rbt values might be improveable
Urad=18.3
b=70
t=80

template_maker () {

	tsavedir=$1; anatlist=$2; brainvol=$3; Urad=$4; b=$5; t=$6 conv=$7; twoblur=$8

	echo tsavedir $tsavedir
	echo anatlist $anatlist
	echo brainvol $brainvol
	echo Urad $Urad
	echo b $b
	echo t $t
	echo conv $conv
	echo twoblur $twoblur
	
	#orientation correctors for use later
	#https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Euler_angles_.28_z-y.E2.80.99-x.E2.80.B3_intrinsic.29_.E2.86.92_Rotation_matrix
	echo 1 0 0 0 0 0 1 0 0 -1 0 0 > $tsavedir/x270.1d
	echo -1 0 0 0 0 1 0 0 0 0 -1 0 > $tsavedir/y180.1d

	nlines=$(wc -l $anatlist | awk '{print $1}')

	for ((a=1; a<=$nlines; a++)); do

		rawdir=$(head -n$a < $anatlist | tail -n1)
		stem=$(basename $rawdir)
		prefix=$(basename $(dirname $rawdir))
		newniftidir=$tsavedir/"$prefix"_"$stem"
		mkdir -p $newniftidir
		cp $rawdir/3DRARE.nii.gz $newniftidir/anat.nii.gz
		nifti_tool -disp_hdr -field srow_x -field srow_y -field srow_z -infiles $newniftidir/anat.nii.gz -quiet > $newniftidir/sform.txt
		fsl5.0-fslorient -setsform $(cat_matvec -ONELINE $newniftidir/sform.txt $tsavedir/x270.1d $tsavedir/y180.1d) 0 0 0 1 $newniftidir/anat.nii.gz
		fsl5.0-fslorient -copysform2qform $newniftidir/anat.nii.gz
	
	done

	3dTcat -prefix $tsavedir/raw_video.nii.gz $(find $tsavedir -name anat.nii.gz)
	3dTstat -prefix $tsavedir/raw_mean.nii.gz $tsavedir/raw_video.nii.gz
	
	bash MRIT2_extrcen.bash $tsavedir $brainvol $Urad $b $t
	bash MRIT3_shr.bash $tsavedir $tsavedir/UnBmBeCC_mean.nii.gz 1 $conv $twoblur
	bash MRIT3_shr.bash $tsavedir $tsavedir/shr1_mean.nii.gz 2 $conv $twoblur
	bash MRIT4_aff.bash $tsavedir $tsavedir/shr2_meanhead.nii.gz 2 3 $conv $twoblur $tsavedir/shr2_count.nii.gz

	3dmask_tool -union -inputs $tsavedir/aff3_video.nii.gz -prefix $tsavedir/aff3_unionmask.nii.gz
	3dmask_tool -dilate_inputs 4 -inputs $tsavedir/aff3_unionmask.nii.gz -prefix $tsavedir/aff3_unionmaskdil4.nii.gz
	weight=$tsavedir/aff3_unionmaskdil4.nii.gz

	bash MRIT5_Qw.bash $tsavedir $tsavedir/aff3_meanhead.nii.gz $weight UnCC.nii.gz UnBmBeCCAl2UnCCAl3.aff12.1D 0 4 1
	bash MRIT5_Qw.bash $tsavedir $tsavedir/Qw1_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw1_WARP.nii.gz 5 7 2
	bash MRIT6_Qw.bash $tsavedir $tsavedir/Qw2_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw2_WARP.nii.gz 8 13 3
	bash MRIT6_Qw.bash $tsavedir $tsavedir/Qw3_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw3_WARP.nii.gz 11 9 4
	
	for result in raw_mean UnCC_mean shr1_meanhead shr2_meanhead aff3_meanhead	 \
		Qw1_meanhead Qw2_meanhead Qw3_meanhead Qw4_meanhead; do
		meanSI=$(3dmaskave -q $tsavedir/$result.nii.gz)
		3dcalc -a $tsavedir/$result.nii.gz -expr "a*1000/$meanSI" -prefix $tsavedir/"meanSIfixed"_$result.nii.gz
	done

	3dTcat -prefix $tsavedir/progression.nii.gz									 \
				   $tsavedir/"meanSIfixed"_raw_mean.nii.gz						 \
				   $tsavedir/"meanSIfixed"_UnCC_mean.nii.gz						 \
				   $tsavedir/"meanSIfixed"_shr*_meanhead.nii.gz					 \
				   $tsavedir/"meanSIfixed"_aff3_meanhead.nii.gz					 \
				   $tsavedir/"meanSIfixed"_Qw*_meanhead.nii.gz

	rm -f $tsavedir/meanSIfixed_*
	
}

template_maker $tsavedir1 $projectdir/goodforbaselinetemplate.txt $brainvol $Urad $b $t $conv $twoblur
template_maker $tsavedir2 $projectdir/goodforalltemplate.txt $brainvol $Urad $b $t $conv $twoblur
template_maker $tsavedir3 $projectdir/allbaseline.txt $brainvol $Urad $b $t $conv $twoblur
template_maker $tsavedir4 $projectdir/goodforposttemplate.txt $brainvol $Urad $b $t $conv $twoblur

MIC40C57=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MRIatlases/MIC40C57/20170223
dorrbrain=$MIC40C57/brain100.nii.gz
dorratlas=$MIC40C57/labels100.nii.gz
3dAllineate -source $dorrbrain -base $tsavedir4/aff3_mean.nii.gz -prefix $tsavedir4/dorrbrainAl.nii.gz -1Dmatrix_save $tsavedir4/dorrbrainAl.aff12.1D
3dQwarp -base $tsavedir4/Qw4_meanhead.nii.gz -source $tsavedir4/dorrbrainAl.nii.gz -prefix $tsavedir4/dorrbrainAlQw.nii.gz -noneg -blur 0 -nmi
3dNwarpCat -prefix $tsavedir4/dorrbrainAlQwcat_WARP.nii.gz $tsavedir4/dorrbrainAlQw_WARP.nii.gz $tsavedir4/dorrbrainAl.aff12.1D
3dNwarpApply -nwarp $tsavedir4/dorrbrainAlQwcat_WARP.nii.gz -source $MIC40C57/labels100.nii.gz -master $tsavedir4/Qw4_meanhead.nii.gz -ainterp NN -prefix $tsavedir4/dorratlasAlQw.nii.gz
3dresample -dxyz 0.2 0.2 0.2 -prefix $tsavedir4/Qw4_meanhead200.nii.gz -input $tsavedir4/Qw4_meanhead.nii.gz
3dNwarpApply -nwarp $tsavedir4/dorrbrainAlQwcat_WARP.nii.gz -iwarp -source $tsavedir4/Qw4_meanhead.nii.gz -master $dorrbrain -prefix $tsavedir4/Qw4_meanhead_Na.nii.gz
3dresample -dxyz 0.2 0.2 0.2 -rmode Cu -prefix $tsavedir4/Qw4_meanhead200_Na.nii.gz -input $tsavedir4/Qw4_meanhead200.nii.gz
3dresample -dxyz 0.2 0.2 0.2 -rmode Cu -prefix $tsavedir4/Qw4_meanhead_Na_200.nii.gz -input $tsavedir4/Qw4_meanhead_Na.nii.gz
3dresample -master $tsavedir4/Qw4_meanhead_Na_200.nii.gz -prefix $tsavedir4/labels200.nii.gz -input $dorratlas
3dNwarpApply -nwarp $tsavedir4/dorrbrainAlQwcat_WARP.nii.gz -iwarp -source $tsavedir4/aff3_unionmaskdil4.nii.gz -master $tsavedir4/Qw4_meanhead_Na.nii.gz -prefix $tsavedir4/aff3_unionmaskdil4_Na.nii.gz
