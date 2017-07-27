#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE #saves a lot of grief

savedir=$(readlink -e $1)
template=$(readlink -e $2)
n1=$3
n2=$4
conv=$5
twoblur=$6
weight=$7

cd $savedir

for dir in $(find -mindepth 1 -maxdepth 1 -type d); do
	
	3dAllineate -base $template -source $dir/UnCCAa"$n1".nii.gz -prefix $dir/UnCCAl"$n2".nii.gz -1Dmatrix_save $dir/UnCCAl"$n2".aff12.1D -conv $conv -onepass -twoblur $twoblur -weight $weight
	cat_matvec -ONELINE $dir/UnBmBeCCAl"$n1".aff12.1D $dir/UnCCAl"$n2".aff12.1D > $dir/UnBmBeCCAl"$n1"UnCCAl"$n2".aff12.1D
	3dAllineate -input $dir/UnCC.nii.gz -master $template -prefix $dir/UnCCAa"$n2".nii.gz -1Dmatrix_apply $dir/UnBmBeCCAl"$n1"UnCCAl"$n2".aff12.1D
	3dAllineate -input $dir/UnBmBeCC.nii.gz -master $template -prefix $dir/UnBmBeCCAa"$n2".nii.gz -1Dmatrix_apply $dir/UnBmBeCCAl"$n1"UnCCAl"$n2".aff12.1D
	
done

3dTcat -prefix $savedir/aff"$n2"_video.nii.gz $(find $savedir -name UnBmBeCCAa"$n2".nii.gz)
3dTstat -prefix $savedir/aff"$n2"_mean.nii.gz $savedir/aff"$n2"_video.nii.gz
3dTcat -prefix $savedir/aff"$n2"_videohead.nii.gz $(find $savedir -name UnCCAa"$n2".nii.gz)
3dTstat -prefix $savedir/aff"$n2"_meanhead.nii.gz $savedir/aff"$n2"_videohead.nii.gz
3dmask_tool -count -inputs $savedir/aff"$n2"_video.nii.gz -prefix $savedir/aff"$n2"_count.nii.gz
