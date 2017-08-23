#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE #saves a lot of grief

savedir=$(readlink -e $1)
template=$(readlink -e $2)
n=$3
conv=$4
twoblur=$5

cd $savedir

for dir in $(find -mindepth 1 -maxdepth 1 -type d); do
	
	3dAllineate -base $template -source $dir/UnBmBeCC.nii.gz -prefix $dir/UnBmBeCCAl$n.nii.gz -1Dmatrix_save $dir/UnBmBeCCAl$n.aff12.1D -conv $conv -twoblur $twoblur -warp shr
	3dAllineate -input $dir/UnCC.nii.gz -master $template -prefix $dir/UnCCAa$n.nii.gz -1Dmatrix_apply $dir/UnBmBeCCAl$n.aff12.1D
	
done

3dTcat -prefix $savedir/shr"$n"_video.nii.gz $(find $savedir -name UnBmBeCCAl"$n".nii.gz)
3dTstat -prefix $savedir/shr"$n"_mean.nii.gz $savedir/shr"$n"_video.nii.gz
3dTcat -prefix $savedir/shr"$n"_videohead.nii.gz $(find $savedir -name UnCCAa"$n".nii.gz)
3dTstat -prefix $savedir/shr"$n"_meanhead.nii.gz $savedir/shr"$n"_videohead.nii.gz
3dmask_tool -count -inputs $savedir/shr"$n"_video.nii.gz -prefix $savedir/shr"$n"_count.nii.gz
