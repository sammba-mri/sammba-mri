#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE

savedir=$(readlink -e $1)
brainvol=$2
R=$3
b=$4
t=$5

cd $savedir

for dir in $(find -mindepth 1 -maxdepth 1 -type d); do

	3dUnifize -prefix $dir/Un.nii.gz -input $dir/anat.nii.gz -rbt $R $b $t

	thresh=$(3dClipLevel $dir/Un.nii.gz)
	printf -v thresh %.0f "$thresh"
	echo "brain/background threshold=$thresh"

	RATS_MM $dir/Un.nii.gz $dir/UnBm.nii.gz -v $brainvol -t $thresh

	fsl5.0-fslmaths $dir/Un.nii.gz -mas $dir/UnBm.nii.gz $dir/UnBmBe.nii.gz
	
	3dCM -set 0 0 0 $dir/UnBmBe.nii.gz
	3drefit -duporigin $dir/UnBmBe.nii.gz $dir/anat.nii.gz
	3drefit -duporigin $dir/UnBmBe.nii.gz $dir/Un.nii.gz
	3drefit -duporigin $dir/UnBmBe.nii.gz $dir/UnBm.nii.gz
	
done

3dTcat -prefix $savedir/Un_video.nii.gz $(find $savedir -name 'Un.nii.gz')
3dTstat -prefix $savedir/Un_mean.nii.gz Un_video.nii.gz
3dTcat -prefix $savedir/UnBmBe_video.nii.gz $(find $savedir -name 'UnBmBe.nii.gz')
3dTstat -prefix $savedir/UnBmBe_mean.nii.gz UnBmBe_video.nii.gz

3dUndump -master $savedir/UnBmBe_mean.nii.gz -prefix $savedir/emptytemplate.nii.gz
3drefit -xorigin cen -yorigin cen -zorigin cen $savedir/emptytemplate.nii.gz

for dir in $(find -mindepth 1 -maxdepth 1 -type d); do

	3dresample -master $savedir/emptytemplate.nii.gz -prefix $dir/anatCC.nii.gz -inset $dir/anat.nii.gz
	3dresample -master $savedir/emptytemplate.nii.gz -prefix $dir/UnCC.nii.gz -inset $dir/Un.nii.gz
	3dresample -master $savedir/emptytemplate.nii.gz -prefix $dir/UnBmCC.nii.gz -inset $dir/UnBm.nii.gz
	3dresample -master $savedir/emptytemplate.nii.gz -prefix $dir/UnBmBeCC.nii.gz -inset $dir/UnBmBe.nii.gz
	
done

3dTcat -prefix $savedir/UnCC_video.nii.gz $(find $savedir -name 'UnCC.nii.gz')
3dTstat -prefix $savedir/UnCC_mean.nii.gz UnCC_video.nii.gz
3dTcat -prefix $savedir/UnBmBeCC_video.nii.gz $(find $savedir -name 'UnBmBeCC.nii.gz')
3dTstat -prefix $savedir/UnBmBeCC_mean.nii.gz UnBmBeCC_video.nii.gz
