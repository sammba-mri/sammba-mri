#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE

savedir=$(readlink -e $1)
n=$2
atlas=$3
atlaswarp=$4

cd $savedir

for dir in $(ls -d */ | sed 's/\///g'); do
	
	3dNwarpApply -nwarp $dir/UnCCQw"$n"_WARP.nii.gz -source $dir/T2anatCC.nii.gz -master $savedir/Qw"$n"_meanhead.nii.gz -prefix $dir/Na.nii.gz
	3dNwarpCat -prefix $dir/atlas_WARPINV.nii.gz "INV("$dir"/UnCCQw"$n"_WARP.nii.gz)" $atlaswarp
	3dNwarpApply -nwarp $dir/atlas_WARPINV.nii.gz -source $atlas -master $dir/T2anatCC.nii.gz -ainterp NN -prefix $dir/atlas_Na.nii.gz
		
done

3dTcat -prefix $savedir/Na_videohead.nii.gz $(find $savedir -name 'Na.nii.gz')
3dTstat -prefix $savedir/Na_meanhead.nii.gz $savedir/Na_videohead.nii.gz
