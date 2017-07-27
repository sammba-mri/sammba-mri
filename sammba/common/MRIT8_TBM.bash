#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE

savedir=$(readlink -e $1)
warp=$2

cd $savedir

for dir in $(find -mindepth 1 -maxdepth 1 -type d); do
	
	3dNwarpFuncs -nwarp $dir/$warp -prefix $dir/bulk.nii.gz -bulk
	3dNwarpFuncs -nwarp $dir/$warp -prefix $dir/shear.nii.gz -shear
	3dNwarpFuncs -nwarp $dir/$warp -prefix $dir/vorticity.nii.gz -vorticity
	
done
