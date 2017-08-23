#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE #saves a lot of grief

rawdatadir=$(readlink -e $1)
savedir=$(readlink -e $2)
resampleres=$3

cd $rawdatadir

for dir in $(find -mindepth 1 -maxdepth 1 -type d) ; do

	NIfTIdir=$savedir/$dir
	mkdir $NIfTIdir
	python -m PVEnDCMtoNIfTI.py /usr/bin/dcmdump $dir $NIfTIdir yes
	IDsequencetypes.bash $NIfTIdir
	
	3dresample -dxyz $resampleres $resampleres $resampleres -rmode Cubic -prefix $NIfTIdir/anat.nii.gz -inset $NIfTIdir/anat_n0.nii.gz
	
done

#create raw data video and mean
3dTcat -prefix $savedir/raw_video.nii.gz $(find $savedir -name 'anat.nii.gz')
3dTstat -prefix $savedir/raw_mean.nii.gz $savedir/raw_video.nii.gz
