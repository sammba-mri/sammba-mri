#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE #saves a lot of grief

dicomdirlist=$(readlink -e $1)
savedir=$(readlink -e $2)
resampleres=$3

nlines=$(wc -l $dicomdirlist | awk '{print $1}')

for ((a=1; a<=$nlines; a++)); do

	dicomdir=$(head -n$a < $dicomdirlist | tail -n1)
	NIfTIdir=$savedir/$(basename $dicomdir)
	mkdir $NIfTIdir
	python -m PVEnDCMtoNIfTI.py /usr/bin/dcmdump $dicomdir $NIfTIdir yes
	IDsequencetypes.bash $NIfTIdir
	
	3dresample -dxyz $resampleres $resampleres $resampleres -rmode Cubic -prefix $NIfTIdir/anat.nii.gz -inset $NIfTIdir/anat_n0.nii.gz
	3dZeropad -RL 150 -AP 200 -IS 100 -prefix $NIfTIdir/anat.nii.gz $NIfTIdir/anat.nii.gz

	rm -f $(find $NIfTIdir -mindepth 1 ! -name anat.nii.gz)
	
done

#create raw data video and mean
3dTcat -prefix $savedir/raw_video.nii.gz $(find $savedir -name 'anat.nii.gz')
3dTstat -prefix $savedir/raw_mean.nii.gz $savedir/raw_video.nii.gz
