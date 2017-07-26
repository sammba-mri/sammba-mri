#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE #saves a lot of grief

savedir=$(readlink -e $1)
template=$(readlink -e $2)
weight=$(readlink -e $3)
source=$4
iniwarp=$5
inilev=$6
minpatch=$7
n=$8

cd $savedir

for dir in $(find -mindepth 1 -maxdepth 1 -type d); do

	if [[ $iniwarp == IDENT ]]; then
		iniwarp2="IDENT($template)"
	else
		iniwarp2=$dir/$iniwarp
	fi
	
	3dQwarp -base $template -source $dir/$source -prefix $dir/UnCCQw$n.nii.gz -noneg -weight $weight -iniwarp $iniwarp2 -inilev $inilev -minpatch $minpatch
	
done

3dTcat -prefix $savedir/Qw"$n"_videohead.nii.gz $(find $savedir -name UnCCQw"$n".nii.gz)
3dTstat -prefix $savedir/Qw"$n"_meanhead.nii.gz $savedir/Qw"$n"_videohead.nii.gz
