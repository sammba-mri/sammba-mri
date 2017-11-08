#!/bin/bash

#bash -xef /home/nadkarni/git/sammba-mri/sammba/projects/MIRCen_sammba-mri_mouse/master_MIRCen_sammba-mri_mouse.bash 2>&1 | tee /home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MIRCen_sammba-mri_mouse/$(date +%Y%m%d_%H%M%S).log

export AFNI_DECONFLICT=OVERWRITE

PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/common
export PATH

PYTHONPATH=$PYTHONPATH:$PATH
export PYTHONPATH

projectdir=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MIRCen_sammba-mri_mouse
savedir=$projectdir/analysis20171012
mkdir $savedir

dicomdirlist=$projectdir/dicomdirs.txt
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
