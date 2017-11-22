#!/bin/bash

#bash -xef /home/nadkarni/git/sammba-mri/sammba/projects/DZNE/pilot_DZNE.bash 2>&1 | tee /home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/DZNE/20170804/$(date +%Y%m%d_%H%M%S).log

startdir=$(pwd)

PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/common
export PATH

PYTHONPATH=$PYTHONPATH:$PATH
export PYTHONPATH

projectdir=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/DZNE/20170804
python -m PVEnDCMtoNIfTI.py /usr/bin/dcmdump $projectdir/rawdata $projectdir/processed yes || true

cd $projectdir/processed

mkdir 21
mv 21__20170317__085458__T1_FLASH_3D_iso__fixedSIAP__bf20170804.nii.gz 21/anat_n0.nii.gz
mv 21__20170317__085733__T1_FLASH__fixedSIAP__bf20170804.nii.gz 21/rs_n0.nii.gz
mkdir 8619
mv 8619__20170317__104039__T1_FLASH_3D_iso__fixedSIAP__bf20170804.nii.gz 8619/anat_n0.nii.gz
mv 8619__20170317__101702__T1_FLASH__fixedSIAP__bf20170804.nii.gz 8619/rs_n0.nii.gz 

cd $startdir

MIC40C57=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MRIatlases/MIC40C57/20170223
brain=$MIC40C57/brain100.nii.gz
atlas=$MIC40C57/labels100.nii.gz
mask=$MIC40C57/mask100.nii.gz
head=$MIC40C57/head100.nii.gz
3dmask_tool -dilate_inputs 7 -inputs $mask -prefix $projectdir/processed/maskdil7.nii.gz
headweight=$projectdir/processed/maskdil7.nii.gz
basetype=brain
Urad=18.3
brainvol=400
scale=0.1
biascorrector=3dUnifize
tmpdir=/volatile
registerfunctional=no
rminterfiles=yes
subpipeline=perslice_registration_subpipeline.bash

for animal in 21 8619; do
	cd $projectdir/processed/$animal
	anattotemplate.bash . yes $biascorrector $brainvol $scale $brain $atlas $mask $head $headweight $basetype $Urad
	3dTcat -prefix expt.nii.gz rs_n0.nii.gz rs_n0.nii.gz rs_n0.nii.gz rs_n0.nii.gz rs_n0.nii.gz rs_n0.nii.gz rs_n0.nii.gz rs_n0.nii.gz rs_n0.nii.gz rs_n0.nii.gz rs_n0.nii.gz rs_n0.nii.gz rs_n0.nii.gz rs_n0.nii.gz rs_n0.nii.gz
	rm -f rs_n0.nii.gz
	mv expt.nii.gz rs_n0.nii.gz
	rs.bash . $subpipeline $brain $tmpdir $registerfunctional $rminterfiles $brainvol
done

cd $startdir
