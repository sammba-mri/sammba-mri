#!/bin/bash

#bash -xef /home/nadkarni/git/sammba-mri/sammba/projects/AmylNet/master_AmylNet.bash 2>&1 | tee /home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/AmylNet/$(date +%Y%m%d_%H%M%S).log

PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/spatialregistration
export PATH

PYTHONPATH=$PYTHONPATH:$PATH
export PYTHONPATH

pipeline=pipeline_all.bash #1
MRIsessionslist=/home/Pmamobipet/Tvx-Manips-MD_/MD_1602-AmylNet-Garin/MRIanalyses/texteAmylNet/MRISessions6.txt #2
analysisdir=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/AmylNet/analysis20170725/MRIsessions #3
MIC40C57=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MRIatlases/MIC40C57/20170223
brain=$MIC40C57/brain100.nii.gz #4
atlas=$MIC40C57/labels100.nii.gz #5
mask=$MIC40C57/mask100.nii.gz #6
head=$MIC40C57/head100.nii.gz #7
3dmask_tool -dilate_inputs 7 -inputs $mask -prefix $analysisdir/maskdil7.nii.gz
headweight=$analysisdir/maskdil7.nii.gz #8
basetype=brain #9
dofolderoverwrite=no #10
tmpdir=/volatile #11
registerfunctional=no #12
subpipeline=perslice_registration_subpipeline.bash #13
Urad=18.3 #14
brainvol=400 #15
scale=0.1 #16

#					  1			2				 3			  4		 5		6	  7		8			9		  10				 11		 12					 13			  14	15		  16
pipeline_wrapper.bash $pipeline $MRIsessionslist $analysisdir $brain $atlas $mask $head $headweight $basetype $dofolderoverwrite $tmpdir $registerfunctional $subpipeline $Urad $brainvol $scale
