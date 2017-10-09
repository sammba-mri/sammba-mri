#!/bin/bash

#bash -xef /home/nadkarni/git/sammba-mri/sammba/projects/MIRCen_sammba-mri_mouselemur/master_MIRCen_sammba-mri_mouselemur.bash 2>&1 | tee /home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouselemur/MIRCen_sammba-mri_mouselemur/$(date +%Y%m%d_%H%M%S).log

PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/common
export PATH

PYTHONPATH=$PYTHONPATH:$PATH
export PYTHONPATH

pipeline=pipeline_conversiononly.bash #1
projectdir=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouselemur/MIRCen_sammba-mri_mouselemur
MRIsessionslist=$projectdir/MRISessions.txt #2
analysisdir=$projectdir/analysis20170928
mkdir $analysisdir #3
templatedir=$projectdir/20170928_template
brain=$templatedir/MLAbrainAlQw.nii.gz #4
atlas=$templatedir/MLAatlasAlQw.nii.gz #5
mask=$templatedir/MLAatlasAlQw.nii.gz #6
head=$templatedir/Qw4_meanhead.nii.gz #7
headweight=$templatedir/aff3_unionmaskdil3.nii.gz #8
basetype=head #9
dofolderoverwrite=no #10
tmpdir=/volatile #11
registerfunctional=no #12
subpipeline=perslice_registration_subpipeline.bash #13
Urad=18.3 #14
brainvol=1600 #15
scale=0.2 #16
T1blood=2800 #17
lambda=0.9 #18
multiplier=6000000 #19
T1guess=1600 #20
mccores=16 #21

#					  1			2				 3			  4		 5		6	  7		8			9		  10				 11		 12					 13			  14	15		  16
pipeline_wrapper.bash $pipeline $MRIsessionslist $analysisdir $brain $atlas $mask $head $headweight $basetype $dofolderoverwrite $tmpdir $registerfunctional $subpipeline $Urad $brainvol $scale

cd $analysisdir

mv 20170425_115121_208CBF_1_1/anat_n0.nii.gz 20170425_115121_208CBF_1_1/badanat_n1.nii.gz
mv 20170425_115121_208CBF_1_1/208CBF__20170425__132441__MSME_MIRCen_RatBrain_atlas__fixedSIAP__bf11.nii.gz 20170425_115121_208CBF_1_1/anat_n0.nii.gz
mv 20170425_164314_289BB_1_1/anat_n0.nii.gz 20170425_164314_289BB_1_1/badanat_n2.nii.gz
mv 20170425_164314_289BB_1_1/289BB__20170425__171045__MSME_200um__fixedSIAP__bf6.nii.gz 20170425_164314_289BB_1_1/anat_n0.nii.gz
mv 20170425_192223_310C_1_1/anat_n0.nii.gz 20170425_192223_310C_1_1/badanat_n1.nii.gz
mv 20170425_192223_310C_1_1/310C__20170425__194739__MSME_200um__fixedSIAP__bf6.nii.gz 20170425_192223_310C_1_1/anat_n0.nii.gz

rawimdir=NA

for dir in $(find $analysisdir -mindepth 1 -maxdepth 1 -type d); do
	NIfTIdir=$(readlink -e $dir)
	pipeline_noconversion.bash $rawimdir $NIfTIdir $brain $atlas $mask $head $headweight $basetype $tmpdir $registerfunctional $subpipeline $Urad $brainvol $scale
done
