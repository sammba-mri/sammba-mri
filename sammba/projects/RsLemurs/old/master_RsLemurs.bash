#!/bin/bash

#bash -xef /home/nadkarni/git/sammba-mri/sammba/projects/RsLemurs/master_RsLemurs.bash 2>&1 | tee /home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouselemur/$(date +%Y%m%d_%H%M%S).log

PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/common
export PATH

PYTHONPATH=$PYTHONPATH:$PATH
export PYTHONPATH

pipeline=pipeline_conversiononly.bash #1
MRIsessionslist=/home/Pmamobipet/Tvx-Manips-MD_/MD_1704_RsLemurs/analyses/acqlist.txt #2
analysisdir=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouselemur/analysis20170825 #3
mkdir $analysisdir
MLA=/home/Pmamobipet/Tvx-Manips-MD_/MD_1701-Microcebe-Creation-Atlas/Final-Images-For-Atlas/20170711
brain=$MLA/brain_200.nii.gz #4
atlas=$MLA/labels_200.nii.gz #5
mask=$MLA/mask_200.nii.gz #6
head=$MLA/head_200.nii.gz #7
headweight=$MLA/aff3_unionmaskdil7_200.nii.gz #8
basetype=head #9
dofolderoverwrite=no #10
tmpdir=/volatile #11
registerfunctional=no #12
subpipeline=perslice_registration_subpipeline.bash #13
Urad=18.3 #14
brainvol=1500 #15
scale=0.2 #16

#					  1			2				 3			  4		 5		6	  7		8			9		  10				 11		 12					 13			  14	15		  16
pipeline_wrapper.bash $pipeline $MRIsessionslist $analysisdir $brain $atlas $mask $head $headweight $basetype $dofolderoverwrite $tmpdir $registerfunctional $subpipeline $Urad $brainvol $scale

cd $analysisdir

mv 20170425_115121_208CBF_1_1/anat_n0.nii.gz 20170425_115121_208CBF_1_1/badanat_n1.nii.gz
mv 20170425_115121_208CBF_1_1/208CBF__20170425__132441__MSME_MIRCen_RatBrain_atlas__fixedSIAP__bf11.nii.gz 20170425_115121_208CBF_1_1/anat_n0.nii.gz
mv 20170425_164314_289BB_1_1/anat_n0.nii.gz 20170425_164314_289BB_1_1/badanat_n2.nii.gz
mv 20170425_164314_289BB_1_1/289BB__20170425__171045__MSME_200um__fixedSIAP__bf6.nii.gz 20170425_164314_289BB_1_1/anat_n0.nii.gz
mv 20170425_192223_310C_1_1/anat_n0.nii.gz 20170425_192223_310C_1_1/badanat_n1.nii.gz
mv 20170425_192223_310C_1_1/310C__20170425__194739__MSME_200um__fixedSIAP__bf6.nii.gz 20170425_192223_310C_1_1/anat_n0.nii.gz
mv 20170426_175644_MD1704_Mc276BC_P01_1_1/anat_n0.nii.gz 20170426_175644_MD1704_Mc276BC_P01_1_1/badanat_n2.nii.gz
mv 20170426_175644_MD1704_Mc276BC_P01_1_1/MD1704-Mc276BC-P01__20170426__182213__MSME_200um__fixedSIAP__bf8.nii.gz 20170426_175644_MD1704_Mc276BC_P01_1_1/anat_n0.nii.gz
mv 20170427_141151_MD1704_Mc283CCA_P01_1_1/anat_n0.nii.gz 20170427_141151_MD1704_Mc283CCA_P01_1_1/badanat_n1.nii.gz
mv 20170427_141151_MD1704_Mc283CCA_P01_1_1/MD1704-Mc283CCA-P01__20170427__143852__MSME_200um_ZF__fixedSIAP__bf10.nii.gz 20170427_141151_MD1704_Mc283CCA_P01_1_1/anat_n0.nii.gz
mv 20170427_155021_MD1704_Mc283EA_P01_1_1/anat_n0.nii.gz 20170427_155021_MD1704_Mc283EA_P01_1_1/badanat_n1.nii.gz
mv 20170427_155021_MD1704_Mc283EA_P01_1_1/MD1704-Mc283EA-P01__20170427__161214__MSME_200um__fixedSIAP__bf5.nii.gz 20170427_155021_MD1704_Mc283EA_P01_1_1/anat_n0.nii.gz
mv 20170427_171736_MD1704_Mc263BCE_P01_1_1/anat_n0.nii.gz 20170427_171736_MD1704_Mc263BCE_P01_1_1/badanat_n1.nii.gz
mv 20170427_171736_MD1704_Mc263BCE_P01_1_1/MD1704-Mc263BCE-P01__20170427__174410__MSME_200um_ZF__fixedSIAP__bf6.nii.gz 20170427_171736_MD1704_Mc263BCE_P01_1_1/anat_n0.nii.gz
mv 20170427_185244_MD1704_Mc285D_P01_1_1/anat_n0.nii.gz 20170427_185244_MD1704_Mc285D_P01_1_1/badanat_n1.nii.gz
mv 20170427_185244_MD1704_Mc285D_P01_1_1/MD1704-Mc285D-P01__20170427__191707__MSME_200um_ZF__fixedSIAP__bf6.nii.gz 20170427_185244_MD1704_Mc285D_P01_1_1/anat_n0.nii.gz

rawimdir=NA

for dir in $(find $analysisdir -mindepth 1 -maxdepth 1 -type d); do
	NIfTIdir=$(readlink -e $dir)
	pipeline_noconversion.bash $rawimdir $NIfTIdir $brain $atlas $mask $head $headweight $basetype $tmpdir $registerfunctional $subpipeline $Urad $brainvol $scale
done
