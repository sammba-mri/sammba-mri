#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE

PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/common
PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/projects/XNAT_CSD_MRI_MOUSE
export PATH

PYTHONPATH=$PYTHONPATH:$PATH
export PYTHONPATH

projectdir=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/XNAT_CSD_MRI_MOUSE
savedir=$projectdir/analysis20171018
tsavedir1=$savedir/goodbaselinetemplate; mkdir -p $tsavedir1
tsavedir2=$savedir/allgoodtemplate; mkdir -p $tsavedir2
tsavedir3=$savedir/allbaselinetemplate; mkdir -p $tsavedir3
tsavedir4=$savedir/goodposttemplate; mkdir -p $tsavedir4

conv=0.005
twoblur=1
brainvol=400
#-rbt values might be improveable
Urad=18.3
b=70
t=80

MIC40C57=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MRIatlases/MIC40C57/20170223
dorrbrain=$MIC40C57/brain100.nii.gz
dorratlas=$MIC40C57/labels100.nii.gz

pipeline () {

	pdldirlist=$1 psavedir=$2 templatedir=$3 dorratlas=$4

	nlines=$(wc -l $pdldirlist | awk '{print $1}')

	for ((a=1; a<=$nlines; a++)); do

		dldir=$(head -n$a < $pdldirlist | tail -n1)
		analdir=$psavedir/$(basename $dldir)
		mkdir -p $analdir
		3dcopy $dldir/3DRARE.nii $analdir/anat_n0.nii.gz
		3dcopy $dldir/rsfMRI.nii $analdir/rs_n0.nii.gz
		3drefit -xyzscale 0.1 $analdir/rs_n0.nii.gz
		3dCM -set 0 0 0 $analdir/rs_n0.nii.gz

		#orientation correctors for use later
		#https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Euler_angles_.28_z-y.E2.80.99-x.E2.80.B3_intrinsic.29_.E2.86.92_Rotation_matrix
		echo 1 0 0 0 0 0 1 0 0 -1 0 0 > $analdir/x270.1d
		echo -1 0 0 0 0 1 0 0 0 0 -1 0 > $analdir/y180.1d

		for niiname in anat_n0; do
			nifti_tool -disp_hdr -field srow_x -field srow_y -field srow_z -infiles $analdir/$niiname.nii.gz -quiet > $analdir/sform.txt
			fsl5.0-fslorient -setsform $(cat_matvec -ONELINE $analdir/sform.txt $analdir/x270.1d $analdir/y180.1d) 0 0 0 1 $analdir/$niiname.nii.gz
			fsl5.0-fslorient -copysform2qform $analdir/$niiname.nii.gz
		done

		for niiname in rs_n0; do
			nifti_tool -disp_hdr -field srow_x -field srow_y -field srow_z -infiles $analdir/$niiname.nii.gz -quiet > $analdir/sform.txt
			fsl5.0-fslorient -setsform $(cat_matvec -ONELINE $analdir/sform.txt $analdir/x270.1d $analdir/y180.1d) 0 0 0 1 $analdir/$niiname.nii.gz
			fsl5.0-fslorient -copysform2qform $analdir/$niiname.nii.gz
			3dLRflip -prefix $analdir/$niiname.nii.gz $analdir/$niiname.nii.gz #not 100% sure about this
		done

		rm -f $analdir/sform.txt $analdir/*.1d
		
	done
	
	pipeline=pipeline_noconversion.bash #1
	mkdir -p $psavedir
	brain=$templatedir/Qw4_meanhead_Na.nii.gz #4
	atlas=$dorratlas #5
	mask=$dorratlas #6
	head=$templatedir/Qw4_meanhead_Na.nii.gz #7
	headweight=$templatedir/aff3_unionmaskdil4_Na.nii.gz #8
	basetype=head #9
	dofolderoverwrite=no #10
	tmpdir=/volatile #11
	registerfunctional=yes #12
	subpipeline=perslice_registration_subpipeline.bash #13
	Urad=18.3 #14
	brainvol=400 #15
	scale=0.1 #16

	cd $psavedir

	rawimdir=NA

	for dir in $(find $psavedir -mindepth 1 -maxdepth 1 -type d); do
		NIfTIdir=$(readlink -e $dir)
		pipeline_noconversion.bash $rawimdir $NIfTIdir $brain $atlas $mask $head $headweight $basetype $tmpdir $registerfunctional $subpipeline $Urad $brainvol $scale
	done

}

pipeline $projectdir/allbaseline.txt $savedir/baseline $tsavedir4 $dorratlas
pipeline $projectdir/allpost.txt $savedir/post $tsavedir4 $dorratlas
