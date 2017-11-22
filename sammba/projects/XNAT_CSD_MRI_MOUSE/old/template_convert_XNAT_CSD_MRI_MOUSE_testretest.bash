#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE

savedir=$(readlink -e $1)
anatlist=$(readlink -e $2)

#orientation correctors for use later
#https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Euler_angles_.28_z-y.E2.80.99-x.E2.80.B3_intrinsic.29_.E2.86.92_Rotation_matrix
echo 1 0 0 0 0 0 1 0 0 -1 0 0 > $savedir/x270.1d
echo -1 0 0 0 0 1 0 0 0 0 -1 0 > $savedir/y180.1d

nlines=$(wc -l $anatlist | awk '{print $1}')

for ((a=1; a<=$nlines; a++)); do

	rawdir=$(head -n$a < $anatlist | tail -n1)
	stem=$(basename $rawdir)
	prefix=$(basename $(dirname $rawdir))
	newniftidir=$savedir/"$prefix"_"$stem"
	mkdir -p $newniftidir
	cp $rawdir/3DRARE.nii.gz $newniftidir/anat.nii.gz
	nifti_tool -disp_hdr -field srow_x -field srow_y -field srow_z -infiles $newniftidir/anat.nii.gz -quiet > $newniftidir/sform.txt
	fsl5.0-fslorient -setsform $(cat_matvec -ONELINE $newniftidir/sform.txt $savedir/x270.1d $savedir/y180.1d) 0 0 0 1 $newniftidir/anat.nii.gz
	fsl5.0-fslorient -copysform2qform $newniftidir/anat.nii.gz
	
done

3dTcat -prefix $savedir/raw_video.nii.gz $(find $savedir -name anat.nii.gz)
3dTstat -prefix $savedir/raw_mean.nii.gz $savedir/raw_video.nii.gz
