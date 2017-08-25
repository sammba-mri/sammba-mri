#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE

dldir=$(readlink -e $1)
NIfTIdir=$(readlink -e $2)

echo "1 dldir" $dldir
echo "2 NIfTIdir" $NIfTIdir

3dcopy $dldir/3DRARE.nii $NIfTIdir/anat_n0.nii.gz
3dcopy $dldir/rsfMRI.nii $NIfTIdir/rs_n0.nii.gz
3drefit -xyzscale 0.1 $NIfTIdir/rs_n0.nii.gz
3dCM -set 0 0 0 $NIfTIdir/rs_n0.nii.gz

#orientation correctors for use later
#https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Euler_angles_.28_z-y.E2.80.99-x.E2.80.B3_intrinsic.29_.E2.86.92_Rotation_matrix
echo 1 0 0 0 0 0 1 0 0 -1 0 0 > $NIfTIdir/x270.1d
echo -1 0 0 0 0 1 0 0 0 0 -1 0 > $NIfTIdir/y180.1d

for niiname in anat_n0; do
	nifti_tool -disp_hdr -field srow_x -field srow_y -field srow_z -infiles $NIfTIdir/$niiname.nii.gz -quiet > $NIfTIdir/sform.txt
	fsl5.0-fslorient -setsform $(cat_matvec -ONELINE $NIfTIdir/sform.txt $NIfTIdir/x270.1d $NIfTIdir/y180.1d) 0 0 0 1 $NIfTIdir/$niiname.nii.gz
	fsl5.0-fslorient -copysform2qform $NIfTIdir/$niiname.nii.gz
done

for niiname in rs_n0; do
	nifti_tool -disp_hdr -field srow_x -field srow_y -field srow_z -infiles $NIfTIdir/$niiname.nii.gz -quiet > $NIfTIdir/sform.txt
	fsl5.0-fslorient -setsform $(cat_matvec -ONELINE $NIfTIdir/sform.txt $NIfTIdir/x270.1d $NIfTIdir/y180.1d) 0 0 0 1 $NIfTIdir/$niiname.nii.gz
	fsl5.0-fslorient -copysform2qform $NIfTIdir/$niiname.nii.gz
	3dLRflip -prefix $NIfTIdir/$niiname.nii.gz $NIfTIdir/$niiname.nii.gz #not 100% sure about this
done

rm -f $NIfTIdir/sform.txt $NIfTIdir/*.1d
