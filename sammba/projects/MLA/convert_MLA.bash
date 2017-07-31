#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE

rawdatadir=$(readlink -e $1)
savedir=$(readlink -e $2)
resampleres=$3

echo 1 0 0 0 0 0 1 0 0 -1 0 0 > $savedir/x90.1d
echo -1 0 0 0 0 1 0 0 0 0 -1 0 > $savedir/y180.1d

find -L $rawdatadir -type d -name '*T2Lemur*.img' | grep -v Mc147BCBB-01-s_2011042110/ptk_server_task_1 | sort > $savedir/rawdatadirs.txt
nlines=$(wc -l $savedir/rawdatadirs.txt | awk '{print $1}')

for ((a=1; a<=$nlines; a++)); do

	fdfdir=$(head -n$a < $savedir/rawdatadirs.txt | tail -n1)
	newnifti=$(echo $fdfdir | sed "s&$rawdatadir/&&g" | sed 's/\//__/g' | sed 's/.img//g' | tr -cd [:alnum:]-_)
	newniftidir=$savedir/$newnifti
	mkdir $newniftidir

	AimsFileConvert -i "$fdfdir"/slice001image001echo001.fdf -o $newniftidir/$newnifti.nii.gz -r 1 --omin 0 --omax 1000

	fsl5.0-fslswapdim $newniftidir/$newnifti.nii.gz -x -y z $newniftidir/anat.nii.gz

	nifti_tool -disp_hdr -field srow_x -field srow_y -field srow_z -infiles $newniftidir/anat.nii.gz -quiet > $newniftidir/sform.txt

	fsl5.0-fslorient -setsform $(cat_matvec -ONELINE $newniftidir/sform.txt $savedir/x90.1d $savedir/y180.1d) 0 0 0 1 $newniftidir/anat.nii.gz
	fsl5.0-fslorient -copysform2qform $newniftidir/anat.nii.gz
	3dresample -dxyz $resampleres $resampleres $resampleres -rmode Cubic -prefix $newniftidir/anat.nii.gz -inset $newniftidir/anat.nii.gz
	3dZeropad -RL 256 -AP 256 -IS 256 -prefix $newniftidir/anat.nii.gz $newniftidir/anat.nii.gz
	
done

3dTcat -prefix $savedir/raw_video.nii.gz $(find $savedir -name 'anat.nii.gz')
3dTstat -prefix $savedir/raw_mean.nii.gz $savedir/raw_video.nii.gz
