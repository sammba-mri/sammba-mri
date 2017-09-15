#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE

rawdatadir=$(readlink -e $1)
savedir=$(readlink -e $2)

find -L $rawdatadir -type d -name '*T2Lemur*.fid' | grep -v Mc147BCBB-01-s_2011042110/ptk_server_task_1 | sort > $savedir/rawdatadirs.txt
nlines=$(wc -l $savedir/rawdatadirs.txt | awk '{print $1}')

for ((a=1; a<=$nlines; a++)); do

	fiddir=$(head -n$a < $savedir/rawdatadirs.txt | tail -n1)
	newnifti=$(echo $fiddir | sed "s&$rawdatadir/&&g" | sed 's/\//__/g' | sed 's/.img//g' | tr -cd [:alnum:]-_)
	newniftidir=$savedir/$newnifti
	mkdir $newniftidir
	
	python -m VarianAgilentFIDtoNIfTI.py $fiddir $newniftidir/anat.nii.gz	
	
done

3dTcat -prefix $savedir/raw_video.nii.gz $(find $savedir -name 'anat.nii.gz')
3dTstat -prefix $savedir/raw_mean.nii.gz $savedir/raw_video.nii.gz
