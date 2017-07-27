#!/bin/bash

#https://afni.nimh.nih.gov/afni/community/board/read.php?1,149385,149385#msg-149385
#https://afni.nimh.nih.gov/afni/community/board/read.php?1,74225,74225#msg-74225

export AFNI_DECONFLICT=OVERWRITE

startdir=$(pwd)

workdir=/tmp/tmp_fixobliquity_"$RANDOM"
mkdir $workdir

goodoblique=$1
badoblique=$2
goodobliquename=$(basename -s .nii.gz $1)
badobliquename=$(basename -s .nii.gz $2)

3dcopy $goodoblique $workdir/"$goodobliquename"+orig
3dcopy $badoblique $workdir/"$badobliquename"+orig
3drefit -atrcopy $workdir/"$goodobliquename"+orig IJK_TO_DICOM_REAL $workdir/"$badobliquename"+orig
3dcopy $workdir/"$badobliquename"+orig $badoblique

rm -rfd $workdir