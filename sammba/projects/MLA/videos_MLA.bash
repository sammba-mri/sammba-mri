#!/bin/bash

#bash /home/Pmamobipet/Tvx-Manips-MD_/MD_1701-Microcebe-Creation-Atlas/code4/videos.bash

export AFNI_DECONFLICT=OVERWRITE #saves a lot of grief

cd /home/Pmamobipet/Tvx-Manips-MD_/MD_1701-Microcebe-Creation-Atlas/processed_20170411

#raw
fsl5.0-fslroi raw_video.nii.gz raw_video_x128.nii.gz 128 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim raw_video_x128.nii.gz -z y x raw_video_x128.nii.gz
gunzip raw_video_x128.nii.gz

fsl5.0-fslroi raw_video.nii.gz raw_video_y80.nii.gz 0 -1 80 1 0 -1 0 -1
fsl5.0-fslswapdim raw_video_y80.nii.gz x -z y raw_video_y80.nii.gz
gunzip raw_video_y80.nii.gz

fsl5.0-fslroi raw_video.nii.gz raw_video_z140.nii.gz 0 -1 0 -1 140 1 0 -1
gunzip raw_video_z140.nii.gz

fsl5.0-fslroi raw_video.nii.gz raw_video_z128.nii.gz 0 -1 0 -1 128 1 0 -1
gunzip raw_video_z128.nii.gz

#UnBmBe
fsl5.0-fslroi UnBmBe_video.nii.gz UnBmBe_video_z140.nii.gz 0 -1 0 -1 140 1 0 -1
gunzip UnBmBe_video_z140.nii.gz

#UnBmBeCC
fsl5.0-fslroi UnBmBeCC_video.nii.gz UnBmBeCC_video_z140.nii.gz 0 -1 0 -1 140 1 0 -1
gunzip UnBmBeCC_video_z140.nii.gz

#shr1
fsl5.0-fslroi shr1_video.nii.gz shr1_video_x128.nii.gz 128 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim shr1_video_x128.nii.gz -z y x shr1_video_x128.nii.gz
gunzip shr1_video_x128.nii.gz

fsl5.0-fslroi shr1_video.nii.gz shr1_video_y127.nii.gz 0 -1 127 1 0 -1 0 -1
fsl5.0-fslswapdim shr1_video_y127.nii.gz x -z y shr1_video_y127.nii.gz
gunzip shr1_video_y127.nii.gz

fsl5.0-fslroi shr1_video.nii.gz shr1_video_z140.nii.gz 0 -1 0 -1 140 1 0 -1
gunzip shr1_video_z140.nii.gz

#shr2
fsl5.0-fslroi shr2_video.nii.gz shr2_video_x128.nii.gz 128 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim shr2_video_x128.nii.gz -z y x shr2_video_x128.nii.gz
gunzip shr2_video_x128.nii.gz

fsl5.0-fslroi shr2_video.nii.gz shr2_video_y127.nii.gz 0 -1 127 1 0 -1 0 -1
fsl5.0-fslswapdim shr2_video_y127.nii.gz x -z y shr2_video_y127.nii.gz
gunzip shr2_video_y127.nii.gz

fsl5.0-fslroi shr2_video.nii.gz shr2_video_z140.nii.gz 0 -1 0 -1 140 1 0 -1
gunzip shr2_video_z140.nii.gz

fsl5.0-fslroi shr2_videohead.nii.gz shr2_videohead_z140.nii.gz 0 -1 0 -1 140 1 0 -1
gunzip shr2_videohead_z140.nii.gz

#aff3
fsl5.0-fslroi aff3_video.nii.gz aff3_video_x128.nii.gz 128 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim aff3_video_x128.nii.gz -z y x aff3_video_x128.nii.gz
gunzip aff3_video_x128.nii.gz

fsl5.0-fslroi aff3_video.nii.gz aff3_video_y127.nii.gz 0 -1 127 1 0 -1 0 -1
fsl5.0-fslswapdim aff3_video_y127.nii.gz x -z y aff3_video_y127.nii.gz
gunzip aff3_video_y127.nii.gz

fsl5.0-fslroi aff3_video.nii.gz aff3_video_z140.nii.gz 0 -1 0 -1 140 1 0 -1
gunzip aff3_video_z140.nii.gz

fsl5.0-fslroi aff3_videohead.nii.gz aff3_videohead_z140.nii.gz 0 -1 0 -1 140 1 0 -1
gunzip aff3_videohead_z140.nii.gz

#UnCC
fsl5.0-fslroi UnCC_video.nii.gz UnCC_video_z140.nii.gz 0 -1 0 -1 140 1 0 -1
gunzip UnCC_video_z140.nii.gz

#Qw1
fsl5.0-fslroi Qw1_videohead.nii.gz Qw1_videohead_z140.nii.gz 0 -1 0 -1 140 1 0 -1
gunzip Qw1_videohead_z140.nii.gz

#Qw2
fsl5.0-fslroi Qw2_videohead.nii.gz Qw2_videohead_z140.nii.gz 0 -1 0 -1 140 1 0 -1
gunzip Qw2_videohead_z140.nii.gz

#Qw3
fsl5.0-fslroi Qw3_videohead.nii.gz Qw3_videohead_z140.nii.gz 0 -1 0 -1 140 1 0 -1
gunzip Qw3_videohead_z140.nii.gz

#Qw4
fsl5.0-fslroi Qw4_videohead.nii.gz Qw4_videohead_x128.nii.gz 128 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim Qw4_videohead_x128.nii.gz -z y x Qw4_videohead_x128.nii.gz
gunzip Qw4_videohead_x128.nii.gz

fsl5.0-fslroi Qw4_videohead.nii.gz Qw4_videohead_y80.nii.gz 0 -1 80 1 0 -1 0 -1
fsl5.0-fslswapdim Qw4_videohead_y80.nii.gz x -z y Qw4_videohead_y80.nii.gz
gunzip Qw4_videohead_y80.nii.gz

fsl5.0-fslroi Qw4_videohead.nii.gz Qw4_videohead_y125.nii.gz 0 -1 125 1 0 -1 0 -1
fsl5.0-fslswapdim Qw4_videohead_y125.nii.gz x -z y Qw4_videohead_y125.nii.gz
gunzip Qw4_videohead_y125.nii.gz

fsl5.0-fslroi Qw4_videohead.nii.gz Qw4_videohead_z140.nii.gz 0 -1 0 -1 140 1 0 -1
gunzip Qw4_videohead_z140.nii.gz
