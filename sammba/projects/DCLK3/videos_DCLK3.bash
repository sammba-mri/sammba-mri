#!/bin/bash

#bash /home/nadkarni/git/sammba-mri/sammba/projects/DCLK3/videos_DCLK3.bash

export AFNI_DECONFLICT=OVERWRITE #saves a lot of grief

cd /home/Pplateforme/Plate-forme_RMN/NachiketNadkarni/DCLK3/processed_20170727/MRIsessions

#raw
fsl5.0-fslroi raw_video.nii.gz raw_video_x83.nii.gz 83 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim raw_video_x83.nii.gz -z y x raw_video_x83.nii.gz
gunzip raw_video_x83.nii.gz

fsl5.0-fslroi raw_video.nii.gz raw_video_y50.nii.gz 0 -1 50 1 0 -1 0 -1
fsl5.0-fslswapdim raw_video_y50.nii.gz x -z y raw_video_y50.nii.gz
gunzip raw_video_y50.nii.gz

fsl5.0-fslroi raw_video.nii.gz raw_video_z98.nii.gz 0 -1 0 -1 98 1 0 -1
gunzip raw_video_z98.nii.gz

#Qw4
fsl5.0-fslroi Qw4_videohead.nii.gz Qw4_videohead_x79.nii.gz 79 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim Qw4_videohead_x79.nii.gz -z y x Qw4_videohead_x79.nii.gz
gunzip Qw4_videohead_x79.nii.gz

fsl5.0-fslroi Qw4_videohead.nii.gz Qw4_videohead_y87.nii.gz 0 -1 87 1 0 -1 0 -1
fsl5.0-fslswapdim Qw4_videohead_y87.nii.gz x -z y Qw4_videohead_y87.nii.gz
gunzip Qw4_videohead_y87.nii.gz

fsl5.0-fslroi Qw4_videohead.nii.gz Qw4_videohead_z110.nii.gz 0 -1 0 -1 110 1 0 -1
gunzip Qw4_videohead_z110.nii.gz
