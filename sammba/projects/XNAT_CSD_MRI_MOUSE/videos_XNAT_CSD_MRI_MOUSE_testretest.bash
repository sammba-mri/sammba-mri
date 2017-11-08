#!/bin/bash

#bash /home/nadkarni/git/sammba-mri/sammba/projects/XNAT_CSD_MRI_MOUSE_testretest/videos_XNAT_CSD_MRI_MOUSE_testretest.bash

export AFNI_DECONFLICT=OVERWRITE #saves a lot of grief

cd /home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/XNAT_CSD_MRI_MOUSE_testretest/analysis20170926/baseline

3dTcat -prefix raw_anat_video.nii.gz $(find -name 'anat_n0.nii.gz')
3dTcat -prefix Qw_anat_video.nii.gz $(find -name 'anat_n0_UnAaQw.nii.gz')
3dTcat -prefix raw_func_video.nii.gz $(find -name 'rs_n0_TsAvAv.nii.gz')
3dTcat -prefix Qw_func_video.nii.gz $(find -name 'rs_n0_TsAvAvN3_NaNaMe.nii.gz')

#rawanat
fsl5.0-fslroi raw_anat_video.nii.gz raw_anat_video_x75.nii.gz 75 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim raw_anat_video_x75.nii.gz -z y x raw_anat_video_x75.nii.gz
gunzip raw_anat_video_x75.nii.gz
fsl5.0-fslroi raw_anat_video.nii.gz raw_anat_video_y85.nii.gz 0 -1 85 1 0 -1 0 -1
fsl5.0-fslswapdim raw_anat_video_y85.nii.gz -x z y raw_anat_video_y85.nii.gz
gunzip raw_anat_video_y85.nii.gz
fsl5.0-fslroi raw_anat_video.nii.gz raw_anat_video_z40.nii.gz 0 -1 0 -1 40 1 0 -1
gunzip raw_anat_video_z40.nii.gz

#rawfunc
fsl5.0-fslroi raw_func_video.nii.gz raw_func_video_x44.nii.gz 44 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim raw_func_video_x44.nii.gz -z y x raw_func_video_x44.nii.gz
gunzip raw_func_video_x44.nii.gz
fsl5.0-fslroi raw_func_video.nii.gz raw_func_video_y51.nii.gz 0 -1 51 1 0 -1 0 -1
fsl5.0-fslswapdim raw_func_video_y51.nii.gz -x -z -y raw_func_video_y51.nii.gz
gunzip raw_func_video_y51.nii.gz
fsl5.0-fslroi raw_func_video.nii.gz raw_func_video_z5.nii.gz 0 -1 0 -1 5 1 0 -1
gunzip raw_func_video_z5.nii.gz

#Qwanat
fsl5.0-fslroi Qw_anat_video.nii.gz Qw_anat_video_x75.nii.gz 75 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim Qw_anat_video_x75.nii.gz -z y x Qw_anat_video_x75.nii.gz
gunzip Qw_anat_video_x75.nii.gz
fsl5.0-fslroi Qw_anat_video.nii.gz Qw_anat_video_y105.nii.gz 0 -1 105 1 0 -1 0 -1
fsl5.0-fslswapdim Qw_anat_video_y105.nii.gz -x z y Qw_anat_video_y105.nii.gz
gunzip Qw_anat_video_y105.nii.gz
fsl5.0-fslroi Qw_anat_video.nii.gz Qw_anat_video_z40.nii.gz 0 -1 0 -1 40 1 0 -1
gunzip Qw_anat_video_z40.nii.gz

#Qwfunc
fsl5.0-fslroi Qw_func_video.nii.gz Qw_func_video_x75.nii.gz 75 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim Qw_func_video_x75.nii.gz -z y x Qw_func_video_x75.nii.gz
gunzip Qw_func_video_x75.nii.gz
fsl5.0-fslroi Qw_func_video.nii.gz Qw_func_video_y105.nii.gz 0 -1 105 1 0 -1 0 -1
fsl5.0-fslswapdim Qw_func_video_y105.nii.gz -x z y Qw_func_video_y105.nii.gz
gunzip Qw_func_video_y105.nii.gz
fsl5.0-fslroi Qw_func_video.nii.gz Qw_func_video_z40.nii.gz 0 -1 0 -1 40 1 0 -1
gunzip Qw_func_video_z40.nii.gz
