#!/bin/bash

#bash /home/nadkarni/git/sammba-mri/sammba/projects/MIRCen_sammba-mri_mouse/videos_MIRCen_sammba-mri_mouse.bash

export AFNI_DECONFLICT=OVERWRITE #saves a lot of grief

cd /home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MIRCen_sammba-mri_mouse/analysis20170926

3dTcat -prefix raw_anat_video.nii.gz $(find -name 'anat_n0.nii.gz')
3dTcat -prefix Qw_anat_video.nii.gz $(find -name 'anat_n0_UnAaQw.nii.gz')
3dTcat -prefix raw_func_video.nii.gz $(find -name 'rs_n0_TsAvAv.nii.gz')
3dTcat -prefix Qw_func_video.nii.gz $(find -name 'rs_n0_TsAvAvN3_NaNaMe.nii.gz')

#rawanat
fsl5.0-fslroi raw_anat_video.nii.gz raw_anat_video_x46.nii.gz 46 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim raw_anat_video_x46.nii.gz -z y x raw_anat_video_x46.nii.gz
gunzip raw_anat_video_x46.nii.gz
fsl5.0-fslroi raw_anat_video.nii.gz raw_anat_video_y98.nii.gz 0 -1 98 1 0 -1 0 -1
fsl5.0-fslswapdim raw_anat_video_y98.nii.gz -x z y raw_anat_video_y98.nii.gz
gunzip raw_anat_video_y98.nii.gz
fsl5.0-fslroi raw_anat_video.nii.gz raw_anat_video_z40.nii.gz 0 -1 0 -1 31 1 0 -1
gunzip raw_anat_video_z40.nii.gz

#rawfunc
fsl5.0-fslroi raw_func_video.nii.gz raw_func_video_x48.nii.gz 48 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim raw_func_video_x48.nii.gz -z y x raw_func_video_x48.nii.gz
gunzip raw_func_video_x48.nii.gz
fsl5.0-fslroi raw_func_video.nii.gz raw_func_video_y23.nii.gz 0 -1 23 1 0 -1 0 -1
fsl5.0-fslswapdim raw_func_video_y23.nii.gz -x -z -y raw_func_video_y23.nii.gz
gunzip raw_func_video_y23.nii.gz
fsl5.0-fslroi raw_func_video.nii.gz raw_func_video_z5.nii.gz 0 -1 0 -1 5 1 0 -1
gunzip raw_func_video_z5.nii.gz

#Qwanat
fsl5.0-fslroi Qw_anat_video.nii.gz Qw_anat_video_x75.nii.gz 75 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim Qw_anat_video_x75.nii.gz -z y x Qw_anat_video_x75.nii.gz
gunzip Qw_anat_video_x75.nii.gz
fsl5.0-fslroi Qw_anat_video.nii.gz Qw_anat_video_y93.nii.gz 0 -1 93 1 0 -1 0 -1
fsl5.0-fslswapdim Qw_anat_video_y93.nii.gz -x z y Qw_anat_video_y93.nii.gz
gunzip Qw_anat_video_y93.nii.gz
fsl5.0-fslroi Qw_anat_video.nii.gz Qw_anat_video_z51.nii.gz 0 -1 0 -1 51 1 0 -1
gunzip Qw_anat_video_z51.nii.gz

#Qwfunc
fsl5.0-fslroi Qw_func_video.nii.gz Qw_func_video_x75.nii.gz 75 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim Qw_func_video_x75.nii.gz -z y x Qw_func_video_x75.nii.gz
gunzip Qw_func_video_x75.nii.gz
fsl5.0-fslroi Qw_func_video.nii.gz Qw_func_video_y93.nii.gz 0 -1 93 1 0 -1 0 -1
fsl5.0-fslswapdim Qw_func_video_y93.nii.gz -x z y Qw_func_video_y93.nii.gz
gunzip Qw_func_video_y93.nii.gz
fsl5.0-fslroi Qw_func_video.nii.gz Qw_func_video_z51.nii.gz 0 -1 0 -1 51 1 0 -1
gunzip Qw_func_video_z51.nii.gz

