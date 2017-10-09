#!/bin/bash

#bash /home/nadkarni/git/sammba-mri/sammba/projects/MIRCen_sammba-mri_mouselemur/videos_MIRCen_sammba-mri_mouselemur.bash

export AFNI_DECONFLICT=OVERWRITE #saves a lot of grief

cd /home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouselemur/MIRCen_sammba-mri_mouselemur/analysis20170928

3dTcat -prefix raw_anat_video.nii.gz $(find -name 'anat_n0.nii.gz')
3dTcat -prefix Qw_anat_video.nii.gz $(find -name 'anat_n0_UnAaQw.nii.gz')
3dTcat -prefix raw_func_video.nii.gz $(cat ../rawlist_96x96_06_short.txt)
3dTcat -prefix Qw_func_video.nii.gz $(cat ../Qwlist_96x96_06_short.txt)

#rawanat
fsl5.0-fslroi raw_anat_video.nii.gz raw_anat_video_x79.nii.gz 79 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim raw_anat_video_x79.nii.gz -z y x raw_anat_video_x79.nii.gz
gunzip raw_anat_video_x79.nii.gz
fsl5.0-fslroi raw_anat_video.nii.gz raw_anat_video_y80.nii.gz 0 -1 80 1 0 -1 0 -1
fsl5.0-fslswapdim raw_anat_video_y80.nii.gz -x -z -y raw_anat_video_y80.nii.gz
gunzip raw_anat_video_y80.nii.gz
fsl5.0-fslroi raw_anat_video.nii.gz raw_anat_video_z37.nii.gz 0 -1 0 -1 37 1 0 -1
gunzip raw_anat_video_z37.nii.gz

#rawfunc
fsl5.0-fslroi raw_func_video.nii.gz raw_func_video_x45.nii.gz 45 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim raw_func_video_x45.nii.gz -z y x raw_func_video_x45.nii.gz
gunzip raw_func_video_x45.nii.gz
fsl5.0-fslroi raw_func_video.nii.gz raw_func_video_y38.nii.gz 0 -1 38 1 0 -1 0 -1
fsl5.0-fslswapdim raw_func_video_y38.nii.gz -x z y raw_func_video_y38.nii.gz
gunzip raw_func_video_y38.nii.gz
fsl5.0-fslroi raw_func_video.nii.gz raw_func_video_z14.nii.gz 0 -1 0 -1 14 1 0 -1
gunzip raw_func_video_z14.nii.gz

#Qwanat
fsl5.0-fslroi Qw_anat_video.nii.gz Qw_anat_video_x79.nii.gz 79 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim Qw_anat_video_x79.nii.gz -z y x Qw_anat_video_x79.nii.gz
gunzip Qw_anat_video_x79.nii.gz
fsl5.0-fslroi Qw_anat_video.nii.gz Qw_anat_video_y76.nii.gz 0 -1 76 1 0 -1 0 -1
fsl5.0-fslswapdim Qw_anat_video_y76.nii.gz -x z y Qw_anat_video_y76.nii.gz
gunzip Qw_anat_video_y76.nii.gz
fsl5.0-fslroi Qw_anat_video.nii.gz Qw_anat_video_z47.nii.gz 0 -1 0 -1 47 1 0 -1
gunzip Qw_anat_video_z47.nii.gz

#Qwfunc
fsl5.0-fslroi Qw_func_video.nii.gz Qw_func_video_x79.nii.gz 79 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim Qw_func_video_x79.nii.gz -z y x Qw_func_video_x79.nii.gz
gunzip Qw_func_video_x79.nii.gz
fsl5.0-fslroi Qw_func_video.nii.gz Qw_func_video_y76.nii.gz 0 -1 76 1 0 -1 0 -1
fsl5.0-fslswapdim Qw_func_video_y76.nii.gz -x z y Qw_func_video_y76.nii.gz
gunzip Qw_func_video_y76.nii.gz
fsl5.0-fslroi Qw_func_video.nii.gz Qw_func_video_z47.nii.gz 0 -1 0 -1 47 1 0 -1
gunzip Qw_func_video_z47.nii.gz
