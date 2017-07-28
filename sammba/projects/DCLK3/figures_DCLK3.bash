#!/bin/bash

#bash /home/nadkarni/git/sammba-mri/sammba/projects/DCLK3/figures_DCLK3.bash

#AFAICT x79y87z110 of FSL == x79y87z110 of AFNI (0.050,-1.050,-0.750 L,A,I under RAI=DICOM)

cd /home/Pplateforme/Plate-forme_RMN/NachiketNadkarni/DCLK3/processed_20170727

afni -noplugins -no_detach																	\
				-com 'OPEN_WINDOW axialimage'												\
				-com 'OPEN_WINDOW sagittalimage'											\
				-com 'OPEN_WINDOW coronalimage'												\
				-com 'SET_DICOM_XYZ 0.050 -0.950 -0.450'									\
				-com 'SWITCH_UNDERLAY anatCC.nii.gz' 										\
				-com 'SAVE_PNG axialimage anatCC_axial' 									\
				-com 'SAVE_PNG sagittalimage anatCC_sagittal' 							\
				-com 'SAVE_PNG coronalimage anatCC_coronal' 								\
				-com 'SWITCH_OVERLAY atlas_Na.nii.gz'										\
				-com 'SET_PBAR_ALL +99 1.0 ROI_i256'										\
				-com 'SAVE_PNG axialimage anatCCatlas_axial' 								\
				-com 'SAVE_PNG sagittalimage anatCCatlas_sagittal' 						\
				-com 'SAVE_PNG coronalimage anatCCatlas_coronal' 							\
				-com 'SET_DICOM_XYZ -1.476 -0.438 3.690'									\
				-com 'SWITCH_UNDERLAY raw_mean.nii.gz' 										\
				-com 'SET_FUNC_VISIBLE -'													\
				-com 'SAVE_PNG axialimage raw_mean_axial' 									\
				-com 'SAVE_PNG sagittalimage raw_mean_sagittal' 							\
				-com 'SAVE_PNG coronalimage raw_mean_coronal' 								\
				-com 'SET_DICOM_XYZ 0.050 -1.050 -0.750'									\
				-com 'SWITCH_UNDERLAY Qw4_meanhead.nii.gz' 									\
				-com 'SAVE_PNG axialimage Qw4_meanhead_axial' 								\
				-com 'SAVE_PNG sagittalimage Qw4_meanhead_sagittal' 						\
				-com 'SAVE_PNG coronalimage Qw4_meanhead_coronal' 							\
				-com 'SWITCH_OVERLAY dorrbrainAlQw.nii.gz' 									\
				-com 'SET_FUNC_VISIBLE +'													\
				-com 'SET_PBAR_ALL +99 1.0 gray_scale' 										\
				-com 'SAVE_PNG axialimage dorrbrainAlQw_axial' 								\
				-com 'SAVE_PNG sagittalimage dorrbrainAlQw_sagittal' 						\
				-com 'SAVE_PNG coronalimage dorrbrainAlQw_coronal' 							\
				-com 'SWITCH_OVERLAY dorratlasAlQw.nii.gz'									\
				-com 'SET_PBAR_ALL +99 1.0 ROI_i256'										\
				-com 'SAVE_PNG axialimage dorratlasAlQw_axial' 								\
				-com 'SAVE_PNG sagittalimage dorratlasAlQw_sagittal' 						\
				-com 'SAVE_PNG coronalimage dorratlasAlQw_coronal' 							\
				-com 'SET_FUNC_VISIBLE -'													\
				-com 'SWITCH_UNDERLAY UnCCQw4.nii.gz' 										\
				-com 'SAVE_PNG axialimage UnCCQw4_axial' 									\
				-com 'SAVE_PNG sagittalimage UnCCQw4_sagittal' 								\
				-com 'SAVE_PNG coronalimage UnCCQw4_coronal' 								\
				-com 'SWITCH_UNDERLAY head100.nii.gz' 										\
				-com 'SAVE_PNG axialimage head100_axial' 									\
				-com 'SAVE_PNG sagittalimage head100_sagittal' 								\
				-com 'SAVE_PNG coronalimage head100_coronal' 								\
				-com 'SWITCH_OVERLAY labels100.nii.gz'										\
				-com 'SET_FUNC_VISIBLE +'													\
				-com 'SET_PBAR_ALL +99 1.0 ROI_i256'										\
				-com 'SAVE_PNG axialimage head100labels100_axial' 							\
				-com 'SAVE_PNG sagittalimage head100labels100_sagittal' 					\
				-com 'SAVE_PNG coronalimage head100labels100_coronal' 						\
				-com 'SWITCH_UNDERLAY brain100.nii.gz' 										\
				-com 'SAVE_PNG axialimage brain100labels100_axial' 							\
				-com 'SAVE_PNG sagittalimage brain100labels100_sagittal' 					\
				-com 'SAVE_PNG coronalimage brain100labels100_coronal' 						\
				-com 'SET_FUNC_VISIBLE -'													\
				-com 'SAVE_PNG axialimage brain100_axial' 							\
				-com 'SAVE_PNG sagittalimage brain100_sagittal' 					\
				-com 'SAVE_PNG coronalimage brain100_coronal' 						\
				-com 'QUITT' 																\
				MRIsessions/*.nii.gz																	\
				MRIsessions/030217_DCLK3_0110_WT/*.nii.gz												\
				/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MRIatlases/MIC40C57/20170223/*.nii.gz

convert -size 200x160 canvas:black filler.png #used below, couldn't make null: work in montage
convert -size 164x126 canvas:black filler2.png #used below, couldn't make null: work in montage
				
for stem in 							\
	raw_mean							\
	Qw4_meanhead						\
	dorrbrainAlQw						\
	dorratlasAlQw						\
	anatCC							\
	anatCCatlas						\
	UnCCQw4								\
	; do
	montage "$stem"_sagittal.png "$stem"_coronal.png filler.png "$stem"_axial.png -geometry +0+0 -tile 2x2 -background 'black' "$stem"_montage.png
done

for stem in 							\
	head100								\
	head100labels100					\
	brain100							\
	brain100labels100					\
; do
	montage "$stem"_sagittal.png "$stem"_coronal.png filler2.png "$stem"_axial.png -geometry +0+0 -tile 2x2 -background 'black' "$stem"_montage.png
done
