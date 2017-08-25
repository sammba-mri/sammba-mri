#!/bin/bash

#bash /home/Pmamobipet/Tvx-Manips-MD_/MD_1701-Microcebe-Creation-Atlas/code4/figures4.bash

#AFAICT x128y127z140 of FSL ==x127y127z140 of AFNI (0.057,-1.438,0.057 L,A,S under RAI=DICOM)

cd /home/Pmamobipet/Tvx-Manips-MD_/MD_1701-Microcebe-Creation-Atlas/Final-Images-For-Atlas/20170516

cp ../../processed_20170411/Mc147BCBB-01-s_2011042110__ptk_server_task_2__MD1101_T2Lemur_230ZFx230ZFx230_128cp_02/Un.nii.gz Un1.nii.gz
cp ../../processed_20170411/Mc153FBA-01-s_2011082405__ptk_server_task_1__MD1101_T2Lemur_230ZFx230ZFx230_128cp_01/Un.nii.gz Un2.nii.gz
cp ../../processed_20170411/Mc169ABB-01-s_2012112902__T2Lemur_230ZFx230ZFx230_128cp_New2012_01/Un.nii.gz Un3.nii.gz

afni -noplugins -no_detach																	\
				-com 'OPEN_WINDOW axialimage opacity=3' 									\
				-com 'OPEN_WINDOW sagittalimage opacity=3' 									\
				-com 'OPEN_WINDOW coronalimage opacity=3' 									\
				-com 'SWITCH_UNDERLAY T2anat.nii.gz' 										\
				-com 'SET_DICOM_XYZ 0.120 -1.723 -0.150'									\
				-com 'SAVE_PNG axialimage Mc147BCBA_axial' 									\
				-com 'SAVE_PNG sagittalimage Mc147BCBA_sagittal' 							\
				-com 'SAVE_PNG coronalimage Mc147BCBA_coronal' 								\
				-com 'SWITCH_UNDERLAY Un.nii.gz' 											\
				-com 'SAVE_PNG axialimage Mc147BCBA_Un_axial' 								\
				-com 'SAVE_PNG sagittalimage Mc147BCBA_Un_sagittal' 						\
				-com 'SAVE_PNG coronalimage Mc147BCBA_Un_coronal' 							\
				-com 'SWITCH_UNDERLAY Un1.nii.gz' 											\
				-com 'SAVE_PNG axialimage Un1_axial' 										\
				-com 'SAVE_PNG sagittalimage Un1_sagittal' 									\
				-com 'SAVE_PNG coronalimage Un1_coronal' 									\
				-com 'SWITCH_UNDERLAY Un2.nii.gz' 											\
				-com 'SAVE_PNG axialimage Un2_axial' 										\
				-com 'SAVE_PNG sagittalimage Un2_sagittal' 									\
				-com 'SAVE_PNG coronalimage Un2_coronal' 									\
				-com 'SWITCH_UNDERLAY Un3.nii.gz' 											\
				-com 'SAVE_PNG axialimage Un3_axial' 										\
				-com 'SAVE_PNG sagittalimage Un3_sagittal' 									\
				-com 'SAVE_PNG coronalimage Un3_coronal' 									\
				-com 'SWITCH_UNDERLAY UnBm.nii.gz' 											\
				-com 'SAVE_PNG axialimage Mc147BCBA_UnBm_axial' 							\
				-com 'SAVE_PNG sagittalimage Mc147BCBA_UnBm_sagittal' 						\
				-com 'SAVE_PNG coronalimage Mc147BCBA_UnBm_coronal' 						\
				-com 'SWITCH_UNDERLAY UnBmBe.nii.gz' 										\
				-com 'SAVE_PNG axialimage Mc147BCBA_UnBmBe_axial' 							\
				-com 'SAVE_PNG sagittalimage Mc147BCBA_UnBmBe_sagittal' 					\
				-com 'SAVE_PNG coronalimage Mc147BCBA_UnBmBe_coronal' 						\
				-com 'SWITCH_UNDERLAY raw_mean.nii.gz' 										\
				-com 'SET_DICOM_XYZ 0.115 -1.437 5.405'										\
				-com 'SAVE_PNG axialimage raw_mean_axial' 									\
				-com 'SAVE_PNG sagittalimage raw_mean_sagittal' 							\
				-com 'SAVE_PNG coronalimage raw_mean_coronal' 								\
				-com 'SWITCH_UNDERLAY bons.nii.gz' 											\
				-com 'SET_DICOM_XYZ -0.057 -0.057 -0.057'									\
				-com 'SAVE_PNG axialimage bons_axial' 										\
				-com 'SAVE_PNG sagittalimage bons_sagittal' 								\
				-com 'SAVE_PNG coronalimage bons_coronal' 									\
				-com 'SWITCH_UNDERLAY UnCC_mean.nii.gz' 									\
				-com 'SET_DICOM_XYZ 0.057 -1.438 0.057'										\
				-com 'SAVE_PNG axialimage UnCC_mean_axial' 									\
				-com 'SAVE_PNG sagittalimage UnCC_mean_sagittal' 							\
				-com 'SAVE_PNG coronalimage UnCC_mean_coronal' 								\
				-com 'SWITCH_UNDERLAY UnBmBeCC_mean.nii.gz'									\
				-com 'SAVE_PNG axialimage UnBmBeCC_mean_axial' 								\
				-com 'SAVE_PNG sagittalimage UnBmBeCC_mean_sagittal' 						\
				-com 'SAVE_PNG coronalimage UnBmBeCC_mean_coronal' 							\
				-com 'SWITCH_UNDERLAY UnBmBe_mean.nii.gz'									\
				-com 'SAVE_PNG axialimage UnBmBe_mean_axial' 								\
				-com 'SAVE_PNG sagittalimage UnBmBe_mean_sagittal' 							\
				-com 'SAVE_PNG coronalimage UnBmBe_mean_coronal' 							\
				-com 'SWITCH_UNDERLAY shr1_mean.nii.gz'										\
				-com 'SAVE_PNG axialimage shr1_mean_axial' 									\
				-com 'SAVE_PNG sagittalimage shr1_mean_sagittal' 							\
				-com 'SAVE_PNG coronalimage shr1_mean_coronal' 								\
				-com 'SWITCH_UNDERLAY shr1_meanhead.nii.gz' 								\
				-com 'SAVE_PNG axialimage shr1_meanhead_axial' 								\
				-com 'SAVE_PNG sagittalimage shr1_meanhead_sagittal' 						\
				-com 'SAVE_PNG coronalimage shr1_meanhead_coronal' 							\
				-com 'SWITCH_OVERLAY shr1_count.nii.gz' 									\
				-com 'SET_PBAR_SIGN +' 														\
				-com 'SAVE_PNG axialimage shr1_meanhead_shr1_count_axial' 					\
				-com 'SAVE_PNG sagittalimage shr1_meanhead_shr1_count_sagittal' 			\
				-com 'SAVE_PNG coronalimage shr1_meanhead_shr1_count_coronal' 				\
				-com 'SEE_OVERLAY -' 														\
				-com 'SET_DICOM_XYZ -5.923 4.887 -4.083'									\
				-com 'SAVE_PNG axialimage shr1_meanhead_colliculus_axial' 					\
				-com 'SAVE_PNG sagittalimage shr1_meanhead_colliculus_sagittal' 			\
				-com 'SAVE_PNG coronalimage shr1_meanhead_colliculus_coronal' 				\
				-com 'SEE_OVERLAY +' 														\
				-com 'SAVE_PNG axialimage shr1_meanhead_colliculus_shr1_count_axial' 		\
				-com 'SAVE_PNG sagittalimage shr1_meanhead_colliculus_shr1_count_sagittal' 	\
				-com 'SAVE_PNG coronalimage shr1_meanhead_colliculus_shr1_count_coronal' 	\
				-com 'SEE_OVERLAY -' 														\
				-com 'SET_DICOM_XYZ 0.057 -1.438 0.057'										\
				-com 'SWITCH_UNDERLAY shr2_mean.nii.gz' 									\
				-com 'SAVE_PNG axialimage shr2_mean_axial' 									\
				-com 'SAVE_PNG sagittalimage shr2_mean_sagittal' 							\
				-com 'SAVE_PNG coronalimage shr2_mean_coronal' 								\
				-com 'SWITCH_UNDERLAY shr2_meanhead.nii.gz' 								\
				-com 'SAVE_PNG axialimage shr2_meanhead_axial' 								\
				-com 'SAVE_PNG sagittalimage shr2_meanhead_sagittal' 						\
				-com 'SAVE_PNG coronalimage shr2_meanhead_coronal' 							\
				-com 'SEE_OVERLAY +'														\
				-com 'SWITCH_OVERLAY shr2_count.nii.gz'										\
				-com 'SAVE_PNG axialimage shr2_meanhead_shr2_count_axial' 					\
				-com 'SAVE_PNG sagittalimage shr2_meanhead_shr2_count_sagittal' 			\
				-com 'SAVE_PNG coronalimage shr2_meanhead_shr2_count_coronal' 				\
				-com 'SEE_OVERLAY -'														\
				-com 'SWITCH_UNDERLAY aff3_meanhead.nii.gz' 								\
				-com 'SAVE_PNG axialimage aff3_meanhead_axial' 								\
				-com 'SAVE_PNG sagittalimage aff3_meanhead_sagittal' 						\
				-com 'SAVE_PNG coronalimage aff3_meanhead_coronal' 							\
				-com 'SEE_OVERLAY +'														\
				-com 'SWITCH_OVERLAY aff3_unionmask.nii.gz'									\
				-com 'PBAR_ROTATE -'														\
				-com 'SAVE_PNG axialimage aff3_meanhead_aff3_unionmask_axial' 				\
				-com 'SAVE_PNG sagittalimage aff3_meanhead_aff3_unionmask_sagittal' 		\
				-com 'SAVE_PNG coronalimage aff3_meanhead_aff3_unionmask_coronal' 			\
				-com 'SWITCH_OVERLAY aff3_unionmaskdil7.nii.gz'								\
				-com 'SAVE_PNG axialimage aff3_meanhead_aff3_unionmaskdil7_axial' 			\
				-com 'SAVE_PNG sagittalimage aff3_meanhead_aff3_unionmaskdil7_sagittal' 	\
				-com 'SAVE_PNG coronalimage aff3_meanhead_aff3_unionmaskdil7_coronal' 		\
				-com 'SEE_OVERLAY -'														\
				-com 'SWITCH_UNDERLAY Qw1_meanhead.nii.gz' 									\
				-com 'SAVE_PNG axialimage Qw1_meanhead_axial' 								\
				-com 'SAVE_PNG sagittalimage Qw1_meanhead_sagittal' 						\
				-com 'SAVE_PNG coronalimage Qw1_meanhead_coronal' 							\
				-com 'SEE_OVERLAY +'														\
				-com 'SAVE_PNG axialimage Qw1_meanhead_aff3_unionmaskdil7_axial' 			\
				-com 'SAVE_PNG sagittalimage Qw1_meanhead_aff3_unionmaskdil7_sagittal' 		\
				-com 'SAVE_PNG coronalimage Qw1_meanhead_aff3_unionmaskdil7_coronal' 		\
				-com 'SEE_OVERLAY -'														\
				-com 'SWITCH_UNDERLAY Qw2_meanhead.nii.gz' 									\
				-com 'SAVE_PNG axialimage Qw2_meanhead_axial' 								\
				-com 'SAVE_PNG sagittalimage Qw2_meanhead_sagittal' 						\
				-com 'SAVE_PNG coronalimage Qw2_meanhead_coronal' 							\
				-com 'SEE_OVERLAY +'														\
				-com 'SAVE_PNG axialimage Qw2_meanhead_aff3_unionmaskdil7_axial' 			\
				-com 'SAVE_PNG sagittalimage Qw2_meanhead_aff3_unionmaskdil7_sagittal' 		\
				-com 'SAVE_PNG coronalimage Qw2_meanhead_aff3_unionmaskdil7_coronal' 		\
				-com 'SEE_OVERLAY -'														\
				-com 'SWITCH_UNDERLAY Qw3_meanhead.nii.gz' 									\
				-com 'SAVE_PNG axialimage Qw3_meanhead_axial' 								\
				-com 'SAVE_PNG sagittalimage Qw3_meanhead_sagittal' 						\
				-com 'SAVE_PNG coronalimage Qw3_meanhead_coronal' 							\
				-com 'SWITCH_UNDERLAY Qw4_meanhead.nii.gz' 									\
				-com 'SAVE_PNG axialimage Qw4_meanhead_axial' 								\
				-com 'SAVE_PNG sagittalimage Qw4_meanhead_sagittal' 						\
				-com 'SAVE_PNG coronalimage Qw4_meanhead_coronal' 							\
				-com 'SWITCH_UNDERLAY Na_meanhead.nii.gz' 									\
				-com 'SAVE_PNG axialimage Na_meanhead_axial' 								\
				-com 'SAVE_PNG sagittalimage Na_meanhead_sagittal' 							\
				-com 'SAVE_PNG coronalimage Na_meanhead_coronal' 							\
				-com 'SET_XHAIRS OFF'														\
				-com 'SWITCH_UNDERLAY raw_mean.nii.gz' 										\
				-com 'SET_DICOM_XYZ 0.115 -1.437 5.405'										\
				-com 'SAVE_PNG axialimage raw_mean_noxh_axial' 								\
				-com 'SAVE_PNG sagittalimage raw_mean_noxh_sagittal' 						\
				-com 'SAVE_PNG coronalimage raw_mean_noxh_coronal' 							\
				-com 'SWITCH_UNDERLAY UnCC_mean.nii.gz' 									\
				-com 'SET_DICOM_XYZ 0.057 -1.438 0.057'										\
				-com 'SAVE_PNG axialimage UnCC_mean_noxh_axial' 							\
				-com 'SAVE_PNG sagittalimage UnCC_mean_noxh_sagittal' 						\
				-com 'SAVE_PNG coronalimage UnCC_mean_noxh_coronal' 						\
				-com 'SWITCH_UNDERLAY shr1_meanhead.nii.gz' 								\
				-com 'SAVE_PNG axialimage shr1_meanhead_noxh_axial' 						\
				-com 'SAVE_PNG sagittalimage shr1_meanhead_noxh_sagittal' 					\
				-com 'SAVE_PNG coronalimage shr1_meanhead_noxh_coronal' 					\
				-com 'SWITCH_UNDERLAY shr2_meanhead.nii.gz' 								\
				-com 'SAVE_PNG axialimage shr2_meanhead_noxh_axial' 						\
				-com 'SAVE_PNG sagittalimage shr2_meanhead_noxh_sagittal' 					\
				-com 'SAVE_PNG coronalimage shr2_meanhead_noxh_coronal' 					\
				-com 'SWITCH_UNDERLAY aff3_meanhead.nii.gz' 								\
				-com 'SAVE_PNG axialimage aff3_meanhead_noxh_axial' 						\
				-com 'SAVE_PNG sagittalimage aff3_meanhead_noxh_sagittal' 					\
				-com 'SAVE_PNG coronalimage aff3_meanhead_noxh_coronal' 					\
				-com 'SWITCH_UNDERLAY Qw1_meanhead.nii.gz' 									\
				-com 'SAVE_PNG axialimage Qw1_meanhead_noxh_axial' 							\
				-com 'SAVE_PNG sagittalimage Qw1_meanhead_noxh_sagittal' 					\
				-com 'SAVE_PNG coronalimage Qw1_meanhead_noxh_coronal' 						\
				-com 'SWITCH_UNDERLAY Qw2_meanhead.nii.gz' 									\
				-com 'SAVE_PNG axialimage Qw2_meanhead_noxh_axial' 							\
				-com 'SAVE_PNG sagittalimage Qw2_meanhead_noxh_sagittal' 					\
				-com 'SAVE_PNG coronalimage Qw2_meanhead_noxh_coronal' 						\
				-com 'SWITCH_UNDERLAY Qw3_meanhead.nii.gz' 									\
				-com 'SAVE_PNG axialimage Qw3_meanhead_noxh_axial' 							\
				-com 'SAVE_PNG sagittalimage Qw3_meanhead_noxh_sagittal' 					\
				-com 'SAVE_PNG coronalimage Qw3_meanhead_noxh_coronal' 						\
				-com 'SWITCH_UNDERLAY Qw4_meanhead.nii.gz' 									\
				-com 'SAVE_PNG axialimage Qw4_meanhead_noxh_axial' 							\
				-com 'SAVE_PNG sagittalimage Qw4_meanhead_noxh_sagittal' 					\
				-com 'SAVE_PNG coronalimage Qw4_meanhead_noxh_coronal' 						\
				-com 'SWITCH_UNDERLAY Na_meanhead.nii.gz' 									\
				-com 'SAVE_PNG axialimage Na_meanhead_noxh_axial' 							\
				-com 'SAVE_PNG sagittalimage Na_meanhead_noxh_sagittal' 					\
				-com 'SAVE_PNG coronalimage Na_meanhead_noxh_coronal' 						\
				-com 'QUITT' 																\
				../../processed_20170411/*.nii.gz											\
				../../processed_20170411/Mc147BCBA-01-s_2012101102__T2Lemur_230ZFx230ZFx230_128cp_New2012_01/*.nii.gz \
				*.nii.gz

convert -size 256x256 canvas:black filler.png #used below, couldn't make null: work in montage
convert -size 256x256 canvas:"rgb(0,0,1)" fillerblackblue.png #otherwise noxh montages end up with excessive contrast
				
for stem in 							\
	Mc147BCBA							\
	Un1									\
	Un2									\
	Un3									\
	Mc147BCBA_Un						\
	Mc147BCBA_UnBm						\
	Mc147BCBA_UnBmBe					\
	raw_mean							\
	UnCC_mean							\
	UnBmBeCC_mean						\
	UnBmBe_mean							\
	bons								\
	shr1_mean							\
	shr1_meanhead						\
	shr1_meanhead_shr1_count			\
	shr1_meanhead_colliculus			\
	shr1_meanhead_colliculus_shr1_count	\
	shr2_mean							\
	shr2_meanhead						\
	shr2_meanhead_shr2_count			\
	aff3_meanhead						\
	aff3_meanhead_aff3_unionmask		\
	aff3_meanhead_aff3_unionmaskdil7	\
	Qw1_meanhead						\
	Qw1_meanhead_aff3_unionmaskdil7		\
	Qw2_meanhead						\
	Qw2_meanhead_aff3_unionmaskdil7		\
	Qw3_meanhead						\
	Qw4_meanhead 						\
	Na_meanhead							\
	; do
	montage "$stem"_sagittal.png "$stem"_coronal.png filler.png "$stem"_axial.png -geometry +0+0 -tile 2x2 -background 'black' "$stem"_montage.png
done

for stem in 							\
	raw_mean							\
	UnCC_mean							\
	shr1_meanhead						\
	shr2_meanhead						\
	aff3_meanhead						\
	Qw1_meanhead						\
	Qw2_meanhead						\
	Qw3_meanhead						\
	Qw4_meanhead 						\
	Na_meanhead							\
	; do
	montage "$stem"_noxh_sagittal.png "$stem"_noxh_coronal.png fillerblackblue.png "$stem"_noxh_axial.png -geometry +0+0 -tile 2x2 -background 'black' "$stem"_noxh_montage.png
done
