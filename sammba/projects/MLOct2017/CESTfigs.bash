#!/bin/bash

cd /home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouselemur/MLOct2017/CEST

afni -noplugins -no_detach														 \
				-com 'SWITCH_OVERLAY 2.nii.gz' 									 \
				-com 'SET_THRESHNEW 0' 											 \
				-com 'SET_FUNC_RANGE 0.2' 										 \
				-com 'SET_PBAR_SIGN +'											 \
				-com 'SAVE_PNG axialimage 2_axial' 								 \
				-com 'SAVE_PNG sagittalimage 2_sagittal' 						 \
				-com 'SAVE_PNG coronalimage 2_coronal' 							 \
				-com 'SWITCH_OVERLAY 4.nii.gz'									 \
				-com 'SAVE_PNG axialimage 4_axial' 								 \
				-com 'SAVE_PNG sagittalimage 4_sagittal' 						 \
				-com 'SAVE_PNG coronalimage 4_coronal' 							 \
				-com 'QUITT'													 \
				../../RsLemurs/analysis20171017/Qw4_meanhead_Na_200.nii.gz 		 \
				0.nii.gz 														 \
				1.nii.gz 														 \
				2.nii.gz 														 \
				3.nii.gz 														 \
				4.nii.gz 														 \
				5.nii.gz

convert -size 172x172 canvas:black filler.png

for stem in 																	 \
	2																			 \
	4																			 \
	; do
	montage "$stem"_sagittal.png "$stem"_coronal.png filler.png "$stem"_axial.png -geometry +0+0 -tile 2x2 -background 'black' "$stem"_montage.png
done
