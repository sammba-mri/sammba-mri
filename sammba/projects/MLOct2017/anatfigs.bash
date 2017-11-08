#!/bin/bash

cd /home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouselemur/MLOct2017/analysis20171030

afni -noplugins -no_detach														 \
				-com 'SET_XHAIRS OFF'											 \
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
