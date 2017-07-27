#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE

procdir=$(readlink -e $1)
subpipeline=$2
template=$(readlink -e $3)
tmpdir=$(readlink -e $4)
registerfunctional=$5
rminterfiles=$6
brainvol=$7

echo "1 procdir" $procdir
echo "2 subpipeline" $subpipeline
echo "3 template" $template
echo "4 tmpdir" $tmpdir
echo "5 registerfunctional" $registerfunctional
echo "6 rminterfiles" $rminterfiles
echo "7 brainvol" $brainvol

tmpselprefix="$RANDOM"
echo 57 77 146 > $tmpdir/"$tmpselprefix"_sel.1D
selfile=$tmpdir/"$tmpselprefix"_sel.1D

unset a_rss a_Qwanats

rss=$(find $procdir -name rs*_n[0-9].nii.gz)
readarray -t a_rss <<< "$rss"

Qwanats=$(find $procdir -name anat_n[0-9]_*AaQw.nii.gz)
readarray -t a_Qwanats <<< "$Qwanats"

rminterregfiles="yes"
if [[ $rminterfiles == "no" ]] ; then rminterregfiles="no" ; fi

echo "processing resting state data..."
if [[ -r "${a_rss[0]}" && -r "${a_Qwanats[0]}" ]]; then
	counter=0
	for file in ${a_rss[@]}; do
		
		p=$(dirname $file)/$(basename -s .nii.gz $file)

		if [[ $(3dinfo -nv "$p".nii.gz) -ge 15 ]]; then

			counter=$(expr $counter + 1); echo "processing $counter of ${#a_rss[@]}"

			tpattern.R $(3dinfo -tr -nk "$p".nii.gz) "$p"_tpattern.txt
			3dTshift -prefix "$p"_Ts.nii.gz -tpattern @"$p"_tpattern.txt "$p".nii.gz
			#fsl5.0-fslcpgeom "$p".nii.gz "$p"_Ts.nii.gz -d #remarkably it's not needed	
			thresh=$(3dClipLevel "$p"_Ts.nii.gz)
			fsl5.0-fslmaths "$p"_Ts.nii.gz -thr $thresh "$p"_TsTm.nii.gz
			3dvolreg -prefix "$p"_TsTmVr.nii.gz -dfile "$p"_TsTmVr.dfile.1D -1Dfile "$p"_TsTmVr.1Dfile.1D -1Dmatrix_save "$p"_TsTmVr.aff12.1D "$p"_TsTm.nii.gz
			3dAllineate -input "$p"_Ts.nii.gz -master "$p"_Ts.nii.gz -prefix "$p"_TsAv.nii.gz -1Dmatrix_apply "$p"_TsTmVr.aff12.1D
			fsl5.0-fslcpgeom "$p"_TsTmVr.nii.gz "$p"_TsAv.nii.gz #3dAllineate removes the obliquity. this is not a good way to readd it as removes motion correction info in the header if it were an AFNI file...as it happens it's NIfTI which does not store that so irrelevant!
			3dTstat -mean -prefix "$p"_TsAvAv.nii.gz "$p"_TsAv.nii.gz #create a (hopefully) nice mean image for use in the registration
			N3BiasFieldCorrection 3 "$p"_TsAvAv.nii.gz "$p"_TsAvAvN3.nii.gz #N4 fails for some reason. Not tried 3dUnifize yet.

			br=$(basename -s .nii.gz $file)
			imfilearray=("$br"_TsAvAvN3 "$br"_TsAv)
			echo ${imfilearray[@]} > "$p"_imfilearray.txt
			bash $subpipeline $procdir $(basename -s AaQw.nii.gz ${a_Qwanats[0]}) "$p"_imfilearray.txt atlas_Na1 no $template $rminterregfiles $tmpdir 0.1 $registerfunctional $brainvol

			1d_tool.py -infile "$p"_TsTmVr.1Dfile.1D -derivative -censor_prev_TR -collapse_cols euclidean_norm -moderate_mask -1.2 1.2 -show_censor_count -write_censor "$p"_censor.1D -write_CENSORTR "$p"_CENSORTR.1D -overwrite
			3dToutcount -mask $procdir/atlas_Na1_Op_"$br"_TsAvAvN3.nii.gz -fraction "$p"_TsAv_NaMe.nii.gz > "$p"_outcount.1D
			1deval -a "$p"_censor.1D -b "$p"_outcount.1D -expr 'a*ispositive(0.2-b)' > "$p"_outcount_Ev.1D

			3dROIstats -mask $procdir/atlas_Na1_Op_"$br"_TsAvAvN3.nii.gz -roisel $selfile -quiet "$p"_TsAv_NaMe.nii.gz > "$p"_ROIs.1D
			Rscript -e 'write(rowMeans(read.table("stdin")),commandArgs(TRUE)[1],ncolumns=1)' "$p"_ventricles.1D < "$p"_ROIs.1D

			1d_tool.py -infile "$p"_TsTmVr.1Dfile.1D -set_nruns 1 -demean -overwrite -write "$p"_motion_demean.1D
			1d_tool.py -infile "$p"_TsTmVr.1Dfile.1D -set_nruns 1 -derivative -demean -overwrite -write "$p"_motion_deriv.1D
			1dBport -nodata $(3dinfo -nv -tr "$p"_TsAv_NaMe.nii.gz) -band 0.01 0.1 -invert -nozero > "$p"_bandpass_001-01.1D
			1dBport -nodata $(3dinfo -nv -tr "$p"_TsAv_NaMe.nii.gz) -band 0.01 0.3 -invert -nozero > "$p"_bandpass_001-03.1D

			#3dNwarpApply -nwarp "$T1anat"_"$Bc"AaQw_WARP.nii.gz "$T1anat"_"$Bc"BmBeAl.aff12.1D "$T1anat"_"$Bc"_Op_"$br"_TsAvAvN3.aff12.1D -source "$p"_TsAv_NaMe.nii.gz -master $template -prefix "$p"_TsAv_NaMeNa.nii.gz
			
			if [[ $rminterfiles == "yes" ]] ; then rm -f "$p"_Ts.nii.gz "$p"_TsTm.nii.gz "$p"_TsTmVr.nii.gz "$p"_TsAv.nii.gz ; fi
			
		fi
	done
		else
	echo "missing files"
fi	
