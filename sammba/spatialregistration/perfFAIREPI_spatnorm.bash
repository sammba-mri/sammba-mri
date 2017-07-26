#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE

procdir=$(readlink -e $1)
subpipeline=$2
template=$(readlink -e $3)
tmpdir=$(readlink -e $4)
registerfunctional=$5
scale=$6
brainvol=$7

echo "1 procdir" $procdir
echo "2 subpipeline" $subpipeline
echo "3 template" $template
echo "4 tmpdir" $tmpdir
echo "5 registerfunctional" $registerfunctional
echo "6 scale" $scale

echo "processing perfFAIREPI(extraTIs) for spatial normalization..."

unset a_perfs a_Qwanats

perfs=$(find $procdir -name perfFAIREPI*_n[0-9]_calc.nii.gz)
readarray -t a_perfs <<< "$perfs"

Qwanats=$(find $procdir -name anat_n[0-9]_*AaQw.nii.gz)
readarray -t a_Qwanats <<< "$Qwanats"

if [[ -r "${a_perfs[0]}" && -r "${a_Qwanats[0]}" ]]; then
	for file in ${a_perfs[@]}; do
		p=$(dirname $file)/$(basename -s _calc.nii.gz $file)
		3dTstat -prefix "$p"_M0.nii.gz "$p".nii.gz[2,8] #for use with later registration
		N3BiasFieldCorrection 3 "$p"_M0.nii.gz "$p"_M0_N3.nii.gz #N4 fails for some reason. Not tried 3dUnifize yet.
		bp=$(basename -s _calc.nii.gz $file)
		imfilearray=("$bp"_M0_N3 "$bp" "$bp"_calc)
		echo ${imfilearray[@]} > "$p"_imfilearray.txt
		bash $subpipeline $procdir $(basename -s AaQw.nii.gz ${a_Qwanats[0]}) "$p"_imfilearray.txt atlas_Na1 yes $template yes $tmpdir $scale $registerfunctional $brainvol
	done
else
	echo "missing files"
fi
	