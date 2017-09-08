#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE

procdir=$(readlink -e $1)
T1blood=$2
lambda=$3
multiplier=$4
T1guess=$5
mccores=$6

echo "1 procdir" $procdir
echo "2 T1blood" $T1blood
echo "3 lambda" $lambda
echo "4 multiplier" $multiplier
echo "5 T1guess" $T1guess
echo "6 mccores" $mccores

echo "processing perfFAIREPI for bias, M0, T1, rCBF and CBF..."

unset a_perfs
perfs=$(find $procdir -name perfFAIREPI_n[0-9].nii.gz)
readarray -t a_perfs <<< "$perfs"

if [ -r "${a_perfs[0]}" ]; then
	counter=0
	for perffile in ${a_perfs[@]}; do
		p=$procdir/$(basename -s .nii.gz $perffile)
		counter=$(expr $counter + 1); echo "processing $counter of ${#a_perfs[@]}"
		echo "perffile $perffile"
		perfFAIREPItoNIfTI.R $perffile $T1blood $lambda $multiplier $T1guess $mccores "$p"_calc
		#R script does not produce a good header
		3dcalc -a "$p".nii.gz[0-13] -expr 'a' -prefix "$p"_0-13.nii.gz
		fsl5.0-fslcpgeom "$p"_0-13.nii.gz "$p"_calc.nii.gz
		rm -f "$p"_0-13.nii.gz
	done
else
	echo "no perfFAIREPIs"
fi

unset a_perfs