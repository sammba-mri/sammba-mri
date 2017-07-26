#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE

procdir=$(readlink -e $1)

echo "1 procdir" $procdir

echo "processing perfFAIREPI for bias, M0, T1, rCBF and CBF..."

unset a_perfs
perfs=$(find $procdir -name perfFAIREPI_n[0-9].nii.gz)
readarray -t a_perfs <<< "$perfs"

if [ -r "${a_perfs[0]}" ]; then
	counter=0
	for file in ${a_perfs[@]}; do
		p=$(dirname $file)/$(basename -s .nii.gz $file)
		pnum=$(basename $p | sed s/perfFAIREPI_n//g)
		TIsfile=$(dirname $file)/perfFAIREPI_TIs_n$pnum.txt
		TIs=$(cat $TIsfile)
		counter=$(expr $counter + 1); echo "processing $counter of ${#a_perfs[@]}"
		echo "file $file"
		echo "TIsfile $TIsfile"
		echo "TIs $TIs"
		perfFAIREPItoNIfTI.r "$p".nii.gz 2800 0.9 $TIs 6000000 1600 16 "$p"_calc $(which perfFAIREPI+T1_T2map_RARE_fitters.r)
		#R script does not produce a good header
		3dcalc -a "$p".nii.gz[0-13] -expr 'a' -prefix "$p"_0-13.nii.gz
		fsl5.0-fslcpgeom "$p"_0-13.nii.gz "$p"_calc.nii.gz
		rm -f "$p"_0-13.nii.gz
	done
else
	echo "no perfFAIREPIs"
fi

unset a_perfs