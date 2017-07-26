#!/bin/bash

pipeline=$1
MRIsessionslist=$(readlink -e $2)
analysisdir=$(readlink -e $3)
brain=$(readlink -e $4)
atlas=$(readlink -e $5)
mask=$(readlink -e $6)
head=$(readlink -e $7)
headweight=$(readlink -e $8)
basetype=$9
dofolderoverwrite=${10}
tmpdir=${11}
registerfunctional=${12}
subpipeline=${13}
Urad=${14}
brainvol=${15}
scale=${16}

echo "1 pipeline" $pipeline
echo "2 MRIsessionslist" $MRIsessionslist
echo "3 analysisdir" $analysisdir
echo "4 brain" $brain
echo "5 atlas" $atlas
echo "6 mask" $mask
echo "7 head" $head
echo "8 headweight" $headweight
echo "9 basetype" $basetype
echo "10 dofolderoverwrite" $dofolderoverwrite
echo "11 tmpdir" $tmpdir
echo "12 registerfunctional" $registerfunctional
echo "13 subpipeline" $subpipeline
echo "14 Urad" $Urad
echo "15 brainvol" $brainvol
echo "16 scale" $scale

tmpdirsuffix="$RANDOM"
tmpdir=$tmpdir/tmpanaldir_$tmpdirsuffix
mkdir $tmpdir

nlines=$(wc -l $MRIsessionslist | awk '{print $1}')
for ((a=2; a<=$nlines; a++)); do

	rawimdir=$(head -n$a < $MRIsessionslist | tail -n1 | awk '{print $3}')
	NIfTIdir=$analysisdir/$(basename $rawimdir)
	
	runpipeline=no

	#amongst other things, I think the following is fail-safe against an empty rawimdir string
	if [[ $dofolderoverwrite == "yes" ]]; then
		if [ -r "$NIfTIdir" ]; then
			echo "session NIfTI directory exists, creating session NIfTI directory with yymmdd_HHMMSS suffix..."
			NIfTIdir="$NIfTIdir"_"$(date +%Y%m%d_%H%M%S)"
			runpipeline=yes
		else
			echo "session NIfTI directory does not exist, creating session NIfTI directory..."
			runpipeline=yes
		fi
		mkdir $NIfTIdir
	else
		if [ -r "$NIfTIdir" ]; then
			echo "session NIfTI directory exists, not creating or overwriting session NIfTI directory"
			runpipeline=no
		else
			echo "session NIfTI directory does not exist, creating session NIfTI directory..."
			runpipeline=yes
			mkdir $NIfTIdir
		fi
	fi

	if [[ $runpipeline == "yes" ]]; then
#skips 10 ($dofolderoverwrite)
#				  1			2		  3		 4		5	  6		7			8		  9		  10				  11		   12	 13		   14
		$pipeline $rawimdir $NIfTIdir $brain $atlas $mask $head $headweight $basetype $tmpdir $registerfunctional $subpipeline $Urad $brainvol $scale
	fi
	
done

rm -rfd $tmpdir
