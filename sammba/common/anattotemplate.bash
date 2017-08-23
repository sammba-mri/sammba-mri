#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE

procdir=$(readlink -e $1)
autofind=$2
biascorrector=$3
brainvol=$4
scale=$5
brain=$(readlink -e $6)
atlas=$(readlink -e $7)
mask=$(readlink -e $8)
head=$(readlink -e $9)
headweight=$(readlink -e ${10})
basetype=${11}
Urad=${12}

echo "1 procdir" $procdir
echo "2 autofind" $autofind
echo "3 biascorrector" $biascorrector
echo "4 brainvol" $brainvol
echo "5 scale" $scale
echo "6 brain" $brain
echo "7 atlas" $atlas
echo "8 mask" $mask
echo "9 head" $head
echo "10 headweight" $headweight
echo "11 basetype" $basetype
echo "12 Urad" $Urad

conv=$(echo "$scale * 0.05" | bc -l | sed 's/^\./0./')
twoblur=$(echo "$scale * 11" | bc -l | sed 's/^\./0./')

if [[ $autofind == "yes" ]]; then
	anat=$(find $procdir -name anat_n[0-9].nii.gz)
	anat=$(dirname $anat)/$(basename -s .nii.gz $anat)
else
	anat=$2
fi
	
#bias correct
N4BiasFieldCorrection -i "$anat".nii.gz -o "$anat"_N4.nii.gz
3dUnifize -prefix "$anat"_Un.nii.gz -Urad $Urad "$anat".nii.gz

if [[ $biascorrector == "ANTSN4" ]]; then Bc=N4; fi
if [[ $biascorrector == "3dUnifize" ]]; then Bc=Un; fi

#brain extraction
thresh=$(3dClipLevel "$anat"_"$Bc".nii.gz) #can be included directly in the next command; keep here as want the value displayed.
printf -v thresh %.0f "$thresh" #RATS cannot handle non-integers
echo "brain/background threshold=$thresh"
#use some fancy tool to identify the brain
RATS_MM "$anat"_"$Bc".nii.gz "$anat"_"$Bc"Bm.nii.gz -v $brainvol -t $thresh
fsl5.0-fslmaths "$anat"_"$Bc".nii.gz -mas "$anat"_"$Bc"Bm.nii.gz "$anat"_"$Bc"BmBe.nii.gz

#the actual T1anat to template registration using the brain extracted image
#could do in one 3dQwarp step using allineate flags but will separate as 3dAllineate performs well on brain image, and 3dQwarp well on whole head
3dAllineate -base $brain -source "$anat"_"$Bc"BmBe.nii.gz -prefix "$anat"_"$Bc"BmBeAl.nii.gz -1Dmatrix_save "$anat"_"$Bc"BmBeAl.aff12.1D -cost nmi -conv $conv -twopass -twoblur $twoblur -cmass -maxrot 90 -master $brain
cat_matvec -ONELINE "$anat"_"$Bc"BmBeAl.aff12.1D -I > "$anat"_"$Bc"BmBeAl_INV.aff12.1D

#application to the whole head image. can also be used for a good demonstration of linear vs. non-linear registration quality
3dAllineate -input "$anat"_"$Bc".nii.gz -master $head -prefix "$anat"_"$Bc"Aa.nii.gz -1Dmatrix_apply "$anat"_"$Bc"BmBeAl.aff12.1D

#non-linear registration of affine pre-registered whole head image to template
#don't initiate straight from the original with an iniwarp due to weird errors (like it creating an Allin it then can't find)

if [[ $basetype == "head" ]]; then
	3dQwarp -base $head -source "$anat"_"$Bc"Aa.nii.gz -prefix "$anat"_"$Bc"AaQw.nii.gz -nmi -iwarp -noneg -weight $headweight -blur 0
elif [[ $basetype == "brain" ]]; then
	3dQwarp -base $brain -source "$anat"_"$Bc"Aa.nii.gz -prefix "$anat"_"$Bc"AaQw.nii.gz -nmi -iwarp -noneg -blur 0
fi

#inverted application of registrations to the atlas
3dNwarpApply -nwarp "$anat"_"$Bc"BmBeAl_INV.aff12.1D "$anat"_"$Bc"AaQw_WARPINV.nii.gz -source $atlas -master "$anat"_"$Bc".nii.gz -ainterp NN -prefix $procdir/atlas_Na1.nii.gz
fixobliquity.bash "$anat"_"$Bc".nii.gz $procdir/atlas_Na1.nii.gz #3dNwarpApply removes obliquity, readd
