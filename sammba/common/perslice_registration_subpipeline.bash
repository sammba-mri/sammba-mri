#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE

workdir=$(readlink -e $1)
anatlab=$2
imfilearray=$3
atlaslab=$4 #not used properly (yet) due to single slice atlas inverse warps often failing (see below)
do4Dtotemplate=$5
template=$(readlink -e $6)
rmfiles=$7
tmpdir=$8
resampleres=$9
dorigidbody=${10}
brainvol=${11}
startdir=$(pwd)

echo "1 workdir" $workdir
echo "2 anatlab" $anatlab
echo "3 imfilearray" ${imfilearray[@]}
echo "4 atlaslab" $atlaslab
echo "5 do4Dtotemplate" $do4Dtotemplate
echo "6 template" $template
echo "7 rmfiles" $rmfiles
echo "8 tmpdir" $tmpdir
echo "9 resampleres" $resampleres
echo "10 dorigidbody" $dorigidbody
echo "11 brainvol" $brainvol

tmpanaldir=tmpanaldir_"$RANDOM"
mkdir $tmpdir/$tmpanaldir
cd $tmpdir/$tmpanaldir

read -r -a imfilearray < $imfilearray
base=${imfilearray[0]}

cp $workdir/"$anatlab".nii.gz .
cp $workdir/"$anatlab"BmBe.nii.gz .
for func in ${imfilearray[@]}; do cp $workdir/"$func".nii.gz .; done
cp $workdir/"$atlaslab".nii.gz .

origanatlab=$anatlab
echo 1 0 0 0 0 1 0 0 0 0 1 0 > ident.1d
cat_matvec -ONELINE ident.1d > ident.aff12.1D
suppanatwarp=ident.aff12.1D

#optional prior whole brain rigid body registration.
if [[ $dorigidbody == "yes" ]]; then
	thresh=$(3dClipLevel $base.nii.gz) #can be included directly in the next command; keep here as want the value displayed.
	printf -v thresh %.0f "$thresh" #RATS cannot handle non-integers
	echo "brain/background threshold=$thresh"
	#use some fancy tool to identify the brain
	RATS_MM $base.nii.gz "$base"_Bm.nii.gz -v $brainvol -t $thresh
	3dcalc -a "$base".nii.gz -b "$base"_Bm.nii.gz -expr "a*b" -prefix "$base"_BmBe.nii.gz
	3dAllineate -base "$anatlab"BmBe.nii.gz -source "$base"_BmBe.nii.gz -prefix "$base"_BmBe_shr.nii.gz -1Dmatrix_save "$base"_BmBe_shr.aff12.1D -cmass -warp shr
	cat_matvec -ONELINE "$base"_BmBe_shr.aff12.1D -I > "$base"_BmBe_shr_INV.aff12.1D
	3dAllineate -input "$anatlab".nii.gz -master $base.nii.gz -prefix "$anatlab"_shr_"$base".nii.gz -1Dmatrix_apply "$base"_BmBe_shr_INV.aff12.1D
	3dAllineate -input "$atlaslab".nii.gz -master $base.nii.gz -prefix "$atlaslab"_shr_"$base".nii.gz -1Dmatrix_apply "$base"_BmBe_shr_INV.aff12.1D #not used later (yet)
	anatlab="$anatlab"_shr_"$base"
	suppanatwarp="$base"_BmBe_shr.aff12.1D
fi

#transform anatomical and atlas to functional space. atlas is already in anatomical space, so only need to record matrix once, from the anatomical
3dWarp -verb -oblique_parent $base.nii.gz -quintic -gridset $base.nii.gz -prefix "$anatlab"_Op_"$base".nii.gz $anatlab.nii.gz > "$anatlab"_Op_"$base"_orig.mat
3dWarp -oblique_parent $base.nii.gz -NN -gridset $base.nii.gz -prefix "$atlaslab"_Op_"$base".nii.gz "$atlaslab".nii.gz #not used later (yet)
cat_matvec -ONELINE "$anatlab"_Op_"$base"_orig.mat > "$anatlab"_Op_"$base".aff12.1D #reformat transform matrix

#invert transform matrix (just in case it's needed, not used for the moment)
cat_matvec -ONELINE "$anatlab"_Op_"$base".aff12.1D -I > "$anatlab"_Op_"$base"_INV.aff12.1D

#3dWarp doesn't put the obliquity in the header, so do it manually
fixobliquity.bash $base.nii.gz "$anatlab"_Op_"$base".nii.gz
fixobliquity.bash $base.nii.gz "$atlaslab"_Op_"$base".nii.gz #not used later (yet)

#slice functional, anatomical and atlas up
fsl5.0-fslslice "$anatlab"_Op_"$base".nii.gz "$anatlab"_Op_"$base"
fsl5.0-fslslice "$atlaslab"_Op_"$base".nii.gz "$atlaslab"_Op_"$base" #not used later (yet)

unset anatlab_Op_base_array atlaslab_array atlaslab_Na_array

for func in ${imfilearray[@]}; do

	fsl5.0-fslslice $func.nii.gz $func
	unset Qw_array Qw_Allin_array func_array Na_array NaOp_array NaNa_array #arrays used to keep track of slices to merge and delete
	for ((q=0; q<$(3dinfo -nk $func.nii.gz); q++)); do #using 3dinfo is wasteful; the result is linked to the original perf file, don't need to remeasure so many times

		func_array=(${func_array[@]} "$func"_slice_"$(printf "%04d\n" $q)".nii.gz)
		Na_array=(${Na_array[@]} "$func"_slice_"$(printf "%04d\n" $q)"_Na.nii.gz)
		NaOp_array=(${NaOp_array[@]} "$func"_slice_"$(printf "%04d\n" $q)"_Na_Op_"$anatlab".nii.gz)
		NaNa_array=(${NaNa_array[@]} "$func"_slice_"$(printf "%04d\n" $q)"_NaNa.nii.gz)
		#below line is to deal with slices where there is no signal (for example rostral end of 20160725_161605_MINDt_AmylNet_64849_1_1 T1anat)
		if [[ $(3dClipLevel "$anatlab"_Op_"$base"_slice_"$(printf "%04d\n" $q)".nii.gz) == 0 || $(3dClipLevel "$func"_slice_"$(printf "%04d\n" $q)".nii.gz) == 0 ]]; then		
		
			cp "$func"_slice_"$(printf "%04d\n" $q)".nii.gz "$func"_slice_"$(printf "%04d\n" $q)"_Na.nii.gz
			#copy fake atlas too
			
		else

			if [[ $func == $base ]]; then			
				
				anatlab_Op_base_array=(${anatlab_Op_base_array[@]} "$anatlab"_Op_"$func"_slice_"$(printf "%04d\n" $q)".nii.gz)
				atlaslab_array=(${atlaslab_array[@]} "$atlaslab"_Op_"$func"_slice_"$(printf "%04d\n" $q)".nii.gz) #not used later (yet)
				atlaslab_Na_array=(${atlaslab_Na_array[@]} "$atlaslab"_Op_"$func"_slice_"$(printf "%04d\n" $q)"_Na.nii.gz)
				Qw_array=(${Qw_array[@]} "$func"_slice_"$(printf "%04d\n" $q)"_Na.nii.gz)
				Qw_Allin_array=(${Qw_Allin_array[@]} "$func"_slice_"$(printf "%04d\n" $q)"_Na_Allin.nii)
				
				#single slice non-linear functional to anatomical registration; the iwarp frequently fails
				#resampling can help it work better
				3dresample -dxyz $resampleres $resampleres $(3dinfo -adk "$anatlab"_Op_"$func"_slice_"$(printf "%04d\n" $q)".nii.gz) -prefix "$anatlab"_Op_"$func"_slice_"$(printf "%04d\n" $q)"_res.nii.gz -inset "$anatlab"_Op_"$func"_slice_"$(printf "%04d\n" $q)".nii.gz
				3dresample -dxyz $resampleres $resampleres $(3dinfo -adk "$func"_slice_"$(printf "%04d\n" $q)".nii.gz) -prefix "$func"_slice_"$(printf "%04d\n" $q)"_res.nii.gz -inset "$func"_slice_"$(printf "%04d\n" $q)".nii.gz
				3dQwarp -base "$anatlab"_Op_"$func"_slice_"$(printf "%04d\n" $q)"_res.nii.gz -source "$func"_slice_"$(printf "%04d\n" $q)"_res.nii.gz -prefix "$func"_slice_"$(printf "%04d\n" $q)"_Na.nii.gz -iwarp -noneg -blur 0 -nmi -noXdis -allineate -allineate_opts '-parfix 1 0 -parfix 2 0 -parfix 3 0 -parfix 4 0 -parfix 5 0 -parfix 6 0 -parfix 7 0 -parfix 9 0 -parfix 10 0 -parfix 12 0'
				3dresample -dxyz $(3dinfo -ad3 "$func"_slice_"$(printf "%04d\n" $q)".nii.gz) -prefix "$func"_slice_"$(printf "%04d\n" $q)"_Na.nii.gz -inset "$func"_slice_"$(printf "%04d\n" $q)"_Na.nii.gz
				rm -f "$anatlab"_Op_"$func"_slice_"$(printf "%04d\n" $q)"_res.nii.gz
				rm -f "$func"_slice_"$(printf "%04d\n" $q)"_res.nii.gz
				fixobliquity.bash "$anatlab"_Op_"$func"_slice_"$(printf "%04d\n" $q)".nii.gz "$func"_slice_"$(printf "%04d\n" $q)"_Na.nii.gz #3dQwarp removes obliquity, readd
				#apply inverse warp to atlas. does not work with either method due to frequent failure to generate inverse map
				#3dNwarpApply -nwarp "$func"_slice_"$(printf "%04d\n" $q)"_Na_WARPINV.nii.gz -source "$atlaslab"_Op_"$func"_slice_"$(printf "%04d\n" $q)".nii.gz -master "$func"_slice_"$(printf "%04d\n" $q)".nii.gz -ainterp NN -prefix "$atlaslab"_Op_"$func"_slice_"$(printf "%04d\n" $q)"_Na.nii.gz
				#3dNwarpApply -nwarp "$func"_slice_"$(printf "%04d\n" $q)"_Na_WARP.nii.gz -iwarp -source "$atlaslab"_Op_"$func"_slice_"$(printf "%04d\n" $q)".nii.gz -master "$func"_slice_"$(printf "%04d\n" $q)".nii.gz -ainterp NN -prefix "$atlaslab"_Op_"$func"_slice_"$(printf "%04d\n" $q)"_Na.nii.gz

			else
				
				#apply warp to functional (unnecessary for the first functional which is done above but can act as a paranoia check: should be an exact copy
				3dNwarpApply -nwarp "$base"_slice_"$(printf "%04d\n" $q)"_Na_WARP.nii.gz -source "$func"_slice_"$(printf "%04d\n" $q)".nii.gz -master "$func"_slice_"$(printf "%04d\n" $q)".nii.gz -prefix "$func"_slice_"$(printf "%04d\n" $q)"_Na.nii.gz
				fixobliquity.bash "$func"_slice_"$(printf "%04d\n" $q)".nii.gz "$func"_slice_"$(printf "%04d\n" $q)"_Na.nii.gz #3dNwarpApply removes obliquity, readd
				
			fi

		fi

		if [[ "$(3dinfo -nv $func.nii.gz)" -le 1 || "$do4Dtotemplate" == "yes" ]]; then
			#into template space, both linear (commented out) and non-linear options
			#3dAllineate -input "$func"_Na_Op_"$anatlab"_MeDoReSd.nii.gz -master $template -prefix "$func"_Na_Op_"$anatlab"_MeDoReSdAa.nii.gz -1Dmatrix_apply "$anatlab"_DoReSdBmBeAl.aff12.1D
			3dWarp -oblique_parent "$anatlab".nii.gz -gridset "$anatlab".nii.gz -prefix "$func"_slice_"$(printf "%04d\n" $q)"_Na_Op_"$anatlab".nii.gz "$func"_slice_"$(printf "%04d\n" $q)"_Na.nii.gz
			fixobliquity.bash "$anatlab".nii.gz "$func"_slice_"$(printf "%04d\n" $q)"_Na_Op_"$anatlab".nii.gz
			#to avoid the double-step, maybe cat_matvec the two aff12.1Ds which might avoid the need to intervene with a fixobliquity
			#-ainterp linear (or maybe NN is intellectually neatest, though ugly) avoid slice boundary artefacts during 3dMean below
			3dNwarpApply -nwarp $workdir/"$origanatlab"AaQw_WARP.nii.gz $workdir/"$origanatlab"BmBeAl.aff12.1D $suppanatwarp -source "$func"_slice_"$(printf "%04d\n" $q)"_Na_Op_"$anatlab".nii.gz -master $template -ainterp linear -prefix "$func"_slice_"$(printf "%04d\n" $q)"_NaNa.nii.gz		
		fi
		
	done

	fsl5.0-fslmerge -z "$func"_NaMe.nii.gz $(echo ${Na_array[@]})
	#below line is pointless as we don't have a line anywhere applying the (currently impossible-to-produce) inverse warps to the atlases 
	#fsl5.0-fslmerge -z "$func"_NaMe.nii.gz $(echo ${atlaslab_array[@]})
	
	if [[ "$(3dinfo -nv $func.nii.gz)" -le 1 || "$do4Dtotemplate" == "yes" ]]; then
		#my head hurts thinking about whether -sum and/or -non_zero are good or not, and their interaction with the interpolation of source images
		#deliberately not using -non_zero in order to make sure the interpolation is accurate (if not, slice boundaries will be bright)
		3dMean -prefix "$func"_NaNaMe.nii.gz -sum $(echo ${NaNa_array[@]})
	fi
	
	if [[ $rmfiles == "yes" ]]; then
		rm -f $(echo ${atlaslab_array[@]})
		rm -f $(echo ${Qw_array[@]})
		rm -f $(echo ${Qw_Allin_array[@]})
		rm -f $(echo ${atlaslab_Na_array[@]})
		rm -f $(echo ${func_array[@]})
		rm -f $(echo ${Na_array[@]})
		rm -f $(echo ${NaOp_array[@]})
		rm -f $(echo ${NaNa_array[@]})
	fi
	
	unset Qw_array Qw_Allin_array func_array Na_array NaOp_array NaNa_array
done

if [[ $rmfiles == "yes" ]]; then
	rm -f $(echo ${anatlab_Op_base_array[@]})
	rm -f "$anatlab".nii.gz
	for func in ${imfilearray[@]}; do rm -f "$func".nii.gz; done
	rm -f "$atlaslab".nii.gz
fi

unset anatlab_Op_base_array atlaslab_array atlaslab_Na_array

cp -f * $workdir
cd $startdir
rm -rfd $tmpdir/$tmpanaldir
