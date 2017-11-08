#!/bin/bash

#bash -xf /home/nadkarni/git/sammba-mri/sammba/projects/MIRCen_sammba-mri_mouse/master_MIRCen_sammba-mri_mouse.bash 2>&1 | tee /home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MIRCen_sammba-mri_mouse/$(date +%Y%m%d_%H%M%S).log

export AFNI_DECONFLICT=OVERWRITE

PATH=$PATH:/home/nadkarni/git/sammba-mri/sammba/common
export PATH

PYTHONPATH=$PYTHONPATH:$PATH
export PYTHONPATH

projectdir=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MIRCen_sammba-mri_mouse
savedir=$projectdir/analysis20171012
mkdir $savedir

dicomdirlist=$projectdir/dicomdirs.txt
nlines=$(wc -l $dicomdirlist | awk '{print $1}')

for ((a=1; a<=$nlines; a++)); do
	dicomdir=$(head -n$a < $dicomdirlist | tail -n1)
	NIfTIdir=$savedir/$(basename $dicomdir)
	mkdir $NIfTIdir
	python -m PVEnDCMtoNIfTI.py /usr/bin/dcmdump $dicomdir $NIfTIdir yes
	IDsequencetypes.bash $NIfTIdir
done

for ((a=1; a<=$nlines; a++)); do
	dicomdir=$(head -n$a < $dicomdirlist | tail -n1)
	NIfTIdir=$savedir/$(basename $dicomdir)
	3dZeropad -RL 150 -AP 200 -IS 100 -prefix $NIfTIdir/anat.nii.gz $NIfTIdir/anat_n0.nii.gz
done

#create raw data video and mean
3dTcat -prefix $savedir/raw_video.nii.gz $(find $savedir -name 'anat.nii.gz')
3dTstat -prefix $savedir/raw_mean.nii.gz $savedir/raw_video.nii.gz

finres=0.1
conv=0.005
twoblur=1
brainvol=400
#-rbt values might be improveable
Urad=18.3
b=70
t=80

bash MRIT2_extrcen.bash $savedir $brainvol $Urad $b $t
bash MRIT3_shr.bash $savedir $savedir/UnBmBeCC_mean.nii.gz 1 $conv $twoblur
bash MRIT3_shr.bash $savedir $savedir/shr1_mean.nii.gz 2 $conv $twoblur
bash MRIT4_aff.bash $savedir $savedir/shr2_meanhead.nii.gz 2 3 $conv $twoblur $savedir/shr2_count.nii.gz

3dmask_tool -frac 0.5 -inputs $savedir/aff3_video.nii.gz -prefix $savedir/aff3_frac05mask.nii.gz
3dcalc -a $savedir/aff3_meanhead.nii.gz -expr 'step(0.5-(x+3.95)*(x+3.95)-(y-3.75)*(y-3.75)-(z+0.75)*(z+0.75))' -prefix $savedir/rightcolliculus.nii.gz
3dcalc -a $savedir/aff3_meanhead.nii.gz -expr 'step(0.5-(x-3.75)*(x-3.75)-(y-3.75)*(y-3.75)-(z+0.85)*(z+0.85))' -prefix $savedir/leftcolliculus.nii.gz
3dmask_tool -union -inputs $savedir/aff3_frac05mask.nii.gz $savedir/rightcolliculus.nii.gz $savedir/leftcolliculus.nii.gz -prefix $savedir/aff3_unionmask.nii.gz
3dmask_tool -dilate_inputs 5 -inputs $savedir/aff3_unionmask.nii.gz -prefix $savedir/aff3_unionmaskdil5.nii.gz
weight=$savedir/aff3_unionmaskdil5.nii.gz

bash MRIT5_Qw.bash $savedir $savedir/aff3_meanhead.nii.gz $weight UnCC.nii.gz UnBmBeCCAl2UnCCAl3.aff12.1D 0 4 1
bash MRIT5_Qw.bash $savedir $savedir/Qw1_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw1_WARP.nii.gz 5 7 2
bash MRIT6_Qw.bash $savedir $savedir/Qw2_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw2_WARP.nii.gz 8 13 3
bash MRIT6_Qw.bash $savedir $savedir/Qw3_meanhead.nii.gz $weight UnCC.nii.gz UnCCQw3_WARP.nii.gz 11 9 4

MIC40C57=/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/MRIatlases/MIC40C57/20170223
dorrbrain=$MIC40C57/brain100.nii.gz
dorratlas=$MIC40C57/labels100.nii.gz
3dAllineate -source $dorrbrain -base $savedir/aff3_mean.nii.gz -prefix $savedir/dorrbrainAl.nii.gz -1Dmatrix_save $savedir/dorrbrainAl.aff12.1D
3dQwarp -base $savedir/Qw4_meanhead.nii.gz -source $savedir/dorrbrainAl.nii.gz -prefix $savedir/dorrbrainAlQw.nii.gz -iwarp -noneg -blur 0 -nmi
cat_matvec -ONELINE $savedir/dorrbrainAl.aff12.1D -I > $savedir/dorrbrainAl_INV.aff12.1D
3dNwarpCat -prefix $savedir/dorrbrainAlQwcat_WARP.nii.gz $savedir/dorrbrainAlQw_WARP.nii.gz $savedir/dorrbrainAl.aff12.1D
3dNwarpCat -prefix $savedir/dorrbrainAlQwcat_WARPINV.nii.gz $savedir/dorrbrainAl_INV.aff12.1D $savedir/dorrbrainAlQw_WARPINV.nii.gz
cp $MIC40C57/labels100.nii.gz $savedir
3dNwarpApply -nwarp $savedir/dorrbrainAlQwcat_WARP.nii.gz -source $MIC40C57/labels100.nii.gz -master $savedir/Qw4_meanhead.nii.gz -ainterp NN -prefix $savedir/dorratlasAlQw.nii.gz
3dNwarpApply -nwarp $savedir/dorrbrainAlQwcat_WARPINV.nii.gz -source $savedir/Qw4_meanhead.nii.gz -master $MIC40C57/brain100.nii.gz -prefix $savedir/Qw4_meanhead_Na.nii.gz
3dNwarpApply -nwarp $savedir/dorrbrainAlQwcat_WARPINV.nii.gz -source $savedir/aff3_unionmaskdil5.nii.gz -master $MIC40C57/brain100.nii.gz -prefix $savedir/aff3_unionmaskdil5_Na.nii.gz
3dresample -dxyz 0.2 0.2 0.2 -prefix $savedir/Qw4_meanhead_Na_200.nii.gz -input $savedir/Qw4_meanhead_Na.nii.gz
3dcalc -a $savedir/Qw4_meanhead_Na.nii.gz -b $savedir/labels100.nii.gz -expr 'a*ispositive(b)' -prefix $savedir/Qw4_meanhead_Na_brain.nii.gz

brain=$savedir/Qw4_meanhead_Na_brain.nii.gz #4
atlas=$savedir/labels100.nii.gz #5
mask=$savedir/labels100.nii.gz #6
head=$savedir/Qw4_meanhead_Na.nii.gz #7
headweight=$savedir/aff3_unionmaskdil5_Na.nii.gz #8
basetype=head #9
dofolderoverwrite=no #10
tmpdir=/volatile #11
registerfunctional=no #12
subpipeline=perslice_registration_subpipeline.bash #13
Urad=18.3 #14
brainvol=400 #15
scale=0.1 #16
T1blood=2800 #17
lambda=0.9 #18
multiplier=6000000 #19
T1guess=1600 #20
mccores=16 #21

rawimdir=NA

for dir in $(find $savedir -mindepth 1 -maxdepth 1 -type d); do
	NIfTIdir=$(readlink -e $dir)
	pipeline_noconversion.bash $rawimdir $NIfTIdir $brain $atlas $mask $head $headweight $basetype $tmpdir $registerfunctional $subpipeline $Urad $brainvol $scale
done
