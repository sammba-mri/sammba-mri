#!/bin/bash

export AFNI_DECONFLICT=OVERWRITE #saves a lot of grief

dicomdirlist=$(readlink -e $1)
savedir=$(readlink -e $2)

nlines=$(wc -l $dicomdirlist | awk '{print $1}')

for ((a=1; a<=$nlines; a++)); do
	dicomdir=$(head -n$a < $dicomdirlist | tail -n1)
	NIfTIdir=$savedir/$(basename $dicomdir)
	mkdir $NIfTIdir
	python -m PVEnDCMtoNIfTI.py /usr/bin/dcmdump $dicomdir $NIfTIdir yes
	IDsequencetypes.bash $NIfTIdir
done	

mv $savedir/20170425_115121_208CBF_1_1/anat_n0.nii.gz $savedir/20170425_115121_208CBF_1_1/badanat_n1.nii.gz
mv $savedir/20170425_115121_208CBF_1_1/208CBF__20170425__132441__MSME_MIRCen_RatBrain_atlas__fixedSIAP__bf11.nii.gz $savedir/20170425_115121_208CBF_1_1/anat_n0.nii.gz
mv $savedir/20170425_164314_289BB_1_1/anat_n0.nii.gz $savedir/20170425_164314_289BB_1_1/badanat_n2.nii.gz
mv $savedir/20170425_164314_289BB_1_1/289BB__20170425__171045__MSME_200um__fixedSIAP__bf6.nii.gz $savedir/20170425_164314_289BB_1_1/anat_n0.nii.gz
mv $savedir/20170425_192223_310C_1_1/anat_n0.nii.gz $savedir/20170425_192223_310C_1_1/badanat_n1.nii.gz
mv $savedir/20170425_192223_310C_1_1/310C__20170425__194739__MSME_200um__fixedSIAP__bf6.nii.gz $savedir/20170425_192223_310C_1_1/anat_n0.nii.gz

for ((a=1; a<=$nlines; a++)); do
	dicomdir=$(head -n$a < $dicomdirlist | tail -n1)
	NIfTIdir=$savedir/$(basename $dicomdir)
	3dZeropad -RL 160 -AP 160 -IS 100 -prefix $NIfTIdir/anat.nii.gz $NIfTIdir/anat_n0.nii.gz
	rm -f $(find $NIfTIdir -mindepth 1 ! -name anat.nii.gz)
done

#create raw data video and mean
3dTcat -prefix $savedir/raw_video.nii.gz $(find $savedir -name 'anat.nii.gz')
3dTstat -prefix $savedir/raw_mean.nii.gz $savedir/raw_video.nii.gz
