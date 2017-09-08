cd /home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouselemur/analysis20170825

3dTcat -prefix ../booyahone.nii.gz $(cat ../listone.txt)

cd ..

fsl5.0-fslroi booyahone.nii.gz booyahone_x73.nii.gz 73 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim booyahone_x73.nii.gz -z y x booyahone_x73.nii.gz
gunzip booyahone_x73.nii.gz

fsl5.0-fslroi booyahone.nii.gz booyahone_y73.nii.gz 0 -1 73 1 0 -1 0 -1
fsl5.0-fslswapdim booyahone_y73.nii.gz x -z y booyahone_y73.nii.gz
gunzip booyahone_y73.nii.gz

fsl5.0-fslroi booyahone.nii.gz booyahone_z80.nii.gz 0 -1 0 -1 80 1 0 -1
gunzip booyahone_z80.nii.gz

cd analysis20170825

3dTcat -prefix ../raw.nii.gz $(cat ../listoneraw.txt)

cd ..

fsl5.0-fslroi raw.nii.gz raw_x46.nii.gz 46 1 0 -1 0 -1 0 -1
fsl5.0-fslswapdim raw_x46.nii.gz -z y x raw_x46.nii.gz
gunzip raw_x46.nii.gz

fsl5.0-fslroi raw.nii.gz raw_y40.nii.gz 0 -1 40 1 0 -1 0 -1
fsl5.0-fslswapdim raw_y40.nii.gz x -z y raw_y40.nii.gz
gunzip raw_y40.nii.gz

fsl5.0-fslroi raw.nii.gz raw_z15.nii.gz 0 -1 0 -1 15 1 0 -1
gunzip raw_z15.nii.gz
