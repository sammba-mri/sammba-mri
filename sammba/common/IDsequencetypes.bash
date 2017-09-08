#!/bin/bash

NIfTIdir=$(readlink -e $1)

echo "1 NIfTIdir" $NIfTIdir

#for info:
#T2star_FID_EPI 								APPmasitinib 20141126 mouse 742
#20141126_T2star_FID_EPI						APPmasitinib up to 20150909 mouse 660
#20150909_T2star_FID_EPI						APPmasitinib 20150909 mice 663, 664 and 20150910 669, 670, 671, 672 and 680
#20150910_08_T2star_FID_EPI						APPmasitinib 20150910 681, 683, 689 and 690
#06_20150910_T2star_FID_EPI						InductAlz session 1, a few InductAlz session 2, BECIM, early AmylNet
#09_20160422_T2star_FID_EPI_sat					InductAlz session 2, late April-early July AmylNet
#10_20160712_T2star_FID_EPI_sat					mid-July AmylNet
#SE_EPI											23 August AmylNet mouse 37329 and 29 August AmylNet mouse 40658
#11_20160823_SE_EPI_sat_TR1200					23+29 August AmylNet
#11_20160823_SE_EPI_sat_TR2500					23+29 August AmylNet
#12_20160831_SE_EPI_sat_TR2000					31 August AmylNet
#DTI_EPI_30dir									31 August AmylNet
#13_20160831_DTI_EPI_30dir_sat					31 August AmylNet
#12_20161031_SE_EPI_sat_TR2000_10mins			InductAlz session 3
#04_T1_FLASH_3D									5 December 2016 newly formatted protocols
#06_Perfusion_FAIR_EPI							5 December 2016 newly formatted protocols
#07_SE_EPI_sat_TR2000_10mins					5 December 2016 newly formatted protocols
#08_DTI_EPI_30dir_sat							5 December 2016 newly formatted protocols
#09_SE_EPI_sat_TR2000_17min						? December 2016 newly formatted protocol, actually only 15min
#SE_EPI_sat_TR1000_15min_GOP_TR1000_12slices	1 February 2017 test
#GE_EPI_sat_TR1000_15min_GOP_TR1000_12slices	1 February 2017 test
#MSME_MIRCen_allbrain							March DCLK3 study

#using .nii.gz ensures there is only one file for each .nii.gz .txt pair
#use tail -n1 for some sequences as only want one; occasionally more than one
#was acquired. where so, take the final one, which presumably was the best.
a_B0map=$(find $NIfTIdir \( -name "*B0Map-ADJ_B0MAP*.nii.gz" \) | sort | tail -n1)
a_anat=$(find $NIfTIdir \( -name "*T1_FLASH_3D*.nii.gz" -o -name "*MSME*.nii.gz" \) | sort | tail -n1)
a_FcFLASH=$(find $NIfTIdir \( -name "*FcFLASH*.nii.gz" \) | sort | tail -n1)
a_scoutpilot=$(find $NIfTIdir \( -name "*T1_FLASH*.nii.gz" \) | grep -v 3D | sort | tail -n1)
a_T1T2map=$(find $NIfTIdir \( -name "*T1_T2map_RARE*.nii.gz" \) | sort | tail -n1)
a_perfFAIREPI=$(find $NIfTIdir \( -name "*Perfusion_FAIR_EPI*.nii.gz" \) | sort)
a_rs=$(find $NIfTIdir \( -name "*T2star_FID_EPI_sat*.nii.gz" -o -name "*SE_EPI_sat*.nii.gz" -o -name "*GE_EPI_sat*.nii.gz" \) | sort)

a_protocols=(a_B0map a_anat a_FcFLASH a_scoutpilot a_T1T2map a_perfFAIREPI a_rs)
for a_protocol in ${a_protocols[@]}; do
	##VERY IMPORTANT TO NOTE: INDIRECT REFERENCE USED BELOW TO LOOP THROUGH SEQUENCE ARRAYS!!!
	niinames=$(eval "echo \${$a_protocol[@]}") #extract names for a particular sequence protocol
	read -r -a niinames <<< "$niinames" #turn back into an array
	if [ -r "${niinames[0]}" ]; then
		for n in ${!niinames[@]}; do #putting the if inside the for loop allows the n to be generated (used below) regardless of dointernaloverwrite (which we actually got rid of :-D)
			p=$(basename -s .nii.gz ${niinames[$n]})
			mv $NIfTIdir/$p.nii.gz $NIfTIdir/${a_protocol:2}_n$n.nii.gz || true #|| true is because if conversion failed, there is no file to have its name modified. The earlier if statement only checks the first acquisition in the array
			mv $NIfTIdir/"$p"_ptbl.txt $NIfTIdir/${a_protocol:2}_n"$n"_ptbl.txt || true 
			echo "$p changed to ${a_protocol:2}_n$n for both .nii.gz and _ptbl.txt" >> $NIfTIdir/namechanges.txt
		done
	else
		echo "no ${a_protocol:2}"
	fi
done
