================================
Measuring regional brain volumes
================================

.. contents:: **Contents**
    :local:
    :depth: 1

	
Sammba-MRI allows to measure regional brain volumes in regions outlined in an atlas. This can be useful to measure cerebral atrophy for example.

Step 1 - Registering several MRI images to an existing template
===============================================================

Sammba-MRI allows to register several MR images to an existing template.
	
The function to use is > anats_to_template(images_to_register, Template_name, Saved_Registered_dir, 400, caching=False)
	
	
Step 2 - Aligning the atlas on the raw MR images thanks to the deformation evaluated in step 1
============================================================================================

So far we did not develop any Sammba procedure to perform this task
We recommend to use the > 3dNwarpApply funtion from AFNI::

    3dNwarpApply -nwarp Saved_Registered_dir/Aal1_WARP.nii.gz Saved_Registered_dir/  
    Aal1_masked_shr.aff12.1D -iwarp -source Atlas_Directory/Atlas.nii.gz -master 
    Images-to-Register/Aal1.nii.gz -ainterp NN -prefix Measures/Aal1_atlas.nii.gz
	
	
Step 3 - Measuring the volume of each registered ROI of the atlas
==================================================================
So far we did not develop any Sammba procedure to perform this task
We recommend to use the > 3dhistog from AFNI::

	3dhistog -int Measures/Aal1_atlas.nii.gz > Measures/Aal1-Measured-volumes.txt

	
