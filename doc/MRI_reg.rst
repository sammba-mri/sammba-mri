======================================================
Registering several MRI images and creating a template
======================================================

.. contents:: **Contents**
    :local:
    :depth: 1


Registering several raw MRI images
==============================

    Sammba-MRI allows to register several MR images.
		° First it creates a template from multiple anatomical scans, then register all of them to it. 
		° Initially, registration is of extracted brains. Once these are reasonably aligned, 
		° whole heads are registered, weighted by masks that, if parameters are chosen well, include some scalp. 
		The amount of scalp is hopefully enough to help in differentiating the brain-scalp boundary without including so much head tissue 
		that it starts to contaminate the registration with the highly variable head tissue.

	The function to use is > anats_to_common(Images_to_register, Saved_Registered_dir, 400, caching=True)
	

Registering several MRI images to an existing template
======================================================

    Sammba-MRI allows to register several MR images to an existing template.
	
	The function to use is > anats_to_template(images_to_register, Template_name, Saved_Registered_dir, 400, caching=True)
	
Registering several functional and anatomic MRI images to an existing template
==============================================================================

    The function to use is > fmri_sessions_to_template(session, Atlas_name,Saved_Registered_dir, maxlev=7, t_r=1.0, brain_volume=400)
	
Registering several functional MRI images to an existing template and EPI correction
=====================================================================================

	Sammba allows to align the functional volume to the anatomical images, first with a rigid body registration 
	and then on a per-slice basis (only a fine correction, this is mostly for correction of EPI distortion).

	The first function to use to encapsulate the data is > FMRISession(Anats_to_register, Funcs_to_register, Saved_Registered_dir)
	The second function to use to process the data is > coregister_fmri_session(session, t_r=1.0, Saved_Registered_dir, brain_volume=400, slice_timing=True)

	