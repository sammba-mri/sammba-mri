======================================================
Registering several MRI images and creating a template
======================================================

.. contents:: **Contents**
    :local:
    :depth: 1


Registering several raw MRI images to common space
====================================================

The **anats_to_common** pipeline implements the multi-level, iterative scheme 
to create a fine anatomical template from individual anatomical MRI scans. 
First, head images are centered on their respective brain mask centroid while 
intensity biais corrected is performed. A first rough template is obtained by
rigid-body aligning extracted brain to a digitized version of a previous 
histological atlas and applying the transform to the original heads. 
Thereafter, sucessive averaging and registration process is iterated while 
increasing the number of degrees of freedom of the estimated
transform. The target template updated at each step allowing the 
creation of high quality templates.   
    
    In sammba-MRI, you can call the function by typing ::

	from sammba import registration
        anats_to_common(Images_to_register, Saved_Registered_dir, 400, caching=True)
	
Note that most of registration steps relies on the freely available 
`AFNI <https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/background_install/main_toc.html>`_ 
software.
For skull-stripping, the open source `RATS <http://www.iibi.uiowa.edu/content/rats-overview/>`_
tool is also available by setting the parameter ``use_rats_tool`` to True.

An example of use this function is available called *Template creation*

Registering several MRI images to template space
=================================================

Matching of individual MRI image to standard template is critical to perform
group-wise analysis. Multimodal images preprocessing can be performed in the
template space through the *TemplateRegistrator* class. 
User can fit each anatomical images to a a set of modality (functional or 
perfusion MRI). This **ready-to-use** pipelines align first the 
functional or perfusion volume to the anatomical images through a rigid body
registration. Then, a per-slice basis registration is performed allowing
correction of EPI distortion. Finally, the rigid body transform and the 
per-slice warps are combined and applied to the initial volume to minimize
interpolation errors.

    In sammba-MRI, the class can be call by typing ::

        from sammba import registration
        registrator = TemplateRegistrator(template, brain_volume, caching=True)
        registrator.fit_anat(anat_files)
        registrator.fit_modality(moadality_images, name_of_the_modality

An example of use this function is available.

Note that TemplateRegistrator class encapsulate several functions
that are available as "stand-alone" in sammba-MRI:

* anats_to_template : Allow to registering several MRI images to an existing template ::

    anats_to_template(images_to_register, Template_name, Saved_Registered_dir, 400, caching=True)

* fmri_sessions_to_template : Allow to register several functional and anatomic MRI images to an existing template ::

    fmri_sessions_to_template(session, Atlas_name,Saved_Registered_dir, maxlev=7, t_r=1.0, brain_volume=400)

* FMRISession : Encapsulate the data for performing registration of several functional MRI images to an existing template and EPI correction ::

    FMRISession(Anats_to_register, Funcs_to_register, Saved_Registered_dir)

* coregister_fmri_session : process the data for performing registration of several functional MRI images to an existing template and EPI correction ::

    coregister_fmri_session(session, t_r=1.0, Saved_Registered_dir, brain_volume=400, slice_timing=True)


