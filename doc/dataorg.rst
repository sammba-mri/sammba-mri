=====================================
Data conversion
=====================================

.. contents:: **Contents**
    :local:
    :depth: 1


Data conversion: from DICOM (Bruker-Dicom) to NIfTI-1 format
=================================================================

Sammba-MRI allows to convert Bruker Paravision enhanced multiframe DICOM files 
to the standard NIfTI-1 format. Extensive information such as date is also 
extracted. Depends on the number of file to be converted, two options are 
available:

* To convert a **single** file, we use ::

    io_conversions.dcm_to_nii('/usr/bin/dcmdump', Dicom_input, Saved_Nifti_dir)

* To convert **several** files, we use ::

    io_conversions.recursive_dcm_to_nii('/usr/bin/dcmdump', session_dir, Saved_Nifti_dir)

