=====================================
Data organization
=====================================

.. contents:: **Contents**
    :local:
    :depth: 1


Conversion from Dicom (Bruker-Dicom) to NIfTI-1 data format
===========================================================

    Sammba-MRI provides **easy-to-use** conversion from Dicom (Bruker-Dicom) to NIfTI-1 data format. Two options are available:
	To convert a single file, use the following command:
io_conversions.dcm_to_nii ('/usr/bin/dcmdump', Dicom_input, Saved_Nifty_dir)
	To convert several files, use the following command:
io_conversions.recursive_dcm_to_nii('/usr/bin/dcmdump', session_dir, Saved_Nifty_dir)
	
	
