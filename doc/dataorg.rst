=====================================
Data organization
=====================================

.. contents:: **Contents**
    :local:
    :depth: 1


conversion from Dicom (Bruker-Dicom) to NIfTI-1
===============================================

    Sammba-MRI provides **easy-to-use** conversion from Dicom (Bruker-Dicom) to NIfTI-1. Two options are available:
	For a single file > io_conversions.dcm_to_nii ('/usr/bin/dcmdump', Dicom_input, Saved_Nifty_dir)
	For several files > io_conversions.recursive_dcm_to_nii('/usr/bin/dcmdump', session_dir, Saved_Nifty_dir)
	
	