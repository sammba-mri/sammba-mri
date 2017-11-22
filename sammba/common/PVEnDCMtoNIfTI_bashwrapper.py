# -*- coding: utf-8 -*-

#sys: grab shell flags
#PVEnDCMtoNIfTI_funcdefs: DICOM to NIfTI conversion

#%%

import sys
import PVEnDCMtoNIfTI

#%%

dcmdump_path = sys.argv[1]
sessdir = sys.argv[2]
save_directory = sys.argv[3]
SIAPfix = sys.argv[4]

PVEnDCMtoNIfTI.recursive_EnDCMs_to_NIIs(dcmdump_path, sessdir, save_directory, SIAPfix)
