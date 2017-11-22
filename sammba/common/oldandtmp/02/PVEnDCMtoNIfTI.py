# -*- coding: utf-8 -*-

#recursively searches a given folder for EnIm*.dcm files. for each EnIm*.dcm,
#saves a .nii (named using several .dcm tag values) to a(nother) given folder.

#sys: grab shell flags
#os: path formatting, file search
#PVEnDCMtoNIfTI_funcdefs: DICOM to NIfTI conversion

#%%

import sys
import os
import PVEnDCMtoNIfTI_funcdefs as PVtoNII

#%%

dcmdump_path = sys.argv[1]
sessdir = sys.argv[2]
save_directory = sys.argv[3]
SIAPfix = sys.argv[4]

#%%

if os.name == 'nt':
    valstart = 17
    splt1 = '\\\\'
    splt2 = '\\'
if os.name == 'posix':
    valstart = 15
    splt1 = '\\'
    splt2 = '/'

#%%

for root, dirs, files in os.walk(sessdir):
    for file in files:
        if file.startswith('EnIm'):
            if file.endswith('.dcm'):
                PVtoNII.EnDCM_to_NII(dcmdump_path, os.path.join(root, file),
                             save_directory, SIAPfix, valstart, splt1, splt2)
