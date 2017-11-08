# -*- coding: utf-8 -*-

#%%

import sys
import os

#%%

sessdir = sys.argv[1]

#%%

filelist = []

d={'anat':['T1_FLASH_3D', 'MSME_MIRCen_allbrain', '__MSME_200um__'],
   'perfFAIREPI':['Perfusion_FAIR_EPI'],
   'rs':['T2star_FID_EPI_sat', 'GE_EPI_sat', 'SE_EPI_sat'],
   'DTIDKI':['DKI_EPI'],
   'CESTz':['zspectrum'],
   'CESTw':['WASSR']}

changelog = []

#%%

for root, dirs, files in os.walk(sessdir):
    for filename in files:
        if filename.endswith('.nii.gz'):
            filelist.append(filename)

#%%

for key, value in d.items():
    minifilelist = []
    for file in filelist:
        if any(s in file for s in value):
            minifilelist.append(file)
    if len(minifilelist) > 0:
        minifilelist.sort()
        for n in range(0, len(minifilelist)):
            nosuffix = minifilelist[n].rstrip('.nii.gz')
            oldbasepath = os.path.join(sessdir, nosuffix)
            newbasepath = os.path.join(sessdir, key) + '_n' + str(n)
            os.rename(os.path.join(oldbasepath + '.nii.gz'),
                      os.path.join(newbasepath + '.nii.gz'))
            os.rename(os.path.join(oldbasepath + '_ptbl.txt'),
                      os.path.join(newbasepath + '_ptbl.txt'))
            changelog.append(nosuffix + ' to ' + key + '_n' + str(n))
    else:
        changelog.append('no files converted to ' + key)

#%%

f = open (os.path.join(sessdir, 'namechanges.txt'), 'w')
for line in changelog:
    f.write("%s\n" % line)
f.close()
