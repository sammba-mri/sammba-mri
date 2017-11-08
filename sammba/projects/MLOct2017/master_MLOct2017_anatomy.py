# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 17:26:06 2017

@author: nn241023
"""

#%%

import os
import pandas as pd
import subprocess

#%%

projectdir = '/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouselemur'
projectxlsx = os.path.join(projectdir, 'Tri_data/Trilemur_ok.xlsx')
analysisdir = os.path.join(projectdir, 'MLOct2017/analysis20171030')
templatedir = os.path.join(projectdir, 'RsLemurs/analysis20171017')
head = os.path.join(templatedir, 'Qw4_meanhead_Na_200.nii.gz')
brainmask = os.path.join(templatedir, 'dilated.nii.gz')
atlas = os.path.join(templatedir, 'Lemur-Atlas-Apr2Feb-cortexRL-label2_200.nii.gz')
savedir = os.path.join(projectdir, 'MLOct2017', 'anatomy')

#%%

MRIsessions = pd.read_excel(projectxlsx, sheetname=1)
miniMRIsessions = MRIsessions[(MRIsessions['date'] > '2017-09-26')]
PVdirs = [str(x) for x in miniMRIsessions.PVdir.tolist()]

#%%

for PVdir in PVdirs:
    infotab = miniMRIsessions[(miniMRIsessions.PVdir == PVdir)]
    session = os.path.basename(PVdir)
    NIfTIdir = os.path.join(analysisdir, session)
    try:
        anat_nx = int(infotab.anat.tolist()[0])
    except:
        anat_nx = 0
    sp = os.path.join(NIfTIdir, 'anat_n' + str(anat_nx))
    tp1 = os.path.join(savedir, session + '__brainmask.nii.gz')
    tp2 = os.path.join(savedir, session + '__atlas.nii.gz')
    
    subprocess.call(['3dNwarpApply',
                     '-nwarp',
                         sp + '_UnBmBeAl_INV.aff12.1D',
                         sp + '_UnAaQw_WARPINV.nii.gz',
                     '-source', brainmask,
                     '-master', sp + '_Un.nii.gz',
                     '-ainterp', 'NN',
                     '-prefix', tp1])

    subprocess.call(['3dNwarpApply',
                     '-nwarp',
                         sp + '_UnBmBeAl_INV.aff12.1D',
                         sp + '_UnAaQw_WARPINV.nii.gz',
                     '-source', atlas,
                     '-master', sp + '_Un.nii.gz',
                     '-ainterp', 'NN',
                     '-prefix', tp2])
