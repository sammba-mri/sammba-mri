# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 17:26:06 2017

@author: nn241023
"""

#%%

import os
import pandas as pd
import shutil

#%%

projectdir = '/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouselemur'
projectxlsx = os.path.join(projectdir, 'Tri_data/Trilemur_ok.xlsx')
analysisdir = os.path.join(projectdir, 'MLOct2017')

#%%

IDs = pd.read_excel(projectxlsx, sheetname=0)
MRIsessions = pd.read_excel(projectxlsx, sheetname=1)
alldata = pd.merge(IDs, MRIsessions)

old = alldata[(alldata['date'] > '2017-09-26') & (alldata['DOB'] < '2010-01-01')]
young = alldata[(alldata['date'] > '2017-09-26') & (alldata['DOB'] > '2010-01-01')]

oldPVdirs = [str(x) for x in old.PVdir.tolist()]
youngPVdirs = [str(x) for x in young.PVdir.tolist()]

#%%

for PVdir in oldPVdirs:
    session = os.path.basename(PVdir)
    NIfTIdir = os.path.join(analysisdir, 'analysis20171030', session)
    CESTz = 'CESTz_n0_C1PxMn_NaNaMe.nii.gz'
    src = os.path.join(NIfTIdir, CESTz)
    dest = os.path.join(analysisdir, 'CEST', 'old', session + '__' + CESTz)
    shutil.copy(src, dest) 

#%%

for PVdir in youngPVdirs:
    session = os.path.basename(PVdir)
    NIfTIdir = os.path.join(analysisdir, 'analysis20171030', session)
    CESTz = 'CESTz_n0_C1PxMn_NaNaMe.nii.gz'
    src = os.path.join(NIfTIdir, CESTz)
    dest = os.path.join(analysisdir, 'CEST', 'young', session + '__' + CESTz)
    shutil.copy(src, dest) 
