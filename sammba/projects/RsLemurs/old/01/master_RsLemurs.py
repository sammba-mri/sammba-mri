# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 17:26:06 2017

@author: nn241023
"""

#%%

import os
import random
from time import gmtime, strftime
import pandas as pd

curdir = os.path.abspath(os.curdir)

os.chdir('/home/nadkarni/git/sammba-mri/sammba/common')

import PVEnDCMtoNIfTI
import IDprotocoltypes
import anattotemplate

os.chdir(curdir)

#%%

projectdir = ('/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouselemur')
projectxlsx = os.path.join(projectdir, 'Tri_data/Trilemur_ok.xlsx')
analysisdir = os.path.join(projectdir, 'MLOct2017/analysis20171026')
templatedir = os.path.join(projectdir, 'RsLemurs/analysis20171017')
brain = os.path.join(templatedir, 'Qw4_meanhead_brain_Na_200.nii.gz')
atlas = os.path.join(templatedir, 'Lemur-Atlas-Apr2Feb-cortexR-label_200.nii.gz')
mask = os.path.join(templatedir, 'Lemur-Atlas-Apr2Feb-cortexR-label_200.nii.gz')
head = os.path.join(templatedir, 'Qw4_meanhead_Na_200.nii.gz')
headweight = os.path.join(templatedir, 'aff3_unionmaskdil5_Na_200.nii.gz')
basetype = 'head'
dofolderoverwrite = 'no'
tmpdir = os.path.join('/volatile', 'tmpanaldir_' + str(random.randrange(100000, 1000000)))
registerfunctional = 'no'
Urad = 18.3
brainvol = 1600
scale = 0.2
T1blood = 2800
lambda_blood = 0.9
multiplier = 6000000
T1guess = 1600
mccores = 16
protocol_dict = {
    'anat':['T1_FLASH_3D', 'MSME_MIRCen_allbrain', '__MSME_200um__'],
    'perfFAIREPI':['Perfusion_FAIR_EPI'],
    'rs':['T2star_FID_EPI_sat', 'GE_EPI_sat', 'SE_EPI_sat'],
    'DTIDKI':['DKI_EPI'],
    'CESTz':['zspectrum'],
    'CESTw':['WASSR']}
dcmdump_path = '/usr/bin/dcmdump'
SIAPfix = 'yes'
logfile = os.path.join(projectdir, 'MLOct2017', strftime("%Y%m%d_%H%M%S", gmtime()) + '_log.txt')

#%%
def pipeline(PVdir, NIfTIdir, brain, atlas, mask, head, headweight, basetype,
             tmpdir, registerfunctional, Urad, brainvol, scale, T1blood,
             lambda_blood, multiplier, T1guess, mccores, infotab, logfile):

    PVEnDCMtoNIfTI.recursive_EnDCMs_to_NIIs(dcmdump_path, PVdir, NIfTIdir, SIAPfix)
    IDprotocoltypes.ID_protocol_types(protocol_dict, NIfTIdir)
    anattotemplate.anat_to_template(NIfTIdir, anat_nx, tmpdir, biascorrector,
                                    brainvol, scale, brain, atlas, mask, head,
                                    headweight, basetype, Urad, logfile)
    
    try:
        anat_nx = int(infotab.anat.tolist()[0])
    except:
        anat_nx = 0

    

#%%

if not os.path.exists(analysisdir):
    os.mkdir(analysisdir)
os.mkdir(tmpdir)

#%%

MRIsessions = pd.read_excel(projectxlsx, sheetname=1)
MLAOct2017 = MRIsessions[(MRIsessions['date'] > '2017-09-26')]
PVdirs = [str(x) for x in MLAOct2017.PVdir.tolist()]

#%%

for PVdir in PVdirs:
    infotab = MLAOct2017[(MLAOct2017.PVdir == PVdir)]
    NIfTIdir = os.path.join(analysisdir, os.path.basename(PVdir))
    if not os.path.exists(NIfTIdir):
        os.mkdir(NIfTIdir)
        pipeline(PVdir, NIfTIdir, brain, atlas, mask, head, headweight,
                 basetype, tmpdir, registerfunctional, Urad, brainvol, scale,
                 T1blood, lambda_blood, multiplier, T1guess, mccores, infotab)

#%%
