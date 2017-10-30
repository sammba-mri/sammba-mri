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
import shutil

startdir = os.path.abspath(os.path.curdir)
#commonPATH = '/mnt/c/Users/nadka/OneDrive/Downloads/git/git/sammba-mri/sammba/common'
commonPATH = '/home/nadkarni/git/sammba-mri/sammba/common'
os.chdir(commonPATH)

import PVEnDCMtoNIfTI
import IDprotocoltypes
import pipelinedefs

os.chdir(startdir)
#afniPATH = '/usr/lib/afni/bin'
#antsPATH = '/usr/lib/ants'
#RATSPATH = '/mnt/c/Users/nadka/OneDrive/Downloads/gittmp/RATS/extracted/distribution'

#os.environ['PATH'] += ':' + afniPATH + ':' + antsPATH + ':' + RATSPATH + ':' + commonPATH
os.environ['PATH'] += ':' + commonPATH
os.environ['AFNI_DECONFLICT'] = 'OVERWRITE'

#%%

projectdir = '/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouselemur'
projectxlsx = os.path.join(projectdir, 'Tri_data/Trilemur_ok.xlsx')
analysisdir = os.path.join(projectdir, 'MLOct2017/analysis20171030')
templatedir = os.path.join(projectdir, 'RsLemurs/analysis20171017')
brain = os.path.join(templatedir, 'Qw4_meanhead_brain_Na_200.nii.gz')
atlas = os.path.join(templatedir, 'Lemur-Atlas-Apr2Feb-cortexRL-label2_200.nii.gz')
mask = os.path.join(templatedir, 'Lemur-Atlas-Apr2Feb-cortexRL-label2_200.nii.gz')
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
    'anat':['_MSME_'],
    'perfFAIREPI':['_Perfusion_FAIR_EPI_'],
    'rs':['_GE_EPI_sat_'],
    'DTIDKI':['_DKI_EPI_'],
    'CESTz':['_zspectrum_']}
dcmdump_path = '/usr/bin/dcmdump'
SIAPfix = 'yes'
logfile = os.path.join(projectdir, 'MLOct2017', strftime("%Y%m%d_%H%M%S", gmtime()) + '_log.txt')
rminterfiles = 'yes'

#%%
def pipeline(PVdir, NIfTIdir, brain, atlas, mask, head, headweight, basetype,
             tmpdir, registerfunctional, Urad, brainvol, scale, T1blood,
             lambda_blood, multiplier, T1guess, mccores, infotab, logfile):

    PVEnDCMtoNIfTI.recursive_EnDCMs_to_NIIs(dcmdump_path, PVdir, NIfTIdir, SIAPfix)
    IDprotocoltypes.ID_protocol_types(protocol_dict, NIfTIdir)
    
    try:
        anat_nx = int(infotab.anat.tolist()[0])
    except:
        anat_nx = 0
    biascorrector = 'Un'

    pipelinedefs.anat_to_template(NIfTIdir, anat_nx, tmpdir, biascorrector,
                                  brainvol, scale, brain, atlas, mask, head,
                                  headweight, basetype, Urad, logfile)
    pipelinedefs.rs_to_template_looper(NIfTIdir, anat_nx, head, tmpdir, scale,
                                       registerfunctional, rminterfiles,
                                       brainvol, logfile)



#%%

if not os.path.exists(analysisdir):
    os.mkdir(analysisdir)
os.mkdir(tmpdir)

#%%

MRIsessions = pd.read_excel(projectxlsx, sheetname=1)
miniMRIsessions = MRIsessions[(MRIsessions['date'] > '2017-09-26')]
PVdirs = [str(x) for x in miniMRIsessions.PVdir.tolist()]

#%%

for PVdir in PVdirs:
    infotab = miniMRIsessions[(miniMRIsessions.PVdir == PVdir)]
    NIfTIdir = os.path.join(analysisdir, os.path.basename(PVdir))
    if not os.path.exists(NIfTIdir):
        os.mkdir(NIfTIdir)
        pipeline(PVdir, NIfTIdir, brain, atlas, mask, head, headweight,
                 basetype, tmpdir, registerfunctional, Urad, brainvol, scale,
                 T1blood, lambda_blood, multiplier, T1guess, mccores, infotab,
                 logfile)

#%%

shutil.rmtree(tmpdir)

