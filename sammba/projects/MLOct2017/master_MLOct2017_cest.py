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

import cesttotemplate

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
head = os.path.join(templatedir, 'Qw4_meanhead_Na_200.nii.gz')
tmpdir = os.path.join('/volatile', 'tmpanaldir_' + str(random.randrange(100000, 1000000)))
registerfunctional = 'no'
Urad = 18.3
brainvol = 1600
scale = 0.2
logfile = os.path.join(projectdir, 'MLOct2017', strftime("%Y%m%d_%H%M%S", gmtime()) + '_log.txt')
rminterfiles = 'yes'

#%%
def pipeline(NIfTIdir, head, tmpdir, scale, registerfunctional, rminterfiles,
             brainvol, Urad, infotab, logfile):
    
    try:
        anat_nx = int(infotab.anat.tolist()[0])
    except:
        anat_nx = 0

    cesttotemplate.cest_to_template_looper(NIfTIdir, anat_nx, head, tmpdir,
                                           scale, registerfunctional,
                                           rminterfiles, brainvol, Urad,
                                           logfile)

#%%

os.mkdir(tmpdir)

#%%

MRIsessions = pd.read_excel(projectxlsx, sheetname=1)
miniMRIsessions = MRIsessions[(MRIsessions['date'] > '2017-09-26')]
PVdirs = [str(x) for x in miniMRIsessions.PVdir.tolist()]

#%%

for PVdir in PVdirs:
    infotab = miniMRIsessions[(miniMRIsessions.PVdir == PVdir)]
    NIfTIdir = os.path.join(analysisdir, os.path.basename(PVdir))
    pipeline(NIfTIdir, head, tmpdir, scale, registerfunctional, rminterfiles,
             brainvol, Urad, infotab, logfile)

#%%

shutil.rmtree(tmpdir)

