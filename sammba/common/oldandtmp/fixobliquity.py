# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 17:17:06 2017

@author: nn241023
"""

#%%

import os
import random
import subprocess
import shutil

#https://afni.nimh.nih.gov/afni/community/board/read.php?1,149385,149385#msg-149385
#https://afni.nimh.nih.gov/afni/community/board/read.php?1,74225,74225#msg-74225

#%%

def fix_obliquity(badfile, goodfile, tmpdir, logfile):

#%%

    os.environ["AFNI_DECONFLICT"] = "OVERWRITE"

#%%

    tmpdir = os.path.join(tmpdir, 'tmp_fixobliquity_' + str(random.randrange(100000, 1000000)))
    os.mkdir(tmpdir)    
    goodobliquename = os.path.join(tmpdir, os.path.basename(goodfile)[:-7])
    badobliquename = os.path.join(tmpdir, os.path.basename(badfile)[:-7])

#%%

    sink1 = subprocess.check_output(['3dcopy', goodfile, goodobliquename + '+orig'], stderr=subprocess.STDOUT)
    sink2 = subprocess.check_output(['3dcopy', badfile, badobliquename + '+orig'], stderr=subprocess.STDOUT)

#%%

    sink3 = subprocess.check_output(['3drefit',
        '-atrcopy', goodobliquename + '+orig', 'IJK_TO_DICOM_REAL', badobliquename + '+orig'], stderr=subprocess.STDOUT)
    sink4 = subprocess.check_output(['3dcopy', badobliquename + '+orig', badfile], stderr=subprocess.STDOUT)

#%%

    #logging
    with open(logfile, 'a') as myfile:
        myfile.write("%s\n" % sink1)
        myfile.write("%s\n" % sink2)
        myfile.write("%s\n" % sink3)
        myfile.write("%s\n" % sink4)

    shutil.rmtree(tmpdir)