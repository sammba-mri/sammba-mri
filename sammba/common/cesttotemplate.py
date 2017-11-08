# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 14:55:30 2017

@author: nadka
"""
#%%

import os
import subprocess
import glob

#%%

#get a proper pythoner to check this, am sure it's a disaster waiting to happen
def spco(*popenargs, **kwargs):
    x = subprocess.check_output(stderr=subprocess.STDOUT, *popenargs, **kwargs)
    return x


#%%

def cest_to_template(cestz_nx, cestw_nx, zref, anat_nx, template, tmpdir,
                     scale, registerfunctional, rminterfiles, brainvol, Urad,
                     logfile):

#%%

    cestz_nx = cestz_nx[:-7]
    cestz_nx_C1 = cestz_nx + '_C1'
    sink1 = spco(['3dcalc',
                  '-float',
                  '-a', cestz_nx + '.nii.gz',
                  '-b', cestz_nx + '.nii.gz[0]',
                  '-expr', 'a/b',
                  '-prefix', cestz_nx_C1 + '.nii.gz'])
    
#%%    

    cestz_nx_C1P5 = cestz_nx_C1 + 'P5'
    sink2 = spco(['3dcalc',
                  '-a', cestz_nx_C1 + '.nii.gz[0]',
                  '-b', cestz_nx_C1 + '.nii.gz[10]',
                  '-expr', 'a-b',
                  '-prefix', cestz_nx_C1P5 + '.nii.gz'])
    cestz_nx_C1P4 = cestz_nx_C1 + 'P4'
    sink3 = spco(['3dcalc',
                  '-a', cestz_nx_C1 + '.nii.gz[1]',
                  '-b', cestz_nx_C1 + '.nii.gz[9]',
                  '-expr', 'a-b',
                  '-prefix', cestz_nx_C1P4 + '.nii.gz'])
    cestz_nx_C1P3 = cestz_nx_C1 + 'P3'
    sink4 = spco(['3dcalc',
                  '-a', cestz_nx_C1 + '.nii.gz[2]',
                  '-b', cestz_nx_C1 + '.nii.gz[8]',
                  '-expr', 'a-b',
                  '-prefix', cestz_nx_C1P3 + '.nii.gz'])
    cestz_nx_C1P2 = cestz_nx_C1 + 'P2'
    sink5 = spco(['3dcalc',
                  '-a', cestz_nx_C1 + '.nii.gz[3]',
                  '-b', cestz_nx_C1 + '.nii.gz[7]',
                  '-expr', 'a-b',
                  '-prefix', cestz_nx_C1P2 + '.nii.gz'])
    cestz_nx_C1P1 = cestz_nx_C1 + 'P1'
    sink6 = spco(['3dcalc',
                  '-a', cestz_nx_C1 + '.nii.gz[4]',
                  '-b', cestz_nx_C1 + '.nii.gz[6]',
                  '-expr', 'a-b',
                  '-prefix', cestz_nx_C1P1 + '.nii.gz'])
    cestz_nx_C1P0 = cestz_nx_C1 + 'P0'
    sink7 = spco(['3dcalc',
                  '-a', cestz_nx_C1 + '.nii.gz[5]',
                  '-b', cestz_nx_C1 + '.nii.gz[5]',
                  '-expr', 'a-b',
                  '-prefix', cestz_nx_C1P0 + '.nii.gz'])
    
#%%

    cest_nx_C1PxMn = cestz_nx_C1 + 'PxMn'
    sink8 = spco(['3dMean',
                  '-prefix', cest_nx_C1PxMn + '.nii.gz',
                  cestz_nx_C1P3 + '.nii.gz', cestz_nx_C1P2 + '.nii.gz'])

#%%

    cestz_nx_Un = cestz_nx + '_Un'
    sink9 = spco(['3dUnifize',
                  '-prefix', cestz_nx_Un + '.nii.gz',
                  '-Urad', str(Urad),
                  cestz_nx + '.nii.gz[0]'])
    
#%%    

    #logging
    with open(logfile, 'a') as myfile:
        myfile.write("%s\n" % sink1)
        myfile.write("%s\n" % sink2)
        myfile.write("%s\n" % sink3)
        myfile.write("%s\n" % sink4)
        myfile.write("%s\n" % sink5)
        myfile.write("%s\n" % sink6)
        myfile.write("%s\n" % sink7)
        myfile.write("%s\n" % sink8)
        myfile.write("%s\n" % sink9)

#%%

    #bash script
    br = os.path.basename(cestz_nx)
    imfilearray = cestz_nx + '_imfilearray.txt'
    with open(imfilearray, 'w') as myfile:
        myfile.write(br + '_Un ' + br + '_C1PxMn')

    startdir = os.path.abspath(os.path.curdir)
    sink10 = spco(['perslice_registration_subpipeline.bash',
                   os.path.dirname(cestz_nx), 'anat_n' + str(anat_nx) + '_Un',
                   imfilearray, 'atlas_Na1', 'yes', template, rminterfiles,
                   tmpdir, str(scale), registerfunctional, str(brainvol),
                   startdir])

#%%

    with open(logfile, 'a') as myfile:
        myfile.write("%s\n" % sink10)

#%%

def cest_to_template_looper(NIfTIdir, anat_nx, template, tmpdir, scale,
                            registerfunctional, rminterfiles, brainvol, Urad,
                            logfile):

    all_cestz_nx = glob.glob(os.path.join(NIfTIdir, 'CESTz_n[0-9].nii.gz'))
    cestw_nx = 'none'
    zref = 'none'

    for cestz_nx in all_cestz_nx:
        cest_to_template(cestz_nx, cestw_nx, zref, anat_nx, template, tmpdir,
                         scale, registerfunctional, rminterfiles, brainvol,
                         Urad, logfile)
