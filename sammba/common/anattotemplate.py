# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 19:51:28 2017

@author: nn241023
"""

#%%

import os
import subprocess
import fixobliquity

#get a proper pythoner to check this, am sure it's a disaster waiting to happen
def spco(*popenargs, **kwargs):
    x = subprocess.check_output(stderr=subprocess.STDOUT, *popenargs, **kwargs)
    return x

#%%

def anat_to_template(NIfTIdir, anat_nx, tmpdir, biascorrector, brainvol, scale,
                     brain, atlas, mask, head, headweight, basetype, Urad,
                     logfile):

#%%
   
    conv = str(scale * 0.05)
    twoblur = str(scale * 11)
    
    anatstem = os.path.join(NIfTIdir, 'anat_n' + str(anat_nx))

#%%
    
    #bias correct
    if biascorrector == 'N4':
        anatstemBc = anatstem + '_N4'
        sink1 = spco(['N4BiasFieldCorrection',
                      '-i', anatstem + '.nii.gz',
                      '-o', anatstemBc + '.nii.gz' ])
    elif biascorrector == 'Un':
        anatstemBc = anatstem + '_Un'
        sink1 = spco(['3dUnifize',
                      '-prefix', anatstemBc + '.nii.gz',
                      '-Urad', str(Urad),
                      anatstem + '.nii.gz' ])

#%%

    #brain extraction
    sink2 = spco(['3dClipLevel', anatstemBc + '.nii.gz'],
                 stderr=subprocess.STDOUT)
    anatstemBcBm = anatstemBc + 'Bm'
    anatstemBcBmBe = anatstemBc + 'BmBe'
    
    #identify the brain. RATS cannot handle non-integers
    sink3 = spco(['RATS_MM',
                  anatstemBc + '.nii.gz',
                  anatstemBcBm + '.nii.gz',
                  '-v', str(int(round(brainvol, 0))),
                  '-t', str(int(round(float(sink2.splitlines()[-1]), 0)))])
    sink4 = spco(['3dcalc',
                  '-a', anatstemBc + '.nii.gz',
                  '-b', anatstemBcBm + '.nii.gz',
                  '-expr', 'a*b',
                  '-prefix', anatstemBcBmBe + '.nii.gz'])

#%%

    #logging
    with open(logfile, 'a') as myfile:
        myfile.write("%s\n" % sink1)
        myfile.write("%s\n" % sink2)
        myfile.write("%s\n" % sink3)
        myfile.write("%s\n" % sink4)

#%%

    #the actual T1anat to template registration using the brain extracted image
    #could do in one 3dQwarp step using allineate flags but will separate as
    #3dAllineate performs well on brain image, and 3dQwarp well on whole head
    anatstemBcBmBeAl = anatstemBcBmBe + 'Al'
    anatstemBcAa = anatstemBc + 'Aa'
    sink5 = spco(['3dAllineate',
                  '-base', brain,
                  '-source', anatstemBcBmBe + '.nii.gz',
                  '-prefix', anatstemBcBmBeAl + '.nii.gz',
                  '-1Dmatrix_save', anatstemBcBmBeAl + '.aff12.1D',
                  '-cost', 'nmi', '-conv', conv, '-twopass',
                  '-twoblur', twoblur, '-cmass', '-maxrot', '90',
                  '-master', brain])
    AlInvmat = spco(['cat_matvec',
                     '-ONELINE', anatstemBcBmBeAl + '.aff12.1D', '-I'])
    with open(anatstemBcBmBeAl + '_INV.aff12.1D', 'w') as myfile:
        myfile.write(AlInvmat)
    #application to the whole head image. can also be used for a good
    #demonstration of linear vs. non-linear registration quality
    sink6 = spco(['3dAllineate',
                  '-input', anatstemBc + '.nii.gz',
                  '-master', head,
                  '-prefix', anatstemBcAa + '.nii.gz',
                  '-1Dmatrix_apply', anatstemBcBmBeAl + '.aff12.1D'])

#%%

    #logging
    with open(logfile, 'a') as myfile:
        myfile.write("%s\n" % sink5)
        myfile.write("%s\n" % sink6)

#%%

    #non-linear registration of affine pre-registered whole head image to template
    #don't initiate straight from the original with an iniwarp due to weird errors
    #(like it creating an Allin it then can't find)
    anatstemBcAaQw = anatstemBcAa + 'Qw'
    if basetype == 'head':
        sink7 = spco(['3dQwarp',
                      '-base', head,
                      '-source', anatstemBcAa + '.nii.gz',
                      '-prefix', anatstemBcAaQw + '.nii.gz',
                      '-nmi', '-iwarp', '-noneg',
                      '-weight', headweight,
                      '-blur', '0'])
    elif basetype == 'brain':
        sink7 = spco(['3dQwarp',
                      '-base', brain,
                      '-source', anatstemBcAa + '.nii.gz',
                      '-prefix', anatstemBcAaQw + '.nii.gz',
                      '-nmi', '-iwarp', '-noneg', '-blur', '0'])
    #inverted application of registrations to the atlas
    atlas_Na1 = os.path.join(NIfTIdir, 'atlas_Na1.nii.gz')
    sink8 = spco(['3dNwarpApply',
                  '-nwarp', anatstemBcBmBeAl + '_INV.aff12.1D',
                            anatstemBcAaQw + '_WARPINV.nii.gz',
                  '-source', atlas,
                  '-master', anatstemBc + '.nii.gz',
                  '-ainterp', 'NN',
                  '-prefix', atlas_Na1])

#%%

    #logging
    with open(logfile, 'a') as myfile:
        myfile.write("%s\n" % sink7)
        myfile.write("%s\n" % sink8)

#%%

    fixobliquity.fix_obliquity(atlas_Na1, anatstemBc + '.nii.gz', tmpdir, logfile)
