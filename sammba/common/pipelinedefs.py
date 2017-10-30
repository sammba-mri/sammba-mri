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
import glob

#%%

#get a proper pythoner to check this, am sure it's a disaster waiting to happen
def spco(*popenargs, **kwargs):
    x = subprocess.check_output(stderr=subprocess.STDOUT, *popenargs, **kwargs)
    return x

#%%

def fix_obliquity(badfile, goodfile, tmpdir, logfile):

    #https://afni.nimh.nih.gov/afni/community/board/read.php?1,149385,149385#msg-149385
    #https://afni.nimh.nih.gov/afni/community/board/read.php?1,74225,74225#msg-74225
    
#%%

    os.environ["AFNI_DECONFLICT"] = "OVERWRITE"

#%%

    tmpdir = os.path.join(tmpdir, 'tmp_fixobliquity_'
                          + str(random.randrange(100000, 1000000)))
    os.mkdir(tmpdir)    
    goodobliquename = os.path.join(tmpdir, os.path.basename(goodfile)[:-7])
    badobliquename = os.path.join(tmpdir, os.path.basename(badfile)[:-7])

#%%

    sink1 = spco(['3dcopy', goodfile, goodobliquename + '+orig'])
    sink2 = spco(['3dcopy', badfile, badobliquename + '+orig'])

#%%

    sink3 = spco(['3drefit',
                  '-atrcopy',
                  goodobliquename + '+orig',
                  'IJK_TO_DICOM_REAL',
                  badobliquename + '+orig'])
    sink4 = spco(['3dcopy', badobliquename + '+orig', badfile])

#%%

    #logging
    with open(logfile, 'a') as myfile:
        myfile.write("%s\n" % sink1)
        myfile.write("%s\n" % sink2)
        myfile.write("%s\n" % sink3)
        myfile.write("%s\n" % sink4)

    shutil.rmtree(tmpdir)

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
    sink2 = spco(['3dClipLevel', anatstemBc + '.nii.gz'])
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

    #non-linear registration of affine pre-registered whole head image to
    #template. don't initiate straight from the original with an iniwarp due to
    #weird errors (like it creating an Allin it then can't find)
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

    fix_obliquity(atlas_Na1, anatstemBc + '.nii.gz', tmpdir, logfile)

#%%

def rs_to_template(rs_nx, anat_nx, template, tmpdir, scale, registerfunctional,
                   rminterfiles, brainvol, logfile):

#%%

    rs_nx = rs_nx[:-7]
    sink1 = spco(['3dinfo', '-nv', rs_nx + '.nii.gz'])

    if int(sink1.splitlines()[0]) >= 15:

#%%

        rs_nx_Ts = rs_nx + '_Ts'
        sink2 = spco(['3dTshift',
                      '-prefix', rs_nx_Ts + '.nii.gz',
                      '-tpattern', 'altplus',
                      rs_nx + '.nii.gz'])
        #remarkably, obliquity is preserved so no need to use fix_obliquity or
        #fsl5.0-fslcpgeom
    
    #%%
    
        sink3 = spco(['3dClipLevel',
                      rs_nx_Ts + '.nii.gz'])
        rs_nx_TsTm = rs_nx_Ts + 'Tm'
        sink4 = spco(['3dcalc',
                      '-a', rs_nx_Ts + '.nii.gz',
                      '-expr', 'a*ispositive(a-' + sink3.splitlines()[-1] +')',
                      '-prefix', rs_nx_TsTm + '.nii.gz'])
    
    #%%
    
        rs_nx_TsTmVr = rs_nx_TsTm + 'Vr'
        sink5 = spco(['3dvolreg',
                      '-prefix', rs_nx_TsTmVr + '.nii.gz',
                      '-dfile', rs_nx_TsTmVr + '.dfile.1D',
                      '-1Dfile', rs_nx_TsTmVr + '.1Dfile.1D',
                      '-1Dmatrix_save', rs_nx_TsTmVr + '.aff12.1D',
                      rs_nx_TsTm + '.nii.gz'])
    
    #%%
    
        rs_nx_TsAv = rs_nx_Ts + 'Av'
        sink6 = spco(['3dAllineate',
                      '-input', rs_nx_Ts + '.nii.gz',
                      '-master', rs_nx_Ts + '.nii.gz',
                      '-prefix', rs_nx_TsAv + '.nii.gz',
                      '-1Dmatrix_apply', rs_nx_TsTmVr + '.aff12.1D'])
    
    #%%
    
        #3dAllineate removes the obliquity. this is not a good way to readd it
        #as it removes motion correction info in the header if it were an AFNI
        #file; as it happens it's NIfTI which does not store that so irrelevant!
        sink7 = spco(['fsl5.0-fslcpgeom',
                      rs_nx_TsTmVr + '.nii.gz',
                      rs_nx_TsAv + '.nii.gz'])
        
        #create a (hopefully) nice mean image for use in the registration
        rs_nx_TsAvAv = rs_nx_TsAv + 'Av'
        sink8 = spco(['3dTstat',
                      '-mean',
                      '-prefix', rs_nx_TsAvAv + '.nii.gz',
                      rs_nx_TsAv + '.nii.gz'])

        #N4 fails for some reason. Not tried 3dUnifize yet
        rs_nx_TsAvAvN3 = rs_nx_TsAvAv + 'N3'
        sink9 = spco(['N3BiasFieldCorrection',
                      '3',
                      rs_nx_TsAvAv + '.nii.gz',
                      rs_nx_TsAvAvN3 + '.nii.gz'])

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
        br = os.path.basename(rs_nx)
        imfilearray = rs_nx + '_imfilearray.txt'
        with open(imfilearray, 'w') as myfile:
            myfile.write(br +'_TsAvAvN3' + ' ' + br + '_TsAv')

        startdir = os.path.abspath(os.path.curdir)
        sink10 = spco(['perslice_registration_subpipeline.bash',
                       os.path.dirname(rs_nx), 'anat_n' + str(anat_nx) + '_Un',
                       imfilearray, 'atlas_Na1', 'no', template, rminterfiles,
                       tmpdir, str(scale), registerfunctional, str(brainvol),
                       startdir])

#%%

        with open(logfile, 'a') as myfile:
            myfile.write("%s\n" % sink10)

#%%

        if rminterfiles == "yes":
            os.remove(rs_nx_Ts + '.nii.gz')
            os.remove(rs_nx_TsTm + '.nii.gz')
            os.remove(rs_nx_TsTmVr + '.nii.gz')
            os.remove(rs_nx_TsAv + '.nii.gz')

#%%

def rs_to_template_looper(NIfTIdir, anat_nx, template, tmpdir, scale,
                          registerfunctional, rminterfiles, brainvol, logfile):

    all_rs_nx = glob.glob(os.path.join(NIfTIdir, 'rs_n[0-9].nii.gz'))

    for rs_nx in all_rs_nx:
        rs_to_template(rs_nx, anat_nx, template, tmpdir, scale,
                       registerfunctional, rminterfiles, brainvol, logfile)

#%%

