# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 13:42:44 2017

@author: nadka
"""

import os
import subprocess
import glob

#get a proper pythoner to check this, am sure it's a disaster waiting to happen
def spco(*popenargs, **kwargs):
    x = subprocess.check_output(stderr=subprocess.STDOUT, *popenargs, **kwargs)
    return x

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

