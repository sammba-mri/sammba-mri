# -*- coding: utf-8 -*-

import os
import glob
import subprocess
import shutil

#%%

def func_to_template(sources_directory, master_3dNwarpApply, target_directory,
                     transform_list, byfolder):

#%%

    input_paths = [sources_directory, master_3dNwarpApply, target_directory]
    input_paths = [os.path.abspath(input_path) for input_path in input_paths]
    sources_directory, master_3dNwarpApply, target_directory = input_paths
    source_directories = [os.path.join(sources_directory, x) for x in
                                           next(os.walk(sources_directory))[1]]

#%%

    for source_directory in source_directories:

#%%

        exist = glob.glob(os.path.join(source_directory, 'rs_n*_TsTmVr.1Dfile.1D'))
        exist = [x[(len(source_directory) + 5):-17] for x in exist]
        target_individual = os.path.join(target_directory,
                                         os.path.basename(source_directory))

#%%

        if byfolder == 'yes':
            os.mkdir(target_individual)            

#%%

        for n in exist:

#%%
            #sp and tp mean source and target prefixes
            p = 'rs_n' + n
            new_transform_list = [transform.replace('rs_nx', p) for transform in transform_list]
            new_transform_list = [os.path.join(source_directory, transform) for transform in new_transform_list]            
            sp = os.path.join(source_directory, p)
            if byfolder == 'no':
                tp = target_individual + '__' + p
            elif byfolder == 'yes':
                tp = os.path.join(target_individual, p)
            func = sp + '_TsAv_NaMe.nii.gz'
            prefix_3dNwarpApply = tp + '_TsAv_NaMeNa.nii.gz'

#%%

            subprocess.call(['3dNwarpApply', '-nwarp'] + new_transform_list +
                            ['-source', func, '-master', master_3dNwarpApply,
                             '-prefix', prefix_3dNwarpApply])
            shutil.copyfile(sp + '_TsTmVr.1Dfile.1D', tp + '_TsTmVr.1Dfile.1D')
