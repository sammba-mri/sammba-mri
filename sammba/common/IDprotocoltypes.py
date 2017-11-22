# -*- coding: utf-8 -*-

#%%

import os

#%%

def ID_protocol_types(protocol_dict, target_dir):
    """
    to fill
    """

    target_dir = os.path.abspath(target_dir)
    filelist = []
    changelog = []

    #%%

    for root, dirs, files in os.walk(target_dir):
        for filename in files:
            if filename.endswith('.nii.gz'):
                filelist.append(filename)

    #%%

    for key, value in protocol_dict.items():
        minifilelist = []
        for file in filelist:
            if any(s in file for s in value):
                minifilelist.append(file)
        if len(minifilelist) > 0:
            minifilelist.sort()
            for n in range(0, len(minifilelist)):
                nosuffix = minifilelist[n].rstrip('.nii.gz')
                oldbasepath = os.path.join(target_dir, nosuffix)
                newbasepath = os.path.join(target_dir, key) + '_n' + str(n)
                os.rename(os.path.join(oldbasepath + '.nii.gz'),
                          os.path.join(newbasepath + '.nii.gz'))
                os.rename(os.path.join(oldbasepath + '_ptbl.txt'),
                          os.path.join(newbasepath + '_ptbl.txt'))
                changelog.append(nosuffix + ' to ' + key + '_n' + str(n))
        else:
            changelog.append('no files converted to ' + key)

    #%%

    f = open (os.path.join(target_dir, 'namechanges.txt'), 'w')
    for line in changelog:
        f.write("%s\n" % line)
    f.close()
