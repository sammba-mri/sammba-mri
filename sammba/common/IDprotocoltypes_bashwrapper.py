# -*- coding: utf-8 -*-

#%%

import sys
import IDprotocoltypes

#%%

target_dir = sys.argv[1]

#%%

protocol_dict = {
    'anat':['T1_FLASH_3D', 'MSME_MIRCen_allbrain', '__MSME_200um__'],
    'perfFAIREPI':['Perfusion_FAIR_EPI'],
    'rs':['T2star_FID_EPI_sat', 'GE_EPI_sat', 'SE_EPI_sat'],
    'DTIDKI':['DKI_EPI'],
    'CESTz':['zspectrum'],
    'CESTw':['WASSR']}

IDprotocoltypes.ID_protocol_types(protocol_dict, target_dir)
