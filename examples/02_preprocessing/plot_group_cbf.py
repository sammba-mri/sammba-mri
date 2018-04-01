"""Compute CBF example
"""
import glob
import os
dcm_dirs = glob.glob('/home/salma/nilearn_data/mlrs_perfusion/MLRS4-PERF-3/20170712/*')

# Convert dicom to niftis
from sammba.io_conversions import recursive_dcm_to_nii
perf_filenames = []
for dcm_dir in dcm_dirs:
    save_dir = os.path.join(dcm_dir, 'niftis')
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)

    converted_files = recursive_dcm_to_nii('/usr/bin/dcmdump', dcm_dir, save_dir)
    perf_filenames.extend(fnmatch.filter(converted_files, '*perfFAIR*.nii.gz'))

from sammba.modality_processors import perfusion_fair
processed_files = []
for perf_filename in perf_filenames:
    perfusion_fair.perf_fair_niiptbl_proc(perf_filename, 2800)
    processed_files.append(fname_presuffix(perf_filename, suffix='_proc'))

from nilearn import image
for proc_file in processed_files:
    m0_imgs = image.index_img(proc_file, [2, 8])
    m0_img = image.mean_img(m0_imgs)
    m0_img.to_filename(proc_file.replace('_proc', '_m0'))

from sammba.registration import perfusion

perf_session = perfusion.PerfSession(anat=anat, perf=perf, m0=m0_file)
perf_session.coregister(verbose=3, caching=True)
perf_session.coreg_perf_
nibabel.load(perf_session.coreg_perf_).shape
import nibabel
nibabel.load(perf_session.coreg_perf_).shape
nibabel.load(perf_session.coreg_m0_).shape
perf = '/home/salma/nilearn_data/mlrs_perfusion/MLRS4-PERF-3/20170712/20170712_152027_MINDt_MLRS_28944_1_2/niftis/MINDt-MLRS-28944_20170712_154434_06_Perfusion_FAIR_EPI_fixedSIAP_bf9_proc.nii.gz'
perf_session = perfusion.PerfSession(anat=anat, perf=perf, m0=m0_file)
perf_session.coregister(verbose=3, caching=True)
nibabel.load(perf_session.coreg_perf_).shape
cbf = image.index_img(perf_session.coreg_perf_, 14)
cbf = image.index_img(perf_session.coreg_perf_, 12)
