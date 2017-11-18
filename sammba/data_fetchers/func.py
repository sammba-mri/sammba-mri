import os
import numpy as np
import warnings
import json
from ..preprocessors.utils import _reset_affines
from sklearn.datasets.base import Bunch
from nilearn.datasets.utils import _fetch_files, _fetch_file, _get_dataset_dir
from .utils import _get_dataset_descr


def fetch_zurich_test_retest(subjects=range(15), sessions=[1], data_dir=None,
                             url=None, resume=True, verbose=1,
                             correct_headers=False):
    """Download and loads the ETH-Zurich test-retest dataset.

    Parameters
    ----------
    subjects : sequence of int or None, optional
        ids of subjects to load, default to loading all subjects.

    sessions : iterable of int, optional
        The sessions to load. Load only the first session by default.

    data_dir : string, optional
        Path of the data directory. Used to force data storage in a specified
        location. Default: None

    resume : bool, optional (default True)
        If true, try resuming download if possible.

    verbose : int, optional (default 0)
        Defines the level of verbosity of the output.

    Returns
    -------
    data : sklearn.datasets.base.Bunch
        Dictionary-like object, the interest attributes are :

        - 'func': string list. Paths to functional images.
        - 'anat': string list. Paths to anatomic images.
        - 'session': numpy array. List of ids corresponding to images sessions.

    Notes
    ------
    This dataset is composed of 2 sessions of 15 male mice.
    For each mice, 2 resting-state scans of continuous EPI
    functional volumes were collected, both with their anatomical scan.
    Session 2  was collected 15-20 days after Session 1.

    References
    ----------
    :Download:
        https://central.xnat.org

    :Reference:
        `Mapping the Mouse Brain with Rs-fMRI: An Optimized Pipeline for
        Functional Network Identification
        <http://dx.doi.org/10.1016/j.neuroimage.2015.07.090>`_
        NeuroImage 123 (2015): 11-21.
        V. Zerbi, J. Grandjean, M. Rudin, and N. Wenderoth.

    """
    if url is None:
        url = 'https://central.xnat.org'

    dataset_name = 'zurich_retest'
    data_dir = _get_dataset_dir(dataset_name, data_dir=data_dir,
                                verbose=verbose)

    # First, fetch the file that references all individual URLs
    json_file = _fetch_file(url + '/data/experiments.html', data_dir,
                            verbose=verbose)

    # Return the json file contents as a dictionary
    with open(json_file) as json_data:
        rows = json.load(json_data).values()[0]['Result']

    names = [name for name in rows[0].keys()]
    projects = {}
    for name in names:
        projects[name] = np.array([row[name] for row in rows])

    # Collect directories for all mice in the test-restest dataset
    subject_ids = ['1366', '1367', '1368', '1369', '1371', '1378', '1380',
                   '1402', '1403', '1404', '1405', '1406', '1407', '1411',
                   '1412']
    baseline_subject_ids = [subject + '_baseline' for subject in subject_ids]
    post_subject_ids = [subject + '_post' for subject in subject_ids]
    baseline_uris = projects['URI'][
        np.in1d(projects['label'], baseline_subject_ids)]
    post_uris = projects['URI'][
        np.in1d(projects['label'], post_subject_ids)]

    # Generate the list of urls by session
    func_file = 'rsfMRI.nii.gz'
    anat_file = '3DRARE.nii.gz'
    func_path = 'scans/rsfMRI/resources/NIFTI/files'
    anat_path = 'scans/anatomical/resources/NIFTI/files'
    func_urls = [
        [os.path.join(url + b, func_path, func_file) for b in baseline_uris],
        [os.path.join(url + p, func_path, func_file) for p in post_uris]
        ]
    anat_urls = [
        [os.path.join(url + b, anat_path, anat_file) for b in baseline_uris],
        [os.path.join(url + p, anat_path, anat_file) for p in post_uris]
        ]

    # Generate the list of target files by session
    func_files = [
        [os.path.join('baseline', sub, func_file) for sub in subject_ids],
        [os.path.join('post', sub, func_file) for sub in subject_ids]
        ]
    anat_files = [
        [os.path.join('baseline', sub, anat_file) for sub in subject_ids],
        [os.path.join('post', sub, anat_file) for sub in subject_ids]
        ]

    # Check arguments
    max_subjects = len(subject_ids)
    if max(subjects) > max_subjects:
        warnings.warn('Warning: there are only {0} subjects'.format(
            max_subjects))
        subjects = range(max_subjects)
    unique_subjects, indices = np.unique(subjects, return_index=True)
    if len(unique_subjects) < len(subjects):
        warnings.warn('Warning: Duplicate subjects, removing them.')
        subjects = unique_subjects[np.argsort(indices)]

    n_subjects = len(subjects)
    target_anat = []
    target_func = []
    source_anat = []
    source_func = []
    session = []
    for i in sessions:
        if not (i in [1, 2]):
            raise ValueError('Zurich dataset session id must be in [1, 2]')
        source_anat += [anat_urls[i - 1][subject] for subject in subjects]
        source_func += [func_urls[i - 1][subject] for subject in subjects]
        target_anat += [anat_files[i - 1][subject] for subject in subjects]
        target_func += [func_files[i - 1][subject] for subject in subjects]
        session += [i] * n_subjects

    # Call fetch_files once per subject.
    func = []
    anat = []
    for anat_u, anat_f, func_u, func_f in zip(source_anat, target_anat,
                                              source_func, target_func):
        a, f = _fetch_files(
            data_dir,
            [(anat_f, anat_u, {'move': anat_f}),
             (func_f, func_u, {'move': func_f})],  verbose=verbose)
        func.append(f)
        anat.append(a)

    fdescr = _get_dataset_descr(dataset_name)

    # This data has wrong sforms and qforms in the headers, so we correct them.
    if correct_headers:
        corrected_anat = []
        for a in anat:
            corrected_a = os.path.join(os.path.dirname(a),
                                       '3DRARE_corrected.nii.gz')
            _reset_affines(a, corrected_a,
                           axes_to_permute=[(1, 2)],
                           axes_to_flip=[0],
                           verbose=0)
            corrected_anat.append(corrected_a)
        corrected_func = []
        for f in func:
            corrected_f = os.path.join(os.path.dirname(f),
                                       'rsfMRI_corrected.nii.gz')
            _reset_affines(f, corrected_f,
                           center_mass=(0, 0, 0),
                           xyzscale=.1,
                           axes_to_permute=[(1, 2)],
                           axes_to_flip=[0],
                           verbose=0)
            corrected_func.append(corrected_f)
        anat = corrected_anat
        func = corrected_func

    return Bunch(anat=anat, func=func, session=session, description=fdescr)


def fetch_zurich_anesthesiant(subjects=range(30), url=None,
                              data_dir=None, resume=True, verbose=1):
    """Download and loads the ETH-Zurich anesthesiant dataset.

    Parameters
    ----------
    subjects : sequence of int or None, optional
        ids of subjects to load, default to loading all subjects.

    data_dir: string, optional
        Path of the data directory. Used to force data storage in a specified
        location. Default: None

    resume: bool, optional (default True)
        If true, try resuming download if possible.

    verbose: int, optional (default 0)
        Defines the level of verbosity of the output.

    Returns
    -------
    data : sklearn.datasets.base.Bunch
        Dictionary-like object, the interest attributes are:
        - 'func': string list. Paths to functional images.
        - 'anesthesiant': string list. Information on used anesthesiant.

    Notes
    ------
    This dataset is composed of 30 male mice with different anesthesia
    protocols.

    References
    ----------
    :Download:
        https://central.xnat.org

    :Reference:
        `Optimization of anesthesia protocol for resting-state fMRI in mice
        based on differential effects of anesthetics on functional connectivity
        patterns.
        <http://dx.doi.org/10.1016/j.neuroimage.2014.08.043>`_
        NeuroImage 102 (2014): 838-847.
        J. Grandjean and A. Schroeter and I. Batata and M. Rudin.
    """
    if url is None:
        url = 'https://central.xnat.org'

    dataset_name = 'zurich_anest'
    data_dir = _get_dataset_dir(dataset_name, data_dir=data_dir,
                                verbose=verbose)

    # First, fetch the file that references all individual URLs
    json_file = _fetch_file(os.path.join(url, 'data', 'experiments.html'),
                            data_dir, verbose=verbose)

    # Return the json file contents as a dictionary
    with open(json_file) as json_data:
        rows = json.load(json_data).values()[0]['Result']

    names = [name for name in rows[0].keys()]
    projects = {}
    for name in names:
        projects[name] = np.array([row[name] for row in rows])

    # Collect directories for all mice in the anesthesiant dataset
    iso_ids = ['iso2273', 'iso2274', 'iso2238', 'iso2239', 'iso2250',
               'iso2270', 'iso3294', 'iso3296']
    med_ids = ['med2259', 'med2241', 'med2247', 'med2251', 'med2256',
               'med2257']
    mi_ids = ['mi272871', 'mi273299', 'mi273457', 'mi273458', 'mi273459',
              'mi273460', 'mi273461', 'mi273300']
    med_half_dose_ids = ['medHalfDose', 'medHalfDose1', 'medHalfDose2',
                         'medHalfDose3']
    iso1_c3_ids = ['iso1c3perc', 'iso1c3perc']
    iso1_c5_ids = ['iso1c5perc', 'iso2870_1c5perc']

    subjects_ids = iso_ids + med_ids + mi_ids + med_half_dose_ids + \
        iso1_c3_ids + iso1_c5_ids
    subjects_labels = ['Iso1'] * len(iso_ids) + ['Med'] * len(med_ids) + \
                      ['Med-Iso'] * len(mi_ids) + \
                      ['Med-half'] * len(med_half_dose_ids) + \
                      ['Iso1pt3'] * len(iso1_c3_ids) + \
                      ['Iso1pt5'] * len(iso1_c5_ids)

    max_subjects = len(subjects_ids)

    # Check arguments
    max_subjects = len(subjects_ids)
    if subjects is None:
        subjects = range(max_subjects)
    elif max(subjects) > max_subjects:
        warnings.warn('Warning: there are only {0} subjects'.format(
            max_subjects))
        subjects = range(max_subjects)
    unique_subjects, indices = np.unique(subjects, return_index=True)
    if len(unique_subjects) < len(subjects):
        warnings.warn('Warning: Duplicate subjects, removing them.')
        subjects = unique_subjects[np.argsort(indices)]

    subjects_ids = [subjects_ids[subject] for subject in subjects]
    subjects_labels = [subjects_labels[subject] for subject in subjects]

    mice_uris = projects['URI'][np.in1d(projects['label'], subjects_ids)]

    # Generate the list of urls by session
    img_file = 'rsfMRI.img'
    hdr_file = 'rsfMRI.hdr'
    func_path = 'scans/rs_fMRI/resources/NULL/files'
    img_urls = [os.path.join(url + b, func_path, img_file) for b in mice_uris]
    hdr_urls = [os.path.join(url + b, func_path, hdr_file) for b in mice_uris]

    # Generate the list of target files by session
    target_img = [os.path.join(label, sub, img_file)
                  for sub, label in zip(subjects_ids, subjects_labels)]
    target_hdr = [os.path.join(label, sub, hdr_file)
                  for sub, label in zip(subjects_ids, subjects_labels)]

    # Call fetch_files once per subject.
    img = []
    for img_u, hdr_u, img_f, hdr_f in zip(img_urls, hdr_urls, target_img,
                                          target_hdr):
        f, _ = _fetch_files(
            data_dir, [(img_f, img_u, {'move': img_f}),
                       (hdr_f, hdr_u, {'move': hdr_f})],  verbose=verbose)
        img.append(f)

    fdescr = _get_dataset_descr(dataset_name)

    return Bunch(func=img, anesthesiant=subjects_labels, description=fdescr)
