import os
import warnings
import csv
import numpy as np
from sklearn.datasets.base import Bunch
from nilearn.datasets.utils import _fetch_files, _fetch_file, _get_dataset_dir
from .utils import _get_dataset_descr


def fetch_lemur_mircen_2019_t2(subjects=[0], data_dir=None, url=None,
                               resume=True, verbose=1):
    """Download and loads the mouse lemur template dataset.

    Parameters
    ----------
    subjects : sequence of int or None, optional
        ids of subjects to load, default to loading one subject.

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

        - 'anat': string list. Paths to T2-weighted images.
        - 'phenotypic': Participants genders, birth dates and MRI scan dates

    References
    ----------
    :Download:
        https://openneuro.org/datasets/ds001945/versions/1.0.0/download

    :Reference:
        `A 3D population-based brain atlas of the mouse lemur primate with
        examples of applications in aging studies and comparative anatomy.
        <http://doi:10.1016/j.neuroimage.2018.10.010>`_
        Neuroimage 185 (2019): 85-95. 
        N. A. Nadkarni, S. Bougacha, C. Garin, M. Dhenain, and J. L. Picq. 

    """
    if url is None:
        url = 'https://openneuro.org/crn/datasets/ds001945/snapshots/1.0.0/files'

    dataset_name = 'mircen2019_t2'
    data_dir = _get_dataset_dir(dataset_name, data_dir=data_dir,
                                verbose=verbose)

    # Check arguments
    max_subjects = 34
    if max(subjects) > max_subjects:
        warnings.warn('Warning: there are only {0} subjects'.format(
            max_subjects))
        subjects = range(max_subjects)

    subject_ids = np.array([
        '"sub-01"', '"sub-02"', '"sub-03"', '"sub-04"', '"sub-05"',
        '"sub-06"', '"sub-07"', '"sub-08"', '"sub-09"', '"sub-10"',
        '"sub-11"', '"sub-12"', '"sub-13"', '"sub-14"', '"sub-15"',
        '"sub-16"', '"sub-17"', '"sub-18"', '"sub-19"', '"sub-20"',
        '"sub-21"', '"sub-22"', '"sub-23"', '"sub-24"', '"sub-25"',
        '"sub-26"', '"sub-27"', '"sub-28"', '"sub-29"', '"sub-30"',
        '"sub-31"', '"sub-32"', '"sub-33"', '"sub-34"'])
    subject_ids = subject_ids[subjects]

    # Generate the list of urls
    json_urls = [os.path.join(url, '{0}:anat:{0}_T2w.json'.format(subject_id))
                 for subject_id in subject_ids]
    anat_urls = [os.path.join(url, '{0}:anat:{0}_T2w.nii.gz'.format(subject_id))
                 for subject_id in subject_ids]

    # Generate the list of target files
    anat_basenames = ['{0}_anat_{0}_T2w.nii.gz'.format(subject_id)
                      for subject_id in subject_ids]
    anat_files = [os.path.join(animal_dir, anat_basename)
        for (animal_dir, anat_basename) in zip(subject_ids, anat_basenames)]

    json_basenames = ['{0}_anat_{0}_T2w.json'.format(subject_id)
                      for subject_id in subject_ids]
    json_files = [os.path.join(animal_dir, json_basename)
        for (animal_dir, json_basename) in zip(subject_ids, json_basenames)]

    # Call fetch_files once per subject.
    anat = []
    json = []
    for anat_u, anat_f, json_u, json_f in zip(anat_urls, anat_files,
                                              json_urls, json_files):
        a, j = _fetch_files(
            data_dir,
            [(anat_f, anat_u, {'move': anat_f}),
             (json_f, json_u, {'move': json_f})],  verbose=verbose)
        json.append(j)
        anat.append(a)

    pheno_url = os.path.join(url, 'lemur_atlas_list_t2_bids.csv')
    pheno_file = _fetch_file(pheno_url, data_dir, verbose=verbose)
    phenotypic = np.recfromcsv(pheno_file, delimiter='\t')
    phenotypic = phenotypic[[np.where(phenotypic['animal_id'] == i)[0][0]
                            for i in subject_ids]]
    fdescr = _get_dataset_descr(dataset_name)

    return Bunch(anat=anat, pheno=phenotypic, description=fdescr)
