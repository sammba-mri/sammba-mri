import numpy as np
from sklearn.datasets.base import Bunch
from nilearn.datasets.utils import _get_dataset_dir, _fetch_files
from .utils import _get_dataset_descr


def fetch_atlas_dorr_2008(data_dir=None, url=None, resume=True, verbose=1,
                          down_sampled=True):
    """Download and load Dorr et al. atlas and average (dated 2008)

    Parameters
    ----------
    data_dir: string, optional
        Path of the data directory. Use to forec data storage in a non-
        standard location. Default: None (meaning: default)
    url: string, optional
        Download URL of the dataset. Overwrite the default URL.

    resume: bool
        whether to resumed download of a partly-downloaded file.

    verbose: int
        verbosity level (0 means no message).

    Returns
    -------
    data: sklearn.datasets.base.Bunch
        dictionary-like object, contains:

        - 't2': str, path to nifti file containing the T2 weighted average.

        - 'maps': str, path to nifti file containing regions definition.

        - 'labels': str list containing the names of the regions.

        - 'indices': int list containing the value of each region.

        - 'description': description about the atlas and the template.

    References
    ----------

    A.E. Dorr, J.P. Lerch, S. Spring, N. Kabani and R.M. Henkelman. "High
    resolution three dimensional brain atlas using an average magnetic
    resonance image of 40 adult C57Bl/6j mice", NeuroImage 42(1):60-69, 2008.

    See http://www.mouseimaging.ca/research/mouse_atlas.html for more
    information on this parcellation.

    Licence: Unknown
    """
    if url is None:
        url = ['http://www.mouseimaging.ca/mnc/C57Bl6j_mouse_atlas/',
               'http://www.mouseimaging.ca/mnc/C57Bl6j_mouse_atlas/',
               'http://www.mouseimaging.ca/research/C57Bl6j_mouse_atlas/']
    files = ['male-female-mouse-atlas.mnc', 'c57_fixed_labels_resized.mnc',
             'c57_brain_atlas_labels.csv']

    files = [(f, u + f, {}) for f, u in zip(files, url)]

    dataset_name = 'dorr_2008'
    data_dir = _get_dataset_dir(dataset_name, data_dir=data_dir,
                                verbose=verbose)
    files_ = _fetch_files(data_dir, files, resume=resume,
                          verbose=verbose)

    fdescr = _get_dataset_descr(dataset_name)
    csv_data = np.recfromcsv(
        files_[2], skip_header=True,
        names=('roi_id', 'roi_label', 'right_index', 'left_index'))

    #TODO try dictionary with their region id as key and name as value
    left_rois = []
    right_rois = []
    lateral_rois = []
    for (idx, label, right_index, left_index) in csv_data:
        if right_index == left_index:
            lateral_rois.append((label, right_index))
        else:
            left_rois.append(('L ' + label, left_index))
            right_rois.append(('R ' + label, right_index))

    rois = lateral_rois + right_rois + left_rois
    labels, indices = zip(*rois)
    t2 = files_[0]
    maps = files_[1]

    params = dict(t2=t2, maps=maps,
                  labels=np.array(labels)[np.argsort(indices)],
                  indices=np.sort(indices),
                  description=fdescr)

    return Bunch(**params)
