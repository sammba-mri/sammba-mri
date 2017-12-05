import json
import numpy as np
from sklearn.datasets.base import Bunch
from nilearn.datasets.utils import _get_dataset_dir, _fetch_files
from .utils import _get_dataset_descr


def fetch_atlas_dorr_2008(image_format='nifti', data_dir=None, url=None,
                          resume=True, verbose=1, down_sampled=True):
    """Download and load Dorr et al. atlas and average (dated 2008)

    Parameters
    ----------
    image_format : one of {'nifti', 'minc'}, optional
        Format to download

    data_dir : str, optional
        Path of the data directory. Use to forec data storage in a non-
        standard location. Default: None (meaning: default)

    url: string, optional
        Download URL of the dataset. Overwrite the default URL.

    resume : bool
        whether to resumed download of a partly-downloaded file.

    verbose : int
        verbosity level (0 means no message).

    Returns
    -------
    data: sklearn.datasets.base.Bunch
        dictionary-like object, contains:

        - 't2' : str, path to nifti file containing the T2 weighted average.

        - 'maps' : str, path to nifti file containing regions definition.

        - 'names' : str list containing the names of the regions.

        - 'values' : int list containing the label value of each region.

        - 'description' : description about the atlas and the template.

    References
    ----------

    A.E. Dorr, J.P. Lerch, S. Spring, N. Kabani and R.M. Henkelman. "High
    resolution three dimensional brain atlas using an average magnetic
    resonance image of 40 adult C57Bl/6j mice", NeuroImage 42(1):60-69, 2008.

    See http://www.mouseimaging.ca/research/mouse_atlas.html for more
    information on this parcellation.

    Licence: Unknown
    """
    if image_format not in ['nifti', 'minc']:
        raise ValueError("Images format must be 'nifti' or 'minc', you "
                         "entered {0}".format(image_format))

    if url is None:
        if image_format == 'minc':
            url = ['http://www.mouseimaging.ca/mnc/C57Bl6j_mouse_atlas/',
                   'http://www.mouseimaging.ca/mnc/C57Bl6j_mouse_atlas/',
                   'http://www.mouseimaging.ca/research/C57Bl6j_mouse_atlas/']
        else:
            url = ['http://repo.mouseimaging.ca/repo/Dorr_2008_nifti/',
                   'http://repo.mouseimaging.ca/repo/Dorr_2008_nifti/',
                   'http://www.mouseimaging.ca/research/C57Bl6j_mouse_atlas/']

    if image_format == 'minc':
        files = ['male-female-mouse-atlas.mnc', 'c57_fixed_labels_resized.mnc',
                 'c57_brain_atlas_labels.csv']
    else:
        files = ['Dorr_2008_average.nii.gz', 'Dorr_2008_labels.nii.gz',
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
                  names=np.array(labels)[np.argsort(indices)],
                  values=np.sort(indices),
                  description=fdescr)

    return Bunch(**params)


def fetch_atlas_waxholm_rat_2014(data_dir=None, url=None, resume=True,
                                 verbose=1, downsample='2',
                                 symmetric_split=False):
    """Download and load Pape et al. rat atlas (dated 2014), downsampled
       by the Scalable Brain Atlas.

    Parameters
    ----------
    data_dir : str, optional
        Path of the data directory. Use to forec data storage in a non-
        standard location. Default: None (meaning: default)

    downsample : one of {'2', '3'}, optional
        Downsampling version, corresponding to resolutions 78 or 117 microns.

    url : string, optional
        Download URL of the dataset. Overwrite the default URL.

    resume : bool, optional
        Whether to resumed download of a partly-downloaded file.

    verbose : int, optional
        Verbosity level (0 means no message).

    Returns
    -------
    data: sklearn.datasets.base.Bunch
        dictionary-like object, contains:

        - 't2star': str, path to nifti file containing the averaged T2* images.

        - 'maps': str, path to nifti file containing regions definition.

        - 'labels': structured numpy.array containing the names of the regions
                    and their label values.

        - 'description': description about the atlas.

    References to cite
    ------------------
    The downsampled versions are provided by the Scalable Brain Atlas.
    The defining citations are

    Papp, Eszter A., Trygve B. Leergaard, Evan Calabrese, G. Allan Johnson,
    and Jan G. Bjaalie.
    `Waxholm Space atlas of the Sprague Dawley rat brain
     <http://dx.doi.org/10.1016/j.neuroimage.2014.04.001>`_
    NeuroImage 97 (2014): 374-86.

    Kjonigsen LJ, Lillehaug S, Bjaalie JG, Witter MP, Leergaard TB.
    `Waxholm Space atlas of the rat brain hippocampal region:
    Three-dimensional delineations based on magnetic resonance
    and diffusion tensor imaging.
    <http://dx.doi.org/10.1016/j.neuroimage.2014.12.080>`_
    NeuroImage 108 (2015):441-449

    Sergejeva M, Papp EA, Bakker R, Gaudnek MA, Okamura-Oho Y, Boline J,
    Bjaalie JG, Hess A. `Anatomical landmarks for registration of experimental
    image data to volumetric rodent brain atlasing templates.
    <http://dx.doi.org/10.1016/j.jneumeth.2014.11.005>`_
    Journal of Neuroscience Methods (2015) 240:161-169.

    See https://scalablebrainatlas.incf.org/rat/PLCJB14 for more
    information on this parcellation.

    Licence
    -------
    Creative Commons Attribution-NonCommercial-ShareAlike 4.0
    International
    """
    if downsample not in ['2', '3']:
        raise ValueError(
            "'downsample' must be '2' or '3', you provided".format(downsample))
    if url is None:
        base_url = 'https://scalablebrainatlas.incf.org/'
        downsampled_atlas = 'WHS_SD_rat_atlas_v1.01_' +\
                            'downsample{0}.nii.gz'.format(downsample)
        downsampled_t2star = 'WHS_SD_rat_T2star_v1.01_' +\
                             'downsample{0}.nii.gz'.format(downsample)
        url = [base_url + 'templates/PLCJB14/source/' + downsampled_t2star,
               base_url + 'templates/PLCJB14/source/' + downsampled_atlas,
               base_url + 'services/labelmapper.php?' +
               'template=PLCJB14&to=acr&format=json']

    files = [downsampled_t2star, downsampled_atlas, 'WHS_SD_rat_labels.json']

    opts = [{}, {}, {'move': 'WHS_SD_rat_labels.json'}]
    files = [(f, u, opt) for (f, u, opt) in zip(files, url, opts)]

    dataset_name = 'waxholm_rat_2014'
    data_dir = _get_dataset_dir(dataset_name, data_dir=data_dir,
                                verbose=verbose)
    files_ = _fetch_files(data_dir, files, resume=resume,
                          verbose=verbose)
    fdescr = _get_dataset_descr(dataset_name)

    # Return the json file contents as a dictionary
    with open(files_[2]) as json_data:
        json_rows = json.load(json_data)

    # Convert it to structured array
    names = ['value', 'name']
    formats = ['|S3', '|S59']
    dtype = dict(names=names, formats=formats)
    labels = np.array(list(json_rows.items()), dtype=dtype)

    # TODO: symmetric_split
    if symmetric_split:
        raise NotImplementedError('Not yet implemented')

    params = dict(t2star=files_[0], maps=files_[1], labels=labels,
                  description=fdescr)

    return Bunch(**params)
