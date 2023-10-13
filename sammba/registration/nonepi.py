import warnings
from sklearn.utils import Bunch
from numpy import VisibleDeprecationWarning
from .base import (_reorient, _rigid_body_register_and_reorient)


def _coregister_nonepi(unifized_anat_file,
                       unbiased_nonepi_file, write_dir,
                       modality_name='nonepi',
                       anat_brain_file=None,
                       nonepi_brain_file=None,
                       prior_rigid_body_registration=None,
                       reorient_only=False, caching=False,
                       verbose=True, **environ_kwargs):
    """
    Coregistration of the subject's non EPI modality and anatomical images.
    The non EPI modality volume is aligned to the anatomical, with a possible
    rigid body registration.

    Parameters
    ----------
    prior_rigid_body_registration : bool, optional
        If True, a rigid-body registration of the anat to the EPI is
        performed prior to the warp. Useful if the images headers have
        missing/wrong information.
        NOTE: prior_rigid_body_registration is deprecated from 0.1 and will be
        removed in next release. Use `reorient_only` instead.
    reorient_only :  bool, optional
        If True, the rigid-body registration of the anat to the EPI is not
        performed and only reorientation is done.
    voxel_size_x : float, optional
        Resampling resolution for the x-axis, in mm.
    voxel_size_y : float, optional
        Resampling resolution for the y-axis, in mm.
    caching : bool, optional
        Wether or not to use caching.
    verbose : bool, optional
        If True, all steps are verbose. Note that caching implies some
        verbosity in any case.
    environ_kwargs : extra arguments keywords
        Extra arguments keywords, passed to interfaces environ variable.

    Returns
    -------
    data : sklearn.utils.Bunch
    Dictionary-like object, the interest attributes are :
        - `coreg_anat_` : str
                          Path to paths to the coregistered anat image.
        - `coreg_transform_` : str
                               Path to the transform from anat to modality.
    """
    if prior_rigid_body_registration is not None:
        warn_str = ("The parameter 'prior_rigid_body_registration' is "
                    "deprecated and will be removed in sammba-mri next "
                    "release. Use parameter 'reorient_only' instead.")
        warnings.warn(warn_str, VisibleDeprecationWarning, stacklevel=2)
        reorient_only = not(prior_rigid_body_registration)

    environ = {'AFNI_DECONFLICT': 'OVERWRITE'}
    for (key, value) in environ_kwargs.items():
        environ[key] = value

    if verbose:
        terminal_output = 'allatonce'
    else:
        terminal_output = 'none'

    if reorient_only:
        registered_anat_oblique_file, transform_file = \
        _reorient(unifized_anat_file, unbiased_nonepi_file, write_dir,
                  modality_name=modality_name,
                  terminal_output=terminal_output,
                  caching=caching, verbose=verbose, environ=environ)
    else:
        if anat_brain_file is None:
            raise ValueError("'anat brain mask file' is needed for "
                             "rigid-body registration")
        if nonepi_brain_file is None:
            raise ValueError("'non EPI modality brain mask file' is needed "
                             "for rigid-body registration")
        registered_anat_oblique_file, transform_file = \
            _rigid_body_register_and_reorient(unifized_anat_file,
                                              unbiased_nonepi_file, write_dir,
                                              anat_brain_file,
                                              nonepi_brain_file,
                                              modality_name=modality_name,
                                              terminal_output=terminal_output,
                                              caching=caching,
                                              verbose=verbose,
                                              environ=environ)

    return Bunch(coreg_anat_=registered_anat_oblique_file,
                 coreg_transform_=transform_file)
