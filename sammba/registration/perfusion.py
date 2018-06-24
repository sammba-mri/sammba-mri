import os
from sklearn.datasets.base import Bunch
from sammba.externals.nipype.caching import Memory
from sammba.externals.nipype.interfaces import afni, ants
from sammba.externals.nipype.utils.filemanip import fname_presuffix
from .base import (_rigid_body_register, _warp, _per_slice_qwarp)


def coregister(unifized_anat_file,
               unbiased_m0_file, write_dir, anat_brain_file=None,
               m0_brain_file=None,
               prior_rigid_body_registration=False,
               apply_to_file=None,
               voxel_size_x=.1, voxel_size_y=.1, caching=False,
               verbose=True, **environ_kwargs):
    """
    Coregistration of the subject's M0 and anatomical images.
    The M0 volume is aligned to the anatomical, first with a
    rigid body registration and then on a per-slice basis (only a fine
    correction, this is mostly for correction of EPI distortion).

    Parameters
    ----------
    use_rats_tool : bool, optional
        If True, brain mask is computed using RATS Mathematical Morphology.
        Otherwise, a histogram-based brain segmentation is used.
    prior_rigid_body_registration : bool, optional
        If True, a rigid-body registration of the anat to the func is
        performed prior to the warp. Useful if the images headers have
        missing/wrong information.
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
    data : sklearn.datasets.base.Bunch
    Dictionary-like object, the interest attributes are :
        - `coreg_perf_` : str
                          Path to paths to the coregistered perfusion
                          image.
        - `coreg_anat_` : str
                          Path to paths to the coregistered functional
                          image.
        - `coreg_transform_` : str
                               Path to the transform from anat to func.
    Notes
    -----
    If `use_rats_tool` is turned on, RATS tool is used for brain extraction
    and has to be cited. For more information, see
    `RATS <http://www.iibi.uiowa.edu/content/rats-overview/>`_
    """
    environ = {'AFNI_DECONFLICT': 'OVERWRITE'}
    for (key, value) in environ_kwargs.items():
        environ[key] = value

    if verbose:
        terminal_output = 'allatonce'
    else:
        terminal_output = 'none'

    if caching:
        memory = Memory(write_dir)
        catmatvec = memory.cache(afni.CatMatvec)
        overwrite = False
    else:
        catmatvec = afni.CatMatvec().run
        overwrite = True

    output_files = []

    #############################################
    # Rigid-body registration anat -> mean func #
    #############################################
    if prior_rigid_body_registration:
        if anat_brain_file is None:
            raise ValueError("'anat_brain_mask_file' is needed for prior "
                             "rigid-body registration")
        if m0_brain_file is None:
            raise ValueError("'m0_brain_mask_file' is needed for prior "
                             "rigid-body registration")
        allineated_anat_file, rigid_transform_file = \
            _rigid_body_register(unifized_anat_file,
                                 unbiased_m0_file,
                                 write_dir, anat_brain_file,
                                 m0_brain_file, 
                                 caching=caching,
                                 terminal_output=terminal_output,
                                 environ=environ)
        output_files.extend([rigid_transform_file,
                             allineated_anat_file])
    else:
        allineated_anat_file = unifized_anat_file

    ############################################
    # Nonlinear registration anat -> mean func #
    ############################################
    registered_anat_oblique_file, mat_file =\
        _warp(allineated_anat_file, unbiased_m0_file, write_dir,
              caching=caching, verbose=verbose,
              terminal_output=terminal_output, overwrite=overwrite,
              environ=environ)

    # Concatenate all the anat to func tranforms
    output_files.append(mat_file)
    transform_file = fname_presuffix(registered_anat_oblique_file,
                                     suffix='_anat_to_func.aff12.1D',
                                     use_ext=False)
    _ = catmatvec(in_file=[(mat_file, 'ONELINE')],
                  oneline=True,
                  out_file=transform_file,
                  environ=environ)

    ##################################################
    # Per-slice non-linear registration func -> anat #
    ##################################################
    warped_m0_file, warp_files, warped_apply_to_file =\
        _per_slice_qwarp(unbiased_m0_file,
                         registered_anat_oblique_file,
                         voxel_size_x, voxel_size_y,
                         apply_to_file=apply_to_file,
                         verbose=verbose,
                         write_dir=write_dir,
                         caching=caching, terminal_output=terminal_output,
                         environ=environ)

    # Remove the intermediate outputs
    if not caching:
        for out_file in output_files:
            os.remove(out_file)

    return Bunch(coreg_apply_to_=warped_apply_to_file,
                  coreg_m0_=warped_m0_file,
                  coreg_anat_=registered_anat_oblique_file,
                  coreg_transform_=transform_file,
                  coreg_warps_=warp_files)
