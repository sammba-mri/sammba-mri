import os
import numpy as np
import nibabel
from sammba.externals.nipype.caching import Memory
from sammba.externals.nipype.interfaces import afni, fsl
from sammba.externals.nipype.utils.filemanip import fname_presuffix
from sammba.interfaces import RatsMM
from .utils import fix_obliquity
from .fmri_session import FMRISession
from .struct import anats_to_template


def coregister_fmri_session(session_data, t_r, write_dir, brain_volume,
                            slice_timing=True,
                            prior_rigid_body_registration=False,
                            caching=False, voxel_size_x=.1, voxel_size_y=.1,
                            verbose=True, **environ_kwargs):
    """
    Coregistration of the subject's functional and anatomical images.
    The functional volume is aligned to the anatomical, first with a rigid body
    registration and then on a per-slice basis (only a fine correction, this is
    mostly for correction of EPI distortion).


    Parameters
    ----------
    session_data : sammba.registration.SessionData
        Single animal data, giving paths to its functional and anatomical
        image, as well as it identifier.

    t_r : float
        Repetition time for the EPI, in seconds.

    write_dir : str
        Directory to save the output and temporary images.

    brain_volume : int
        Volumes of the brain as passed to Rats_MM brain extraction tool.
        Typically 400 for mouse and 1800 for rat.

    prior_rigid_body_registration : bool, optional
        If True, a rigid-body registration of the anat to the func is performed
        prior to the warp. Useful if the images headers have missing/wrong
        information.

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
    the same sequence with each animal_data updated: the following attributes
    are added
        - `output_dir_` : str
                          Path to the output directory.
        - `coreg_func_` : str
                          Path to paths to the coregistered functional image.
        - `coreg_anat_` : str
                          Path to paths to the coregistered functional image.
        - `coreg_transform_` : str
                               Path to the transform from anat to func.
    """
    func_filename = session_data.func
    anat_filename = session_data.anat

    environ = {'AFNI_DECONFLICT': 'OVERWRITE'}
    for (key, value) in environ_kwargs.items():
        environ[key] = value

    if verbose:
        terminal_output = 'allatonce'
    else:
        terminal_output = 'none'

    if caching:
        memory = Memory(write_dir)
        tshift = memory.cache(afni.TShift)
        unifize = memory.cache(afni.Unifize)
        catmatvec = memory.cache(afni.CatMatvec)
        for step in [tshift, unifize]:
            step.interface().set_default_terminal_output(terminal_output)
    else:
        tshift = afni.TShift(terminal_output=terminal_output).run
        unifize = afni.Unifize(terminal_output=terminal_output).run
        catmatvec = afni.CatMatvec().run

    session_data._check_inputs()
    output_dir = os.path.join(os.path.abspath(write_dir),
                              session_data.animal_id)
    session_data._set_output_dir_(output_dir)
    current_dir = os.getcwd()
    os.chdir(output_dir)
    output_files = []

    #######################################
    # Correct functional for slice timing #
    #######################################
    if slice_timing:
        out_tshift = tshift(in_file=func_filename,
                            outputtype='NIFTI_GZ',
                            tpattern='altplus',
                            tr=str(t_r),
                            environ=environ)
        func_filename = out_tshift.outputs.out_file

    ################################################
    # Register functional volumes to the first one #
    ################################################
    allineated_filename, mean_aligned_filename, _ = \
        _realign(func_filename, write_dir, caching=caching,
                 terminal_output=terminal_output, environ=environ)

    ###########################################
    # Corret anat and func for intensity bias #
    ###########################################
    # Correct the functional average for intensities bias
    out_bias_correct = unifize(in_file=mean_aligned_filename,
                               outputtype='NIFTI_GZ', environ=environ)
    unbiased_mean_func_filename = out_bias_correct.outputs.out_file

    # Bias correct the antomical image
    out_unifize = unifize(in_file=anat_filename, outputtype='NIFTI_GZ',
                          environ=environ)
    unbiased_anat_filename = out_unifize.outputs.out_file

    # Update outputs
    output_files.extend([mean_aligned_filename, unbiased_mean_func_filename,
                         unbiased_anat_filename])

    #############################################
    # Rigid-body registration anat -> mean func #
    #############################################
    if prior_rigid_body_registration:
        anat_brain_filename = _compute_brain_mask(
            unbiased_anat_filename, write_dir, brain_volume, caching=False,
            terminal_output=terminal_output, environ=environ)
        func_brain_filename = _compute_brain_mask(
            unbiased_mean_func_filename, write_dir, brain_volume,
            caching=False, terminal_output=terminal_output, environ=environ)
        allineated_anat_filename, rigid_transform_file = \
            _rigid_body_anat_to_func(unbiased_anat_filename,
                                     anat_brain_filename,
                                     unbiased_mean_func_filename,
                                     func_brain_filename, write_dir,
                                     brain_volume, caching=False,
                                     terminal_output=terminal_output,
                                     environ=environ)
        output_files.append(allineated_anat_filename)
    else:
        allineated_anat_filename = unbiased_anat_filename

    ############################################
    # Nonlinear registration anat -> mean func #
    ############################################
    registered_anat_oblique_filename, mat_filename, warp_output_files =\
        _warp_anat_to_func(allineated_anat_filename,
                           unbiased_mean_func_filename,
                           write_dir, caching=caching,
                           terminal_output=terminal_output,
                           environ=environ)
    # Concatenate all the anat to func tranforms
    output_files.extend(warp_output_files)
    transform_filename = fname_presuffix(registered_anat_oblique_filename,
                                         suffix='.aff12.1D',
                                         use_ext=False)
    if prior_rigid_body_registration:
        _ = catmatvec(in_file=[(mat_filename, 'ONELINE'),
                               (rigid_transform_file, 'ONELINE')],
                      oneline=True,
                      out_file=transform_filename)
    else:
        _ = catmatvec(in_file=[(mat_filename, 'ONELINE')],
                      oneline=True,
                      out_file=transform_filename)

    ##################################################
    # Per-slice non-linear registration func -> anat #
    ##################################################
    warped_mean_func_filename, warp_filenames, warped_func_filename =\
        _per_slice_qwarp_func_to_anat(unbiased_mean_func_filename,
                                      registered_anat_oblique_filename,
                                      write_dir,
                                      voxel_size_x, voxel_size_y,
                                      apply_to_filename=allineated_filename,
                                      caching=caching,
                                      terminal_output=terminal_output,
                                      environ=environ)
    # Update the outputs
    output_files.append(warped_mean_func_filename)
    if not caching:
        for out_file in output_files:
            os.remove(out_file)

    # Update the fmri data
    setattr(session_data, "coreg_func_", warped_func_filename)
    setattr(session_data, "coreg_anat_", registered_anat_oblique_filename)
    setattr(session_data, "coreg_transform_", transform_filename)
    os.chdir(current_dir)


def fmri_sessions_to_template(sessions, t_r, head_template_filename,
                              write_dir,
                              brain_volume, brain_template_filename=None,
                              dilated_head_mask_filename=None,
                              prior_rigid_body_registration=False,
                              slice_timing=True,
                              maxlev=2,
                              caching=False, verbose=True):
    """ Registration of subject's functional and anatomical images to
    a given template.

    Parameters
    ----------
    sessions : sequence of sammba.registration.FMRISession
        fMRI sessions to register.

    t_r : float
        Repetition time for the EPI, in seconds.

    head_template_filename : str
        Template to register the functional to.

    brain_volume : int
        Volumes of the brain as passed to Rats_MM brain extraction tool.
        Typically 400 for mouse and 1800 for rat.

    write_dir : str
        Path to the affine 1D transform from anatomical to template space.

    brain_template_filename : str, optional
        Path to a brain template, passed to
        sammba.registration.anats_to_template

    dilated_head_mask_filename : str, optional
        Path to a dilated head mask, passed to
        sammba.registration.anats_to_template

    maxlev : int or None, optional
        Maximal level for the warp when registering anat to template. Passed to
        sammba.registration.anats_to_template

    caching : bool, optional
        Wether or not to use caching.

    verbose : bool, optional
        If True, all steps are verbose. Note that caching implies some
        verbosity in any case.

    Returns
    -------
    the same sequence with each animal_data updated: the following attributes
    are added
        - `template_` : str
                       Path to the given registration template.
        - `output_dir_` : str
                          Path to the output directory for each animal.
        - `coreg_func_` : str
                          Path to the coregistered functional image.
        - `coreg_anat_` : str
                          Path to the coregistered anatomical image.
        - `coreg_transform_` : str
                               Path to the transform from anat to func.
        - `registered_func_` : str
                               Path to the functional registered to template.
        - `registered_anat_` : str
                               Path to the anatomical registered to template.

    See also
    --------
    sammba.registration.coregister_fmri_session,
    sammba.registration.anats_to_template
    """
    if not hasattr(sessions, "__iter__"):
            raise ValueError(
                "'animals_data' input argument must be an iterable. You "
                "provided {0}".format(sessions.__class__))

    for n, data in enumerate(sessions):
        if not isinstance(data, FMRISession):
            raise ValueError('Each animal data must have type '
                             'sammba.registration.Animal. You provided {0} of'
                             ' type {1}'.format(data, type(data)))
        if data.animal_id is None:
            setattr(data, "animal_id", 'animal{0:03d}'.format(n + 1))

    # Check that ids are different
    animals_ids = [data.animal_id for data in sessions]
    if len(set(animals_ids)) != len(animals_ids):
        raise ValueError('Animals ids must be different. You'
                         ' provided {0}'.format(animals_ids))
    for n, animal_data in enumerate(sessions):
        animal_data._check_inputs()
        animal_output_dir = os.path.join(os.path.abspath(write_dir),
                                         animal_data.animal_id)
        animal_data._set_output_dir_(animal_output_dir)
        # XXX do a function for creating new attributes ?
        setattr(animal_data, "template_", head_template_filename)

        coregister_fmri_session(
            animal_data, t_r, write_dir, brain_volume,
            prior_rigid_body_registration=prior_rigid_body_registration,
            slice_timing=slice_timing,
            caching=caching, verbose=verbose)
        sessions[n] = animal_data

    anat_filenames = [animal_data.anat for animal_data in sessions]
    anats_registration = anats_to_template(
        anat_filenames,
        head_template_filename,
        animal_data.output_dir_,
        brain_volume,
        brain_template_filename=brain_template_filename,
        dilated_head_mask_filename=dilated_head_mask_filename,
        maxlev=maxlev,
        caching=caching, verbose=verbose)
    for n, (animal_data, normalized_anat_filename) in enumerate(zip(
            sessions, anats_registration.registered)):
        setattr(animal_data, "registered_anat_", normalized_anat_filename)
        sessions[n] = animal_data

    for n, (animal_data, anat_to_template_oned_filename,
            anat_to_template_warp_filename) in enumerate(zip(sessions,
                anats_registration.pre_transforms,
                anats_registration.transforms)):
        normalized_func_filename = _func_to_template(
            animal_data.coreg_func_,
            head_template_filename,
            animal_data.output_dir_,
            animal_data.coreg_transform_,
            anat_to_template_oned_filename,
            anat_to_template_warp_filename,
            caching=caching, verbose=verbose)

        setattr(animal_data, "registered_func_", normalized_func_filename)
        sessions[n] = animal_data

    return sessions
