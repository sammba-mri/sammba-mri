import os
import numpy as np
import nibabel
from sammba.externals.nipype.caching import Memory
from sammba.externals.nipype.interfaces import afni, ants, fsl
from sammba.externals.nipype.utils.filemanip import fname_presuffix
from sammba.interfaces import segmentation
from .utils import fix_obliquity
from .fmri_session import FMRISession
from .struct import anats_to_template


def coregister_fmri_session(session_data, t_r, write_dir, brain_volume,
                            use_rats_tool=True,
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
        Volume of the brain used for brain extraction.
        Typically 400 for mouse and 1800 for rat.

    use_rats_tool : bool, optional
        If True, brain mask is computed using RATS Mathematical Morphology.
        Otherwise, a histogram-based brain segmentation is used.

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

    Notes
    -----
    If `use_rats_tool` is turned on, RATS tool is used for brain extraction
    and has to be cited. For more information, see
    `RATS <http://www.iibi.uiowa.edu/content/rats-overview/>`_
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

    if use_rats_tool:
        if segmentation.Info().version() is None:
            raise ValueError('Can not locate RATS')
        else:
            ComputeMask = segmentation.MathMorphoMask
    else:
        ComputeMask = segmentation.HistogramMask

    if ants.base.Info().version is None:
        raise ValueError('Can not locate ANTS')

    if caching:
        memory = Memory(write_dir)
        tshift = memory.cache(afni.TShift)
        clip_level = memory.cache(afni.ClipLevel)
        volreg = memory.cache(afni.Volreg)
        allineate = memory.cache(afni.Allineate)
        tstat = memory.cache(afni.TStat)
        compute_mask = memory.cache(ComputeMask)
        calc = memory.cache(afni.Calc)
        allineate = memory.cache(afni.Allineate)
        allineate2 = memory.cache(afni.Allineate)
        unifize = memory.cache(afni.Unifize)
        bias_correct = memory.cache(ants.N4BiasFieldCorrection)
        catmatvec = memory.cache(afni.CatMatvec)
        warp = memory.cache(afni.Warp)
        resample = memory.cache(afni.Resample)
        slicer = memory.cache(afni.ZCutUp)
        warp_apply = memory.cache(afni.NwarpApply)
        qwarp = memory.cache(afni.Qwarp)
        merge = memory.cache(afni.Zcat)
        copy_geom = memory.cache(fsl.CopyGeom)
        overwrite = False
        for step in [tshift, volreg, allineate, allineate2,
                     tstat, compute_mask, calc, unifize, resample,
                     slicer, warp_apply, qwarp, merge]:
            step.interface().set_default_terminal_output(terminal_output)
    else:
        tshift = afni.TShift(terminal_output=terminal_output).run
        clip_level = afni.ClipLevel().run
        volreg = afni.Volreg(terminal_output=terminal_output).run
        allineate = afni.Allineate(terminal_output=terminal_output).run
        allineate2 = afni.Allineate(terminal_output=terminal_output).run  # TODO: remove after fixed bug
        tstat = afni.TStat(terminal_output=terminal_output).run
        compute_mask = ComputeMask(terminal_output=terminal_output).run
        calc = afni.Calc(terminal_output=terminal_output).run
        unifize = afni.Unifize(terminal_output=terminal_output).run
        bias_correct = ants.N4BiasFieldCorrection(terminal_output=terminal_output).run
        catmatvec = afni.CatMatvec().run
        warp = afni.Warp().run
        resample = afni.Resample(terminal_output=terminal_output).run
        slicer = afni.ZCutUp(terminal_output=terminal_output).run
        warp_apply = afni.NwarpApply(terminal_output=terminal_output).run
        qwarp = afni.Qwarp(terminal_output=terminal_output).run
        merge = afni.Zcat(terminal_output=terminal_output).run
        copy_geom = fsl.CopyGeom(terminal_output=terminal_output).run
        overwrite = True

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
        output_files.append(func_filename)

    ################################################
    # Register functional volumes to the first one #
    ################################################
    # XXX why do you need a thresholded image ?
    out_clip_level = clip_level(in_file=func_filename)
    out_calc_threshold = calc(
        in_file_a=func_filename,
        expr='ispositive(a-{0}) * a'.format(out_clip_level.outputs.clip_val),
        outputtype='NIFTI_GZ')
    thresholded_filename = out_calc_threshold.outputs.out_file

    out_volreg = volreg(  # XXX dfile not saved
        in_file=thresholded_filename,
        outputtype='NIFTI_GZ',
        environ=environ,
        oned_file=fname_presuffix(thresholded_filename,
                                  suffix='Vr.1Dfile.1D', use_ext=False),
        oned_matrix_save=fname_presuffix(thresholded_filename,
                                         suffix='Vr.aff12.1D', use_ext=False))

    # Apply the registration to the whole head
    out_allineate = allineate(in_file=func_filename,
                              master=func_filename,
                              in_matrix=out_volreg.outputs.oned_matrix_save,
                              out_file=fname_presuffix(func_filename,
                                                       suffix='Av'),
                              environ=environ)


    # 3dAllineate removes the obliquity. This is not a good way to readd it as
    # removes motion correction info in the header if it were an AFNI file...as
    # it happens it's NIfTI which does not store that so irrelevant!
    out_copy_geom = copy_geom(dest_file=out_allineate.outputs.out_file,
                              in_file=out_volreg.outputs.out_file)

    allineated_filename = out_copy_geom.outputs.out_file

    # Create a (hopefully) nice mean image for use in the registration
    out_tstat = tstat(in_file=allineated_filename, args='-mean',
                      outputtype='NIFTI_GZ', environ=environ)

    # Update outputs
    output_files.extend([thresholded_filename,
                         out_volreg.outputs.oned_matrix_save,
                         out_volreg.outputs.out_file,
                         out_volreg.outputs.md1d_file,
                         allineated_filename,
                         out_tstat.outputs.out_file])

    ###########################################
    # Corret anat and func for intensity bias #
    ###########################################
    # Correct the functional average for intensities bias
    out_bias_correct = bias_correct(input_image=out_tstat.outputs.out_file)
    unbiased_func_filename = out_bias_correct.outputs.output_image

    # Bias correct the antomical image
    out_unifize = unifize(in_file=anat_filename, outputtype='NIFTI_GZ',
                          environ=environ)
    unbiased_anat_filename = out_unifize.outputs.out_file

    # Update outputs
    output_files.extend([unbiased_func_filename, unbiased_anat_filename])

    #############################################
    # Rigid-body registration anat -> mean func #
    #############################################
    if prior_rigid_body_registration:
        # Mask the mean functional volume outside the brain.
        out_clip_level = clip_level(in_file=unbiased_func_filename)
        out_compute_mask_func = compute_mask(
            in_file=unbiased_func_filename,
            volume_threshold=brain_volume,
            intensity_threshold=int(out_clip_level.outputs.clip_val))
        out_cacl_func = calc(in_file_a=unbiased_func_filename,
                             in_file_b=out_compute_mask_func.outputs.out_file,
                             expr='a*b',
                             outputtype='NIFTI_GZ',
                             environ=environ)

        # Mask the anatomical volume outside the brain.
        out_clip_level = clip_level(in_file=unbiased_anat_filename)
        out_compute_mask_anat = compute_mask(
            in_file=unbiased_anat_filename,
            volume_threshold=brain_volume,
            intensity_threshold=int(out_clip_level.outputs.clip_val))
        out_cacl_anat = calc(in_file_a=unbiased_anat_filename,
                             in_file_b=out_compute_mask_anat.outputs.out_file,
                             expr='a*b',
                             outputtype='NIFTI_GZ',
                             environ=environ)

        # Compute the transformation from functional to anatomical brain
        # XXX: why in this sense
        out_allineate = allineate2(
            in_file=out_cacl_func.outputs.out_file,
            reference=out_cacl_anat.outputs.out_file,
            out_matrix=fname_presuffix(out_cacl_func.outputs.out_file,
                                       suffix='_shr.aff12.1D',
                                       use_ext=False),
            center_of_mass='',
            warp_type='shift_rotate',
            out_file=fname_presuffix(out_cacl_func.outputs.out_file,
                                     suffix='_shr'),
            environ=environ)
        rigid_transform_file = out_allineate.outputs.out_matrix
        output_files.extend([out_compute_mask_func.outputs.out_file,
                             out_cacl_func.outputs.out_file,
                             out_compute_mask_anat.outputs.out_file,
                             out_cacl_anat.outputs.out_file,
                             rigid_transform_file,
                             out_allineate.outputs.out_file])

        # apply the inverse transform to register the anatomical to the func
        catmatvec_out_file = fname_presuffix(rigid_transform_file,
                                             suffix='INV')
        out_catmatvec = catmatvec(in_file=[(rigid_transform_file, 'I')],
                                           oneline=True,
                                  out_file=catmatvec_out_file)
        output_files.append(out_catmatvec.outputs.out_file)
        out_allineate = allineate(
            in_file=unbiased_anat_filename,
            master=unbiased_func_filename,
            in_matrix=out_catmatvec.outputs.out_file,
            out_file=fname_presuffix(unbiased_anat_filename,
                                     suffix='_shr_in_func_space'),
            environ=environ)
        allineated_anat_filename = out_allineate.outputs.out_file
        output_files.append(allineated_anat_filename)
    else:
        allineated_anat_filename = unbiased_anat_filename

    ############################################
    # Nonlinear registration anat -> mean func #
    ############################################
    # 3dWarp doesn't put the obliquity in the header, so do it manually
    # This step generates one file per slice and per time point, so we are
    # making sure they are removed at the end
    out_warp = warp(in_file=allineated_anat_filename,
                    oblique_parent=unbiased_func_filename,
                    interp='quintic',
                    gridset=unbiased_func_filename,
                    outputtype='NIFTI_GZ',
                    verbose=True,
                    environ=environ)
    registered_anat_filename = out_warp.outputs.out_file
    registered_anat_oblique_filename = fix_obliquity(
        registered_anat_filename, unbiased_func_filename,
        overwrite=overwrite, verbose=verbose)

    # Concatenate all the anat to func tranforms
    mat_filename = fname_presuffix(registered_anat_filename,
                                   suffix='_warp.mat', use_ext=False)
    if not os.path.isfile(mat_filename):
        np.savetxt(mat_filename, [out_warp.runtime.stdout], fmt='%s')
        output_files.append(mat_filename)

    transform_filename = fname_presuffix(registered_anat_filename,
                                         suffix='_anat_to_func.aff12.1D',
                                         use_ext=False)
    if False:#prior_rigid_body_registration:
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
    # Slice anatomical image
    anat_img = nibabel.load(registered_anat_oblique_filename)
    anat_n_slices = anat_img.header.get_data_shape()[2]
    sliced_registered_anat_filenames = []
    for slice_n in range(anat_n_slices):
        out_slicer = slicer(in_file=registered_anat_oblique_filename,
                            keep='{0} {0}'.format(slice_n),
                            out_file=fname_presuffix(
                                registered_anat_oblique_filename,
                                suffix='Sl%d' % slice_n),
                            environ=environ)
        _ = fix_obliquity(out_slicer.outputs.out_file,
                          registered_anat_oblique_filename,
                          overwrite=overwrite, verbose=verbose)
        sliced_registered_anat_filenames.append(out_slicer.outputs.out_file)

    # Slice mean functional
    sliced_bias_corrected_filenames = []
    img = nibabel.load(func_filename)
    n_slices = img.header.get_data_shape()[2]
    for slice_n in range(n_slices):
        out_slicer = slicer(in_file=unbiased_func_filename,
                            keep='{0} {0}'.format(slice_n),
                            out_file=fname_presuffix(unbiased_func_filename,
                                                     suffix='Sl%d' % slice_n),
                            environ=environ)
        _ = fix_obliquity(out_slicer.outputs.out_file,
                          unbiased_func_filename,
                          overwrite=overwrite, verbose=verbose)
        sliced_bias_corrected_filenames.append(out_slicer.outputs.out_file)

    # Below line is to deal with slices where there is no signal (for example
    # rostral end of some anatomicals)

    # The inverse warp frequently fails, Resampling can help it work better
    # XXX why specifically .1 in voxel_size ?
    voxel_size_z = anat_img.header.get_zooms()[2]
    resampled_registered_anat_filenames = []
    for sliced_registered_anat_filename in sliced_registered_anat_filenames:
        out_resample = resample(in_file=sliced_registered_anat_filename,
                                voxel_size=(voxel_size_x, voxel_size_y,
                                            voxel_size_z),
                                outputtype='NIFTI_GZ',
                                environ=environ)
        resampled_registered_anat_filenames.append(
            out_resample.outputs.out_file)

    resampled_bias_corrected_filenames = []
    for sliced_bias_corrected_filename in sliced_bias_corrected_filenames:
        out_resample = resample(in_file=sliced_bias_corrected_filename,
                                voxel_size=(voxel_size_x, voxel_size_y,
                                            voxel_size_z),
                                outputtype='NIFTI_GZ',
                                environ=environ)
        resampled_bias_corrected_filenames.append(
            out_resample.outputs.out_file)

    # single slice non-linear functional to anatomical registration
    warped_slices = []
    warp_filenames = []
    for (resampled_bias_corrected_filename,
         resampled_registered_anat_filename) in zip(
            resampled_bias_corrected_filenames,
            resampled_registered_anat_filenames):
        warped_slice = fname_presuffix(resampled_bias_corrected_filename,
                                       suffix='_qw')
        out_qwarp = qwarp(in_file=resampled_bias_corrected_filename,
                          base_file=resampled_registered_anat_filename,
                          iwarp=True,  # XXX: is this necessary
                          noneg=True,
                          blur=[0],
                          nmi=True,
                          noXdis=True,
                          allineate=True,
                          allineate_opts='-parfix 1 0 -parfix 2 0 -parfix 3 0 '
                                         '-parfix 4 0 -parfix 5 0 -parfix 6 0 '
                                         '-parfix 7 0 -parfix 9 0 '
                                         '-parfix 10 0 -parfix 12 0',
                          out_file=warped_slice,
                          environ=environ)
        warped_slices.append(out_qwarp.outputs.warped_source)
        warp_filenames.append(out_qwarp.outputs.source_warp)
        output_files.append(out_qwarp.outputs.base_warp)
        # There are files geenrated by the allineate option
        output_files.extend([
            fname_presuffix(out_qwarp.outputs.warped_source, suffix='_Allin'),
            fname_presuffix(out_qwarp.outputs.warped_source,
                            suffix='_Allin.nii', use_ext=False),
            fname_presuffix(out_qwarp.outputs.warped_source,
                            suffix='_Allin.aff12.1D', use_ext=False)])

    # Resample the mean volume back to the initial resolution,
    voxel_size = nibabel.load(
        unbiased_func_filename).header.get_zooms()
    resampled_warped_slices = []
    for warped_slice in warped_slices:
        out_resample = resample(in_file=warped_slice,
                                voxel_size=voxel_size,
                                outputtype='NIFTI_GZ',
                                environ=environ)
        resampled_warped_slices.append(out_resample.outputs.out_file)

    # fix the obliquity
    for (sliced_registered_anat_filename, resampled_warped_slice) in zip(
            sliced_registered_anat_filenames, resampled_warped_slices):
        _ = fix_obliquity(resampled_warped_slice,
                          sliced_registered_anat_filename,
                          overwrite=overwrite, verbose=verbose)

    # slice functional
    sliced_func_filenames = []
    for slice_n in range(n_slices):
        out_slicer = slicer(in_file=allineated_filename,
                            keep='{0} {0}'.format(slice_n),
                            out_file=fname_presuffix(allineated_filename,
                                                     suffix='Sl%d' % slice_n),
                            environ=environ)
        _ = fix_obliquity(out_slicer.outputs.out_file,
                          allineated_filename,
                          overwrite=overwrite, verbose=verbose)
        sliced_func_filenames.append(out_slicer.outputs.out_file)

    # Apply the precomputed warp slice by slice
    warped_func_slices = []
    for (sliced_func_filename, warp_filename) in zip(
            sliced_func_filenames, warp_filenames):
        out_warp_apply = warp_apply(in_file=sliced_func_filename,
                                    master=sliced_func_filename,
                                    warp=warp_filename,
                                    out_file=fname_presuffix(
                                        sliced_func_filename, suffix='_qw'),
                                    environ=environ)
        warped_func_slices.append(out_warp_apply.outputs.out_file)

    # Finally, merge all slices !
    out_merge_func = merge(in_files=warped_func_slices,
                           outputtype='NIFTI_GZ',
                           environ=environ)

    # Fix the obliquity
    _ = fix_obliquity(out_merge_func.outputs.out_file, allineated_filename,
                      overwrite=overwrite, verbose=verbose)

    # Update the fmri data
    setattr(session_data, "coreg_func_", out_merge_func.outputs.out_file)
    setattr(session_data, "coreg_anat_", registered_anat_oblique_filename)
    setattr(session_data, "coreg_transform_", transform_filename)
    os.chdir(current_dir)

    # Collect the outputs
    output_files.extend(sliced_registered_anat_filenames +
                        sliced_bias_corrected_filenames +
                        resampled_registered_anat_filenames +
                        resampled_bias_corrected_filenames +
                        warped_slices + warp_filenames +
                        resampled_warped_slices +
                        sliced_func_filenames +
                        warped_func_slices)
    if not caching:
        for out_file in output_files:
            if os.path.isfile(out_file):
                os.remove(out_file)


def _func_to_template(func_coreg_filename, template_filename, write_dir,
                      func_to_anat_oned_filename,
                      anat_to_template_oned_filename,
                      anat_to_template_warp_filename,
                      voxel_size=None,
                      caching=False, verbose=True):
    """ Applies successive transforms to coregistered functional to put it in
    template space.

    Parameters
    ----------
    coreg_func_filename : str
        Path to functional volume, coregistered to a common space with the
        anatomical volume.

    template_filename : str
        Template to register the functional to.

    func_to_anat_oned_filename : str
        Path to the affine 1D transform from functional to coregistration
        space.

    anat_to_template_oned_filename : str
        Path to the affine 1D transform from anatomical to template space.

    anat_to_template_warp_filename : str
        Path to the warp transform from anatomical to template space.

    voxel_size : 3-tuple of floats, optional
        Voxel size of the registered functional, in mm.

    caching : bool, optional
        Wether or not to use caching.

    verbose : bool, optional
        If True, all steps are verbose. Note that caching implies some
        verbosity in any case.
    """
    environ = {}
    if verbose:
        terminal_output = 'allatonce'
    else:
        terminal_output = 'none'

    if caching:
        memory = Memory(write_dir)
        resample = memory.cache(afni.Resample)
        catmatvec = memory.cache(afni.CatMatvec)
        allineate = memory.cache(afni.Allineate)
        warp_apply = memory.cache(afni.NwarpApply)
        for step in [resample, allineate, warp_apply]:
            step.interface().set_default_terminal_output(terminal_output)
    else:
        resample = afni.Resample(terminal_output=terminal_output).run
        catmatvec = afni.CatMatvec().run
        allineate = afni.Allineate(terminal_output=terminal_output).run
        warp_apply = afni.NwarpApply(terminal_output=terminal_output).run
        environ['AFNI_DECONFLICT'] = 'OVERWRITE'

    current_dir = os.getcwd()
    os.chdir(write_dir)  # XXX to remove
    normalized_filename = fname_presuffix(func_coreg_filename,
                                          suffix='_normalized')
    if voxel_size is None:
        func_template_filename = template_filename
    else:
        out_resample = resample(in_file=template_filename,
                                voxel_size=voxel_size,
                                outputtype='NIFTI_GZ',
                                environ=environ)
        func_template_filename = out_resample.outputs.out_file

    if anat_to_template_warp_filename is None:
        affine_transform_filename = fname_presuffix(func_to_anat_oned_filename,
                                                    suffix='_to_template')
        _ = catmatvec(in_file=[(anat_to_template_oned_filename, 'ONELINE'),
                               (func_to_anat_oned_filename, 'ONELINE')],
                      oneline=True,
                      out_file=affine_transform_filename,
                      environ=environ)
        _ = allineate(
            in_file=func_coreg_filename,
            master=func_template_filename,
            in_matrix=affine_transform_filename,
            out_file=normalized_filename,
            environ=environ)
    else:
        warp = "'{0} {1} {2}'".format(anat_to_template_warp_filename,
                                      anat_to_template_oned_filename,
                                      func_to_anat_oned_filename)
    
        _ = warp_apply(in_file=func_coreg_filename,
                       master=func_template_filename,
                       warp=warp,
                       out_file=normalized_filename,
                       environ=environ)
    os.chdir(current_dir)
    return normalized_filename


def fmri_sessions_to_template(sessions, t_r, head_template_filename,
                              write_dir,
                              brain_volume, use_rats_tool=True,
                              brain_template_filename=None,
                              dilated_head_mask_filename=None,
                              prior_rigid_body_registration=False,
                              slice_timing=True,
                              registration_kind='affine',
                              maxlev=2,
                              func_voxel_size=None,
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
        Volume of the brain used for brain extraction.
        Typically 400 for mouse and 1800 for rat.

    use_rats_tool : bool, optional
        If True, brain mask is computed using RATS Mathematical Morphology.
        Otherwise, a histogram-based brain segmentation is used.

    write_dir : str
        Path to the affine 1D transform from anatomical to template space.

    brain_template_filename : str, optional
        Path to a brain template, passed to
        sammba.registration.anats_to_template

    dilated_head_mask_filename : str, optional
        Path to a dilated head mask, passed to
        sammba.registration.anats_to_template

    registration_kind : one of {'rigid', 'affine', 'nonlinear'}, optional
        The allowed transform kind from anat to template.

    maxlev : int or None, optional
        Maximal level for the warp when registering anat to template. Passed to
        sammba.registration.anats_to_template

    func_voxel_size : 3-tuple of floats, optional
        Voxel size of the registered functional, in mm.

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

    Notes
    -----
    If `use_rats_tool` is turned on, RATS tool is used for brain extraction
    and has to be cited. For more information, see
    `RATS <http://www.iibi.uiowa.edu/content/rats-overview/>`_
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
            use_rats_tool=use_rats_tool,
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
        use_rats_tool=use_rats_tool,
        brain_template_filename=brain_template_filename,
        dilated_head_mask_filename=dilated_head_mask_filename,
        registration_kind=registration_kind,
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
            voxel_size=func_voxel_size,
            caching=caching, verbose=verbose)

        setattr(animal_data, "registered_func_", normalized_func_filename)
        sessions[n] = animal_data

    return sessions
