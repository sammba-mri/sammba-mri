import os
import numpy as np
import nibabel
from nilearn.image.resampling import coord_transform
from ..externals.nipype.caching import Memory
from ..externals.nipype.interfaces import afni, ants, fsl
from ..externals.nipype.utils.filemanip import fname_presuffix
from ..interfaces import segmentation
from .utils import fix_obliquity, compute_n4_max_shrink


def mask_report(mask_file, expected_volume):
    """ Outputs the mask

    Parameters
    ----------
    mask_file : str
        Path to the mask image

    expected_volume : float
        Expected volume in the mask.
    """
    # TODO: symmetry, length and width
    mask_img = nibabel.load(mask_file)
    volume = segmentation.compute_volume(mask_img)
    volume_accuracy = volume / expected_volume * 100.

    mask_data = mask_img.get_data()
    i, j, k = np.where(mask_data != 0)
    voxels_coords = np.array(coord_transform(i, j, k, mask_img.affine)).T
    x_range, y_range, z_range = voxels_coords.max(axis=0) - \
        voxels_coords.min(axis=0)
    return "extracted volume is {0:0.1f}% of expected volume, sizes "\
           "x: {1:0.2f}, y: {2:0.2f}, z: {3:0.2f}".format(volume_accuracy,
                                                          x_range, y_range,
                                                          z_range)


def compute_brain_mask(head_file, brain_volume, write_dir, bias_correct=True,
                       caching=False,
                       terminal_output='allatonce',
                       lower_cutoff=.2, upper_cutoff=.85, closing=0,
                       connected=True, dilation_size=(1, 1, 2), opening=5,
                       use_rats_tool=True, **unifize_kwargs):
    """
    Parameters
    ----------
    brain_volume : int
        Volume of the brain used for brain extraction.
        Typically 400 for mouse and 1800 for rat.

    use_rats_tool : bool, optional
        If True, brain mask is computed using RATS Mathematical Morphology.
        Otherwise, a histogram-based brain segmentation is used.

    caching : bool, optional
        Wether or not to use caching.

    unifize_kwargs : dict, optional
        Is passed to sammba.externals.nipype.interfaces.afni.Unifize.

    Returns
    -------
    path to brain extracted image.

    Notes
    -----
    If `use_rats_tool` is turned on, RATS tool is used for brain extraction
    and has to be cited. For more information, see
    `RATS <http://www.iibi.uiowa.edu/content/rats-overview/>`_
    """
    if use_rats_tool:
        if segmentation.Info().version() is None:
            raise ValueError('Can not locate Rats')
        else:
            ComputeMask = segmentation.MathMorphoMask
            compute_mask_args = {}
    else:
        ComputeMask = segmentation.HistogramMask
        compute_mask_args = {'lower_cutoff': lower_cutoff,
                             'upper_cutoff': upper_cutoff,
                             'closing': closing,
                             'connected': connected,
                             'dilation_size': dilation_size,
                             'opening': opening}

    environ = {'AFNI_DECONFLICT': 'OVERWRITE'}

    if caching:
        memory = Memory(write_dir)
        clip_level = memory.cache(afni.ClipLevel)
        compute_mask = memory.cache(ComputeMask)
        unifize = memory.cache(afni.Unifize)
        for step in [compute_mask, unifize]:
            step.interface().set_default_terminal_output(terminal_output)
    else:
        clip_level = afni.ClipLevel().run
        compute_mask = ComputeMask(terminal_output=terminal_output).run
        unifize = afni.Unifize(terminal_output=terminal_output).run

    if bias_correct:
        if unifize_kwargs is None:
            unifize_kwargs = {}

        file_to_mask = _afni_bias_correct(
            head_file, write_dir,
            out_file=fname_presuffix(head_file,
                                     suffix='_unifized_for_extraction',
                                     newpath=write_dir),
            caching=caching,
            terminal_output=terminal_output,
            environ=environ, **unifize_kwargs)
    else:
        file_to_mask = head_file

    out_clip_level = clip_level(in_file=file_to_mask)
    out_compute_mask = compute_mask(
        in_file=file_to_mask,
        out_file=fname_presuffix(file_to_mask,
                                 suffix='_brain_mask',
                                 newpath=write_dir),
        volume_threshold=brain_volume,
        intensity_threshold=int(out_clip_level.outputs.clip_val),
        **compute_mask_args)

    # Remove intermediate output
    if not caching and bias_correct:
        os.remove(file_to_mask)

    return out_compute_mask.outputs.out_file


def _apply_mask(head_file, mask_file, write_dir=None,
                caching=False, terminal_output='allatonce'):
    """
    Parameters
    ----------
    head_file : str
        Path to the image to mask.

    mask_file : str
        Path to the image mask to apply.

    write_dir : str or None, optional
        Path to the folder to save the output file to. If None, the folder
        of the head file is used.

    caching : bool, optional
        Wether or not to use caching.

    terminal_output : one of {'stream', 'allatonce', 'file', 'none'}
        Control terminal output :
            'stream' : displays to terminal immediately,
            'allatonce' : waits till command is finished to display output,
            'file' : writes output to file
            'none' : output is ignored

    Returns
    -------
    path to brain extracted image.
    """
    if write_dir is None:
        write_dir = os.path.dirname(head_file)

    if caching:
        memory = Memory(write_dir)
        apply_mask = memory.cache(fsl.ApplyMask)
        apply_mask.interface().set_default_terminal_output(terminal_output)
    else:
        apply_mask = fsl.ApplyMask(terminal_output=terminal_output).run

    # Check mask is binary
    mask_img = nibabel.load(mask_file)
    mask = mask_img.get_data()
    values = np.unique(mask)
    if len(values) == 2:
        # If there are 2 different values, one of them must be 0 (background)
        if not 0 in values:
            raise ValueError('Background of the mask must be represented with'
                             '0. Given mask contains: %s.' % values)
    elif len(values) != 2:
        # If there are more than 2 values, the mask is invalid
        raise ValueError('Given mask is not made of 2 values: %s'
                         '. Cannot interpret as true or false' % values)
    try:
        np.testing.assert_array_equal(nibabel.load(mask_file).affine,
                                      nibabel.load(head_file).affine)
    except AssertionError:
        raise ValueError('Given mask {0} and file {1} do not have the same '
                         'affine'.format(mask_file, head_file))

    out_apply_mask = apply_mask(in_file=head_file,
                                mask_file=mask_file,
                                out_file=fname_presuffix(head_file,
                                                         suffix='_masked',
                                                         newpath=write_dir))
    return out_apply_mask.outputs.out_file


def _delete_orientation(in_file, write_dir, min_zoom=.1, caching=False,
                        verbose=True):
    if verbose:
        terminal_output = 'stream'
    else:
        terminal_output = 'none'

    if caching:
        memory = Memory(write_dir)
        copy = memory.cache(afni.Copy)
        refit = memory.cache(afni.Refit)
        center_mass = memory.cache(afni.CenterMass)
        for step in [copy, refit]:
            step.interface().set_default_terminal_output(terminal_output)
    else:
        copy = afni.Copy(terminal_output=terminal_output).run
        refit = afni.Refit(terminal_output=terminal_output).run
        center_mass = afni.CenterMass().run

    out_copy = copy(in_file=in_file,
                    out_file=fname_presuffix(in_file, newpath=write_dir))

    zooms = nibabel.load(in_file).header.get_zooms()[:3]
    out_refit = refit(in_file=out_copy.outputs.out_file,
                      xyzscale=min_zoom / min(zooms))
    out_center_mass = center_mass(in_file=out_refit.outputs.out_file,
                                  cm_file=fname_presuffix(in_file,
                                                          suffix='.txt',
                                                          use_ext=False,
                                                          newpath=write_dir),
                                  set_cm=(0, 0, 0))
    return out_center_mass.outputs.out_file


def _ants_bias_correct(in_file, write_dir=None, caching=False,
                       terminal_output='allatonce', verbose=True,
                       environ=None):
    if write_dir is None:
        write_dir = os.path.dirname(in_file)

    if environ is None:
        environ = {'AFNI_DECONFLICT': 'OVERWRITE'}

    if caching:
        memory = Memory(write_dir)
        bias_correct = memory.cache(ants.N4BiasFieldCorrection)
        copy = memory.cache(afni.Copy)
        copy_geom = memory.cache(fsl.CopyGeom)
        bias_correct.interface().set_default_terminal_output(terminal_output)
        copy.interface().set_default_terminal_output(terminal_output)
    else:
        bias_correct = ants.N4BiasFieldCorrection(
            terminal_output=terminal_output).run
        copy = afni.Copy(terminal_output=terminal_output).run
        copy_geom = fsl.CopyGeom(terminal_output=terminal_output).run

    out_bias_correct = bias_correct(
        input_image=in_file,
        shrink_factor=compute_n4_max_shrink(in_file),
        verbose=verbose,
        output_image=fname_presuffix(in_file, suffix='_n4_deoblique',
                                     newpath=write_dir))
    out_copy = copy(
        in_file=out_bias_correct.outputs.output_image,
        out_file=fname_presuffix(in_file,
                                 suffix='_n4',
                                 newpath=write_dir),
        environ=environ,
        verb=verbose)

    out_copy_geom = copy_geom(dest_file=out_copy.outputs.out_file,
                              in_file=in_file)

    return out_copy_geom.outputs.out_file


def _afni_bias_correct(in_file, write_dir=None, out_file=None, caching=False,
                       terminal_output='allatonce', verbose=True,
                       environ=None,
                       **unifize_kwargs):
    if write_dir is None:
        write_dir = os.path.dirname(in_file)

    if environ is None:
        environ = {'AFNI_DECONFLICT': 'OVERWRITE'}

    if caching:
        memory = Memory(write_dir)
        copy_geom = memory.cache(fsl.CopyGeom)
        unifize = memory.cache(afni.Unifize)
        unifize.interface().set_default_terminal_output(terminal_output)
    else:
        copy_geom = fsl.CopyGeom(terminal_output=terminal_output).run
        unifize = afni.Unifize(terminal_output=terminal_output).run

    if out_file is None:
        out_file = fname_presuffix(in_file,
                                   suffix='_unifized',
                                   newpath=write_dir)
    out_unifize = unifize(in_file=in_file,
                          out_file=out_file,
                          environ=environ,
                          quiet=not(verbose),
                          **unifize_kwargs)
    if True:
        out_copy_geom = copy_geom(dest_file=out_unifize.outputs.out_file,
                                  in_file=in_file)

    return out_copy_geom.outputs.out_file


def _rigid_body_register(moving_head_file, reference_head_file,
                         moving_brain_file, reference_brain_file,
                         write_dir=None, verbose=True,
                         caching=False, terminal_output='allatonce',
                         environ=None):
    # XXX: add verbosity
    if write_dir is None:
        write_dir = os.path.dirname(moving_head_file)

    if environ is None:
        environ = {'AFNI_DECONFLICT': 'OVERWRITE'}

    if caching:
        memory = Memory(write_dir)
        allineate = memory.cache(afni.Allineate)
        allineate2 = memory.cache(afni.Allineate)
        catmatvec = memory.cache(afni.CatMatvec)
        for step in [allineate, allineate2]:
            step.interface().set_default_terminal_output(terminal_output)
    else:
        allineate = afni.Allineate(terminal_output=terminal_output).run
        allineate2 = afni.Allineate(terminal_output=terminal_output).run
        catmatvec = afni.CatMatvec().run

    output_files = [moving_brain_file, reference_brain_file]

    # Compute the transformation from functional to anatomical brain
    # XXX: why in this sense
    out_allineate = allineate2(
        in_file=reference_brain_file,
        reference=moving_brain_file,
        out_matrix=fname_presuffix(reference_brain_file,
                                   suffix='_shr.aff12.1D',
                                   use_ext=False,
                                   newpath=write_dir),
        center_of_mass='',
        warp_type='shift_rotate',
        out_file=fname_presuffix(reference_brain_file, suffix='_shr'),
        environ=environ,
        verb=verbose)
    rigid_transform_file = out_allineate.outputs.out_matrix
    output_files.append(out_allineate.outputs.out_file)

    # apply the inverse transform to register the anatomical to the func
    out_catmatvec = catmatvec(
        in_file=[(rigid_transform_file, 'I')],
        oneline=True,
        out_file=fname_presuffix(rigid_transform_file, suffix='INV',
                                 newpath=write_dir),
        environ=environ)
    output_files.append(out_catmatvec.outputs.out_file)
    out_allineate_apply = allineate(
        in_file=moving_head_file,
        master=reference_head_file,
        in_matrix=out_catmatvec.outputs.out_file,
        out_file=fname_presuffix(moving_head_file, suffix='_shr',
                                 newpath=write_dir),
        environ=environ)

    # Remove intermediate output
    if not caching:
        for output_file in output_files:
            os.remove(output_file)

    return out_allineate_apply.outputs.out_file, rigid_transform_file


def _warp(to_warp_file, reference_file, write_dir=None, caching=False,
          terminal_output='allatonce', verbose=True, environ=None):
    if write_dir is None:
        write_dir = os.path.dirname(to_warp_file)

    if environ is None:
        environ = {'AFNI_DECONFLICT': 'OVERWRITE'}

    if caching:
        memory = Memory(write_dir)
        warp = memory.cache(afni.Warp)
    else:
        warp = afni.Warp().run

    out_warp = warp(in_file=to_warp_file,
                    oblique_parent=reference_file,
                    interp='quintic',
                    gridset=reference_file,
                    out_file=fname_presuffix(to_warp_file, suffix='_warped',
                                             newpath=write_dir),
                    verbose=True,
                    save_matfile=True,
                    environ=environ)

    # 3dWarp doesn't put the obliquity in the header, so do it manually
    warped_oblique_file = fix_obliquity(
        out_warp.outputs.out_file, reference_file,
        verbose=verbose, caching=caching,
        caching_dir=write_dir, environ=environ)

    return warped_oblique_file, out_warp.outputs.mat_file


def _per_slice_qwarp(to_qwarp_file, reference_file,
                     voxel_size_x, voxel_size_y, apply_to_file=None,
                     write_dir=None,
                     caching=False,
                     verbose=True, terminal_output='allatonce', environ=None):
    if write_dir is None:
        write_dir = os.path.dirname(to_qwarp_file),

    if environ is None:
        environ = {'AFNI_DECONFLICT': 'OVERWRITE'}

    if caching:
        memory = Memory(write_dir)
        resample = memory.cache(afni.Resample)
        slicer = memory.cache(fsl.Slice)
        warp_apply = memory.cache(afni.NwarpApply)
        qwarp = memory.cache(afni.Qwarp)
        merge = memory.cache(fsl.Merge)
        for step in [resample, slicer, warp_apply, qwarp, merge]:
            step.interface().set_default_terminal_output(terminal_output)
    else:
        resample = afni.Resample(terminal_output=terminal_output).run
        slicer = fsl.Slice(terminal_output=terminal_output).run
        warp_apply = afni.NwarpApply(terminal_output=terminal_output).run
        qwarp = afni.Qwarp(terminal_output=terminal_output).run
        merge = fsl.Merge(terminal_output=terminal_output).run

    # Slice anatomical image
    reference_img = nibabel.load(reference_file)
    per_slice_dir = os.path.join(write_dir, 'per_slice')
    if not os.path.isdir(per_slice_dir):
        os.makedirs(per_slice_dir)

    out_slicer = slicer(in_file=reference_file,
                        out_base_name=fname_presuffix(reference_file,
                                                      newpath=per_slice_dir,
                                                      use_ext=False))
    sliced_reference_files = out_slicer.outputs.out_files

    # Slice mean functional
    out_slicer = slicer(in_file=to_qwarp_file,
                        out_base_name=fname_presuffix(to_qwarp_file,
                                                      newpath=per_slice_dir,
                                                      use_ext=False))
    sliced_to_qwarp_files = out_slicer.outputs.out_files

    # Below line is to deal with slices where there is no signal (for example
    # rostral end of some anatomicals)

    # The inverse warp frequently fails, Resampling can help it work better
    # XXX why specifically .1 in voxel_size ?
    voxel_size_z = reference_img.header.get_zooms()[2]
    resampled_sliced_reference_files = []
    for sliced_reference_file in sliced_reference_files:
        out_resample = resample(
            in_file=sliced_reference_file,
            voxel_size=(voxel_size_x, voxel_size_y, voxel_size_z),
            out_file=fname_presuffix(sliced_reference_file,
                                     suffix='_resampled'),
            environ=environ)
        resampled_sliced_reference_files.append(out_resample.outputs.out_file)

    resampled_sliced_to_qwarp_files = []
    for sliced_to_qwarp_file in sliced_to_qwarp_files:
        out_resample = resample(
            in_file=sliced_to_qwarp_file,
            voxel_size=(voxel_size_x, voxel_size_y, voxel_size_z),
            out_file=fname_presuffix(sliced_to_qwarp_file,
                                     suffix='_resampled'),
            environ=environ)
        resampled_sliced_to_qwarp_files.append(
            out_resample.outputs.out_file)

    # single slice non-linear functional to anatomical registration
    warped_slices = []
    warp_files = []
    output_files = []
    resampled_sliced_to_qwarp_files_to_remove = []
    for (resampled_sliced_to_qwarp_file,
         resampled_sliced_reference_file) in zip(
            resampled_sliced_to_qwarp_files,
            resampled_sliced_reference_files):
        warped_slice = fname_presuffix(resampled_sliced_to_qwarp_file,
                                       suffix='_qwarped')
        to_qwarp_data = nibabel.load(resampled_sliced_to_qwarp_file).get_data()
        ref_data = nibabel.load(resampled_sliced_reference_file).get_data()

        if to_qwarp_data.max() == 0 or ref_data.max() == 0:
            # deal with slices where there is no signal
            warped_slices.append(resampled_sliced_to_qwarp_file)
            warp_files.append(None)
        else:
            resampled_sliced_to_qwarp_files_to_remove.append(resampled_sliced_to_qwarp_file)
            out_qwarp = qwarp(
                in_file=resampled_sliced_to_qwarp_file,
                base_file=resampled_sliced_reference_file,
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
                environ=environ,
                verb=verbose)
            # XXX fix qwarp bug : out_qwarp.outputs.warped_source extension is
            # +tlrc.HEAD if base_file and in_file are of different extensions
            warped_slices.append(warped_slice)
            warp_files.append(out_qwarp.outputs.source_warp)
            # There are files geenrated by the allineate option
            output_files.extend([
                fname_presuffix(out_qwarp.outputs.warped_source,
                                suffix='_Allin.nii', use_ext=False),
                fname_presuffix(out_qwarp.outputs.warped_source,
                                suffix='_Allin.aff12.1D', use_ext=False)])

    # Resample the mean volume back to the initial resolution,
    voxel_size = nibabel.load(to_qwarp_file).header.get_zooms()[:3]
    resampled_warped_slices = []
    for warped_slice in warped_slices:
        out_resample = resample(in_file=warped_slice,
                                voxel_size=voxel_size,
                                out_file=fname_presuffix(warped_slice,
                                                         suffix='_resampled'),
                                environ=environ)
        resampled_warped_slices.append(out_resample.outputs.out_file)

    # fix the obliquity
    oblique_resampled_warped_slices = []
    for (sliced_reference_file, resampled_warped_slice) in zip(
            sliced_reference_files, resampled_warped_slices):
        oblique_slice = fix_obliquity(resampled_warped_slice,
                                      sliced_reference_file,
                                      verbose=verbose,
                                      caching=caching,
                                      caching_dir=per_slice_dir,
                                      environ=environ)
        oblique_resampled_warped_slices.append(oblique_slice)

    out_merge_func = merge(
        in_files=oblique_resampled_warped_slices,
        dimension='z',
        merged_file=fname_presuffix(to_qwarp_file, suffix='_perslice',
                                    newpath=write_dir),
        environ=environ)

    # Fix the obliquity
    oblique_merged = fix_obliquity(out_merge_func.outputs.merged_file,
                                   reference_file,
                                   verbose=verbose,
                                   caching=caching, caching_dir=per_slice_dir,
                                   environ=environ)

    # Collect the outputs
    output_files.extend(sliced_reference_files +
                        sliced_to_qwarp_files +
                        resampled_sliced_reference_files +
                        resampled_sliced_to_qwarp_files_to_remove +
                        warped_slices + oblique_resampled_warped_slices)

    # Apply the precomputed warp slice by slice
    if apply_to_file is not None:
        # slice functional
        out_slicer = slicer(in_file=apply_to_file,
                            out_base_name=fname_presuffix(apply_to_file,
                                                          newpath=per_slice_dir,
                                                          use_ext=False))
        sliced_apply_to_files = out_slicer.outputs.out_files
        warped_apply_to_slices = []
        sliced_apply_to_files_to_remove = []
        for (sliced_apply_to_file, warp_file) in zip(
                sliced_apply_to_files, warp_files):
            if warp_file is None:
                warped_apply_to_slices.append(sliced_apply_to_file)
            else:
                sliced_apply_to_files_to_remove.append(sliced_apply_to_file)
                out_warp_apply = warp_apply(in_file=sliced_apply_to_file,
                                            master=sliced_apply_to_file,
                                            warp=warp_file,
                                            out_file=fname_presuffix(
                                                sliced_apply_to_file,
                                                suffix='_qwarped'),
                                            environ=environ)
                warped_apply_to_slices.append(out_warp_apply.outputs.out_file)

        # Fix the obliquity
        oblique_warped_apply_to_slices = []
        for (sliced_apply_to_file, warped_apply_to_slice) in zip(
                sliced_apply_to_files, warped_apply_to_slices):
            oblique_slice = fix_obliquity(warped_apply_to_slice,
                                          sliced_apply_to_file,
                                          verbose=verbose,
                                          caching=caching,
                                          caching_dir=per_slice_dir,
                                          environ=environ)
            oblique_warped_apply_to_slices.append(oblique_slice)

        # Finally, merge all slices !
        out_merge_apply_to = merge(
            in_files=oblique_warped_apply_to_slices,
            dimension='z',
            merged_file=fname_presuffix(apply_to_file, suffix='_perslice',
                                        newpath=write_dir),
            environ=environ)

        # Fix the obliquity
        merged_apply_to_file = fix_obliquity(
            out_merge_apply_to.outputs.merged_file, apply_to_file,
            verbose=verbose, caching=caching,
            caching_dir=per_slice_dir, environ=environ)

        # Update the outputs
        output_files.extend(sliced_apply_to_files_to_remove +
                            oblique_warped_apply_to_slices)
    else:
        merged_apply_to_file = None

    if not caching:
        for out_file in output_files:
            os.remove(out_file)

    return (oblique_merged, warp_files,
            merged_apply_to_file)


def _apply_perslice_warp(apply_to_file, warp_files,
                         voxel_size_x, voxel_size_y,
                         write_dir=None,
                         caching=False,
                         verbose=True, terminal_output='allatonce',
                         environ=None):

    # Apply the precomputed warp slice by slice
    if write_dir is None:
        write_dir = os.path.dirname(apply_to_file),

    if environ is None:
        environ = {'AFNI_DECONFLICT': 'OVERWRITE'}

    if caching:
        memory = Memory(write_dir)
        resample = memory.cache(afni.Resample)
        slicer = memory.cache(fsl.Slice)
        warp_apply = memory.cache(afni.NwarpApply)
        qwarp = memory.cache(afni.Qwarp)
        merge = memory.cache(fsl.Merge)
        for step in [resample, slicer, warp_apply, qwarp, merge]:
            step.interface().set_default_terminal_output(terminal_output)
    else:
        resample = afni.Resample(terminal_output=terminal_output).run
        slicer = fsl.Slice(terminal_output=terminal_output).run
        warp_apply = afni.NwarpApply(terminal_output=terminal_output).run
        qwarp = afni.Qwarp(terminal_output=terminal_output).run
        merge = fsl.Merge(terminal_output=terminal_output).run

    apply_to_img = nibabel.load(apply_to_file)
    n_slices = apply_to_img.header.get_data_shape()[2]

    if len(warp_files) != n_slices:
        raise ValueError('number of warp files {0} does not match number of '
                         'slices {1}'.format(len(warp_files), n_slices))
    # Slice anatomical image
    output_files = []
    per_slice_dir = os.path.join(write_dir, 'per_slice')
    if not os.path.isdir(per_slice_dir):
        os.makedirs(per_slice_dir)

    # slice functional
    sliced_apply_to_files = []
    out_slicer = slicer(in_file=apply_to_file,
                        out_base_name=fname_presuffix(apply_to_file,
                                                      newpath=per_slice_dir,
                                                      use_ext=False))
    sliced_apply_to_files = out_slicer.outputs.out_files

    warped_apply_to_slices = []
    sliced_apply_to_files_to_remove = []
    for (sliced_apply_to_file, warp_file) in zip(
            sliced_apply_to_files, warp_files):
        if warp_file is None:
            warped_apply_to_slices.append(sliced_apply_to_file)
        else:
            sliced_apply_to_files_to_remove.append(sliced_apply_to_file)
            out_warp_apply = warp_apply(in_file=sliced_apply_to_file,
                                        master=sliced_apply_to_file,
                                        warp=warp_file,
                                        out_file=fname_presuffix(
                                            sliced_apply_to_file,
                                            suffix='_qwarped'),
                                        environ=environ)
            warped_apply_to_slices.append(out_warp_apply.outputs.out_file)

    # Fix the obliquity
    oblique_warped_apply_to_slices = []
    for (sliced_apply_to_file, warped_apply_to_slice) in zip(
            sliced_apply_to_files, warped_apply_to_slices):
        oblique_slice = fix_obliquity(warped_apply_to_slice,
                                      sliced_apply_to_file,
                                      verbose=verbose,
                                      caching=caching,
                                      caching_dir=per_slice_dir,
                                      environ=environ)
        oblique_warped_apply_to_slices.append(oblique_slice)

    # Finally, merge all slices !
    out_merge_apply_to = merge(
        in_files=oblique_warped_apply_to_slices,
        dimension='z',
        merged_file=fname_presuffix(apply_to_file, suffix='_perslice',
                                    newpath=write_dir),
        environ=environ)

    # Fix the obliquity
    merged_apply_to_file = fix_obliquity(
        out_merge_apply_to.outputs.merged_file, apply_to_file,
        verbose=verbose, caching=caching,
        caching_dir=per_slice_dir, environ=environ)

    # Update the outputs
    output_files.extend(sliced_apply_to_files_to_remove +
                        oblique_warped_apply_to_slices)

    if not caching:
        for out_file in output_files:
            os.remove(out_file)

    return merged_apply_to_file


def _transform_to_template(to_register_filename, template_filename, write_dir,
                           func_to_anat_oned_filename,
                           anat_to_template_oned_filename,
                           anat_to_template_warp_filename,
                           voxel_size=None,
                           caching=False, verbose=True):
    """ Applies successive transforms to a given image to put it in
    template space.

    Parameters
    ----------
    to_register_filename : str
        Path to functional volume, coregistered to a common space with the
        anatomical volume.

    template_filename : str
        Template to register the functional to.

    func_to_anat_oned_filename : str
        Coregistration transform.

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
        warp_apply = memory.cache(afni.NwarpApply)
        resample = memory.cache(afni.Resample)
        warp_apply.interface().set_default_terminal_output(terminal_output)
        resample.interface().set_default_terminal_output(terminal_output)
    else:
        resample = afni.Resample(terminal_output=terminal_output).run
        warp_apply = afni.NwarpApply(terminal_output=terminal_output).run
        environ['AFNI_DECONFLICT'] = 'OVERWRITE'

    current_dir = os.getcwd()
    os.chdir(write_dir)
    normalized_filename = fname_presuffix(to_register_filename,
                                          suffix='_normalized')

    if voxel_size is None:
        resampled_template_filename = template_filename
    else:
        out_resample = resample(in_file=template_filename,
                                voxel_size=voxel_size,
                                outputtype='NIFTI_GZ')
        resampled_template_filename = out_resample.outputs.out_file

    transforms = [anat_to_template_warp_filename,
                  anat_to_template_oned_filename,
                  func_to_anat_oned_filename]
    warp = "'"
    warp += " ".join(transforms)
    warp += "'"
    _ = warp_apply(in_file=to_register_filename,
                   master=resampled_template_filename,
                   warp=warp,
                   out_file=normalized_filename)
    os.chdir(current_dir)
    return normalized_filename


def _apply_transforms(to_register_filename, target_filename,
                      write_dir,
                      transforms,
                      transformed_filename=None,
                      transforms_kind='nonlinear',
                      interpolation='wsinc5',
                      voxel_size=None, inverse=False,
                      caching=False, verbose=True):
    """ Applies successive transforms to a given image to put it in
    template space.

    Parameters
    ----------
    to_register_filename : str
        Path to functional volume, coregistered to a common space with the
        anatomical volume.

    template_filename : str
        Template to register the functional to.

    transforms : list
        List of transforms in order of 3dNWarpApply application: first must
        one must be in the target space and last one must be in
        the source space.

    inverse : bool, optional
        If True, after the transforms composition is computed, invert it.
        If the input transforms would take a dataset from space A to B,
        then the inverted transform will do the reverse.

    interpolation : one of {'nearestneighbour', 'linear', 'cubic', 'quintic',
                            'wsinc5'}, optional
        Interpolation type.

    voxel_size : 3-tuple of floats, optional
        Voxel size of the registered functional, in mm.

    caching : bool, optional
        Wether or not to use caching.

    verbose : bool, optional
        If True, all steps are verbose. Note that caching implies some
        verbosity in any case.
    """
    environ = {'AFNI_DECONFLICT': 'OVERWRITE'}
    if verbose:
        terminal_output = 'allatonce'
    else:
        terminal_output = 'none'

    if caching:
        memory = Memory(write_dir)
        catmatvec = memory.cache(afni.CatMatvec)
        allineate = memory.cache(afni.Allineate)
        warp_apply = memory.cache(afni.NwarpApply)
        resample = memory.cache(afni.Resample)        
        for step in [resample, allineate, warp_apply]:
            step.interface().set_default_terminal_output(terminal_output)
    else:
        resample = afni.Resample(terminal_output=terminal_output).run
        catmatvec = afni.CatMatvec().run
        allineate = afni.Allineate(terminal_output=terminal_output).run
        warp_apply = afni.NwarpApply(terminal_output=terminal_output).run

    if transformed_filename is None:
        target_basename = os.path.basename(target_filename)
        target_basename = os.path.splitext(target_basename)[0]
        target_basename = os.path.splitext(target_basename)[0]
        transformed_filename = fname_presuffix(
            to_register_filename, suffix='_to_' + target_basename,
            newpath=write_dir)

    if voxel_size is None:
        resampled_template_filename = target_filename
    else:
        out_resample = resample(in_file=target_filename,
                                voxel_size=voxel_size,
                                out_file=fname_presuffix(target_filename,
                                                         suffix='_resampled',
                                                         newpath=write_dir),
                                environ=environ)
        resampled_template_filename = out_resample.outputs.out_file
    if transforms_kind is not 'nonlinear':
        affine_transform_filename = fname_presuffix(transformed_filename,
                                                    suffix='.aff12.1D',
                                                    use_ext=False)
        out_catmatvec = catmatvec(in_file=[(transform, 'ONELINE')
                                           for transform in transforms],
                                  oneline=True,
                                  out_file=affine_transform_filename,
                                  environ=environ)
        if inverse:
            affine_transform_filename = fname_presuffix(transformed_filename,
                                                        suffix='_INV.aff12.1D',
                                                        use_ext=False)
            _ = catmatvec(in_file=[(out_catmatvec.outputs.out_file, 'I')],
                          oneline=True,
                          out_file=affine_transform_filename,
                          environ=environ)
        _ = allineate(
            in_file=to_register_filename,
            master=resampled_template_filename,
            in_matrix=affine_transform_filename,
            final_interpolation=interpolation,
            out_file=transformed_filename,
            environ=environ)
    else:
        warp = "'"
        warp += " ".join(transforms)
        warp += "'"
        _ = warp_apply(in_file=to_register_filename,
                       master=resampled_template_filename,
                       warp=warp,
                       inv_warp=inverse,
                       interp=interpolation,
                       out_file=transformed_filename,
                       environ=environ)
    # XXX obliquity information is lost if resampling is done
    transformed_filename = fix_obliquity(transformed_filename,
                                         resampled_template_filename)
    return transformed_filename