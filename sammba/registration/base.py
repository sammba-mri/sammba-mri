import os
import numpy as np
import nibabel
from sklearn.base import _pprint
from sklearn.utils.fixes import signature
from ..externals.nipype.caching import Memory
from ..externals.nipype.interfaces import afni, ants
from ..externals.nipype.utils.filemanip import fname_presuffix
from ..interfaces import segmentation
from .utils import fix_obliquity, _get_afni_output_type, compute_n4_max_shrink


def compute_brain_mask(head_file, brain_volume, write_dir, bias_correct=True,
                       caching=False,
                       terminal_output='allatonce',
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
    else:
        ComputeMask = segmentation.HistogramMask

    environ = {}
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
        environ['AFNI_DECONFLICT'] = 'OVERWRITE'

    if bias_correct:
        if unifize_kwargs is None:
            unifize_kwargs = {}
        out_unifize = unifize(
            in_file=head_file,
            out_file=fname_presuffix(head_file,
                                     suffix='_unifized_for_brain_mask',
                                     newpath=write_dir),
            environ=environ,
            **unifize_kwargs)
        file_to_mask = out_unifize.outputs.out_file
    else:
        file_to_mask = head_file

    out_clip_level = clip_level(in_file=file_to_mask)
    out_compute_mask = compute_mask(
        in_file=file_to_mask,
        out_file=fname_presuffix(file_to_mask,
                                 suffix='_brain_mask',
                                 newpath=write_dir),
        volume_threshold=brain_volume,
        intensity_threshold=int(out_clip_level.outputs.clip_val))

    if not caching and bias_correct:
        os.remove(out_unifize.outputs.out_file)

    return out_compute_mask.outputs.out_file


def _apply_mask(head_file, brain_mask_file, write_dir, bias_correct=True,
                caching=False, terminal_output='allatonce'):
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
    environ = {}
    if caching:
        memory = Memory(write_dir)
        calc = memory.cache(afni.Calc)
        calc.interface().set_default_terminal_output(terminal_output)
    else:
        calc = afni.Calc(terminal_output=terminal_output).run
        environ['AFNI_DECONFLICT'] = 'OVERWRITE'

    out_calc_mask = calc(in_file_a=head_file,
                         in_file_b=brain_mask_file,
                         expr='a*b',
                         out_file=fname_presuffix(head_file,
                                                  suffix='_brain',
                                                  newpath=write_dir))

    if not caching:
        os.remove(brain_mask_file)

    return out_calc_mask.outputs.out_file


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


def _bias_correct(in_file, write_dir, caching=False,
                  terminal_output='allatonce'):
    if caching:
        memory = Memory(write_dir)
        bias_correct = memory.cache(ants.N4BiasFieldCorrection)
        bias_correct.interface().set_default_terminal_output(terminal_output)
    else:
        bias_correct = ants.N4BiasFieldCorrection(
            terminal_output=terminal_output).run

    out_bias_correct = bias_correct(
        input_image=in_file,
        shrink_factor=compute_n4_max_shrink(in_file),
        output_image=fname_presuffix(in_file, suffix='_unbiased',
                                     newpath=write_dir))
    return out_bias_correct.outputs.output_image


def _afni_bias_correct(in_file, write_dir, caching=False,
                       terminal_output='allatonce',
                       **unifize_kwargs):
    environ = {}
    if caching:
        memory = Memory(write_dir)
        unifize = memory.cache(afni.Unifize)
        unifize.interface().set_default_terminal_output(terminal_output)
    else:
        unifize = afni.Unifize(terminal_output=terminal_output).run
        environ['AFNI_DECONFLICT'] = 'OVERWRITE'

    out_unifize = unifize(in_file=in_file,
                          out_file=fname_presuffix(in_file,
                                                   suffix='_unifized',
                                                   newpath=write_dir),
                          environ=environ,
                          **unifize_kwargs)

    return out_unifize.outputs.out_file


def _rigid_body_register(moving_head_file, reference_head_file,
                         moving_brain_mask, reference_brain_mask,
                         write_dir=None,
                         caching=False, terminal_output='allatonce'):
    # XXX: add verbosity
    if write_dir is None:
        write_dir = os.path.dirname(moving_head_file)

    environ = {}
    if caching:
        memory = Memory(write_dir)
        calc = memory.cache(afni.Calc)
        allineate = memory.cache(afni.Allineate)
        allineate2 = memory.cache(afni.Allineate)
        catmatvec = memory.cache(afni.CatMatvec)
        for step in [calc, allineate, allineate2]:
            step.interface().set_default_terminal_output(terminal_output)
    else:
        calc = afni.Calc(terminal_output=terminal_output).run
        allineate = afni.Allineate(terminal_output=terminal_output).run
        allineate2 = afni.Allineate(terminal_output=terminal_output).run  # TODO: remove after fixed bug
        catmatvec = afni.CatMatvec().run
        environ['AFNI_DECONFLICT'] = 'OVERWRITE'

    out_cacl = calc(in_file_a=moving_head_file,
                    in_file_b=moving_brain_mask,
                    expr='a*b',
                    outputtype=_get_afni_output_type(moving_head_file),
                    environ=environ)
    moving_brain_file = out_cacl.outputs.out_file

    out_cacl = calc(in_file_a=reference_head_file,
                    in_file_b=reference_brain_mask,
                    expr='a*b',
                    outputtype=_get_afni_output_type(reference_head_file),
                    environ=environ)
    reference_brain_file = out_cacl.outputs.out_file
    output_files = [moving_brain_file, reference_brain_file]

    # Compute the transformation from functional to anatomical brain
    # XXX: why in this sense
    out_allineate = allineate2(
        in_file=reference_brain_file,
        reference=moving_brain_file,
        out_matrix=fname_presuffix(reference_brain_file,
                                   suffix='_shr.aff12.1D',
                                   use_ext=False,
                                   new_path=write_dir),
        center_of_mass='',
        warp_type='shift_rotate',
        out_file=fname_presuffix(reference_brain_file, suffix='_shr'),
        environ=environ)
    rigid_transform_file = out_allineate.outputs.out_matrix
    output_files.append(out_allineate.outputs.out_file)

    # apply the inverse transform to register the anatomical to the func
    catmatvec_out_file = fname_presuffix(rigid_transform_file, suffix='INV',
                                         new_path=write_dir)
    if not os.path.isfile(catmatvec_out_file):
        _ = catmatvec(in_file=[(rigid_transform_file, 'I')],
                      oneline=True,
                      out_file=catmatvec_out_file)
        # XXX not cached I don't understand why
        output_files.append(catmatvec_out_file)
    out_allineate_apply = allineate(
        in_file=moving_head_file,
        master=reference_head_file,
        in_matrix=catmatvec_out_file,
        out_file=fname_presuffix(moving_head_file, suffix='_shr',
                                 new_path=write_dir),
        environ=environ)

    # Remove intermediate output
    if not caching:
        for output_file in output_files:
            os.remove(output_file)

    return out_allineate_apply.outputs.out_file, rigid_transform_file


def _warp(to_warp_file, reference_file, write_dir=None, caching=False,
          terminal_output='allatonce', verbose=True,
          overwrite=True, environ={}):
    # XXX: add verbosity and handle environ correctly
    environ = {}
    if write_dir is None:
        write_dir = os.path.dirname(to_warp_file)

    if caching:
        memory = Memory(write_dir)
        warp = memory.cache(afni.Warp)
    else:
        warp = afni.Warp().run
        environ['AFNI_DECONFLICT'] = 'OVERWRITE'

    # 3dWarp doesn't put the obliquity in the header, so do it manually
    # This step generates one file per slice and per time point, so we are
    # making sure they are removed at the end
    out_warp = warp(in_file=to_warp_file,
                    oblique_parent=reference_file,
                    interp='quintic',
                    gridset=reference_file,
                    out_file=fname_presuffix(to_warp_file, suffix='_warped',
                                             newpath=write_dir),
                    verbose=True,
                    environ=environ)
    warped_file = out_warp.outputs.out_file
    warped_oblique_file = fix_obliquity(
        warped_file, reference_file,
        overwrite=overwrite, verbose=verbose, caching=caching,
        caching_dir=write_dir)

    # Concatenate all the anat to func tranforms
    mat_file = fname_presuffix(warped_oblique_file, suffix='_warp.mat',
                               use_ext=False,
                               newpath=write_dir)
    if not os.path.isfile(mat_file):
        np.savetxt(mat_file, [out_warp.runtime.stdout], fmt='%s')
    return warped_oblique_file, mat_file


def _per_slice_qwarp(to_qwarp_file, reference_file,
                     voxel_size_x, voxel_size_y, apply_to_file=None,
                     write_dir=None,
                     caching=False, overwrite=True,
                     verbose=True, terminal_output='allatonce', environ={}):
    environ = {}
    if write_dir is None:
        write_dir = os.path.dirname(to_qwarp_file),

    if caching:
        memory = Memory(write_dir)
        resample = memory.cache(afni.Resample)
        slicer = memory.cache(afni.ZCutUp)
        warp_apply = memory.cache(afni.NwarpApply)
        qwarp = memory.cache(afni.Qwarp)
        merge = memory.cache(afni.Zcat)
        for step in [resample, slicer, warp_apply, qwarp, merge]:
            step.interface().set_default_terminal_output(terminal_output)
    else:
        resample = afni.Resample(terminal_output=terminal_output).run
        slicer = afni.ZCutUp(terminal_output=terminal_output).run
        warp_apply = afni.NwarpApply(terminal_output=terminal_output).run
        qwarp = afni.Qwarp(terminal_output=terminal_output).run
        merge = afni.Zcat(terminal_output=terminal_output).run
        environ['AFNI_DECONFLICT'] = 'OVERWRITE'

    # Slice anatomical image
    reference_img = nibabel.load(reference_file)
    reference_n_slices = reference_img.header.get_data_shape()[2]
    per_slice_dir = os.path.join(write_dir, 'per_slice')
    if not os.path.isdir(per_slice_dir):
        os.makedirs(per_slice_dir)
    sliced_reference_files = []
    for slice_n in range(reference_n_slices):
        out_slicer = slicer(in_file=reference_file,
                            keep='{0} {0}'.format(slice_n),
                            out_file=fname_presuffix(
                                reference_file, newpath=per_slice_dir,
                                suffix='_sl%d' % slice_n),
                            environ=environ)
        oblique_slice = fix_obliquity(out_slicer.outputs.out_file,
                          reference_file,
                          overwrite=overwrite, verbose=verbose,
                          caching=caching, caching_dir=per_slice_dir)
        sliced_reference_files.append(oblique_slice)

    # Slice mean functional
    sliced_to_qwarp_files = []
    img = nibabel.load(to_qwarp_file)
    n_slices = img.header.get_data_shape()[2]
    for slice_n in range(n_slices):
        out_slicer = slicer(in_file=to_qwarp_file,
                            keep='{0} {0}'.format(slice_n),
                            out_file=fname_presuffix(
                                to_qwarp_file, newpath=per_slice_dir,
                                suffix='_sl%d' % slice_n),
                            environ=environ)
        oblique_slice = fix_obliquity(out_slicer.outputs.out_file,
                          to_qwarp_file,
                          overwrite=overwrite, verbose=verbose,
                          caching=caching, caching_dir=per_slice_dir)
        sliced_to_qwarp_files.append(oblique_slice)

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
                environ=environ)
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
                                      overwrite=overwrite, verbose=verbose,
                                      caching=caching,
                                      caching_dir=per_slice_dir)
        oblique_resampled_warped_slices.append(oblique_slice)

    out_merge_func = merge(
        in_files=oblique_resampled_warped_slices,
        out_file=fname_presuffix(to_qwarp_file, suffix='_perslice',
                                 newpath=write_dir),
        environ=environ)

    # Fix the obliquity
    oblique_merged = fix_obliquity(out_merge_func.outputs.out_file,
                                   reference_file,
                                   overwrite=overwrite, verbose=verbose,
                                   caching=caching, caching_dir=per_slice_dir)

    # Collect the outputs
    output_files.extend(sliced_reference_files +
                        sliced_to_qwarp_files +
                        resampled_sliced_reference_files +
                        resampled_sliced_to_qwarp_files +
                        warped_slices + resampled_warped_slices)

    # Apply the precomputed warp slice by slice
    if apply_to_file is not None:
        # slice functional
        sliced_apply_to_files = []
        for slice_n in range(n_slices):
            out_slicer = slicer(in_file=apply_to_file,
                                keep='{0} {0}'.format(slice_n),
                                out_file=fname_presuffix(
                                    apply_to_file, newpath=per_slice_dir,
                                    suffix='_sl%d' % slice_n),
                                environ=environ)
            oblique_slice = fix_obliquity(out_slicer.outputs.out_file,
                                          apply_to_file,
                                          overwrite=overwrite, verbose=verbose,
                                          caching=caching,
                                          caching_dir=per_slice_dir)
            sliced_apply_to_files.append(oblique_slice)

        warped_apply_to_slices = []
        for (sliced_apply_to_file, warp_file) in zip(
                sliced_apply_to_files, warp_files):
            if warp_file is None:
                warped_apply_to_slices.append(sliced_apply_to_file)
            else:
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
                                          overwrite=overwrite, verbose=verbose,
                                          caching=caching,
                                          caching_dir=per_slice_dir)
            oblique_warped_apply_to_slices.append(oblique_slice)

        # Finally, merge all slices !
        out_merge_apply_to = merge(
            in_files=oblique_warped_apply_to_slices,
            out_file=fname_presuffix(apply_to_file, suffix='_perslice',
                                     newpath=write_dir),
            environ=environ)

        # Fix the obliquity
        merged_apply_to_file = fix_obliquity(
            out_merge_apply_to.outputs.out_file, apply_to_file,
            overwrite=overwrite, verbose=verbose, caching=caching,
            caching_dir=per_slice_dir)

        # Update the outputs
        output_files.extend(sliced_apply_to_files + warped_apply_to_slices)
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
                         caching=False, overwrite=True,
                         verbose=True, terminal_output='allatonce',
                         environ={}):

    # Apply the precomputed warp slice by slice

    environ = {}
    if write_dir is None:
        write_dir = os.path.dirname(apply_to_file),

    if caching:
        memory = Memory(write_dir)
        resample = memory.cache(afni.Resample)
        slicer = memory.cache(afni.ZCutUp)
        warp_apply = memory.cache(afni.NwarpApply)
        qwarp = memory.cache(afni.Qwarp)
        merge = memory.cache(afni.Zcat)
        for step in [resample, slicer, warp_apply, qwarp, merge]:
            step.interface().set_default_terminal_output(terminal_output)
    else:
        resample = afni.Resample(terminal_output=terminal_output).run
        slicer = afni.ZCutUp(terminal_output=terminal_output).run
        warp_apply = afni.NwarpApply(terminal_output=terminal_output).run
        qwarp = afni.Qwarp(terminal_output=terminal_output).run
        merge = afni.Zcat(terminal_output=terminal_output).run
        environ['AFNI_DECONFLICT'] = 'OVERWRITE'

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
    for slice_n in range(n_slices):
        out_slicer = slicer(in_file=apply_to_file,
                            keep='{0} {0}'.format(slice_n),
                            out_file=fname_presuffix(
                                apply_to_file, newpath=per_slice_dir,
                                suffix='_sl%d' % slice_n),
                            environ=environ)
        oblique_slice = fix_obliquity(out_slicer.outputs.out_file,
                                      apply_to_file,
                                      overwrite=overwrite, verbose=verbose,
                                      caching=caching,
                                      caching_dir=per_slice_dir)
        sliced_apply_to_files.append(oblique_slice)

    warped_apply_to_slices = []
    for (sliced_apply_to_file, warp_file) in zip(
            sliced_apply_to_files, warp_files):
        if warp_file is None:
            warped_apply_to_slices.append(sliced_apply_to_file)
        else:
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
                                      overwrite=overwrite, verbose=verbose,
                                      caching=caching,
                                      caching_dir=per_slice_dir)
        oblique_warped_apply_to_slices.append(oblique_slice)

    # Finally, merge all slices !
    out_merge_apply_to = merge(
        in_files=oblique_warped_apply_to_slices,
        out_file=fname_presuffix(apply_to_file, suffix='_perslice',
                                 newpath=write_dir),
        environ=environ)

    # Fix the obliquity
    merged_apply_to_file = fix_obliquity(
        out_merge_apply_to.outputs.out_file, apply_to_file,
        overwrite=overwrite, verbose=verbose, caching=caching,
        caching_dir=per_slice_dir)

    # Update the outputs
    output_files.extend(sliced_apply_to_files + warped_apply_to_slices)

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
    print('********************************', nibabel.load(to_register_filename).affine)
    _ = warp_apply(in_file=to_register_filename,
                   master=resampled_template_filename,
                   warp=warp,
                   out_file=normalized_filename)
    os.chdir(current_dir)
    return normalized_filename


class BaseRegistrator(object):
    """
    Base class for all registrators.
    """
    @classmethod
    def _get_param_names(cls):
        """Get parameter names for the registrato"""
        # fetch the constructor or the original constructor before
        # deprecation wrapping if any
        init = getattr(cls.__init__, 'deprecated_original', cls.__init__)
        if init is object.__init__:
            # No explicit constructor to introspect
            return []

        # introspect the constructor arguments to find the model parameters
        # to represent
        init_signature = signature(init)
        # Consider the constructor parameters excluding 'self'
        parameters = [p for p in init_signature.parameters.values()
                      if p.name != 'self' and p.kind != p.VAR_KEYWORD]
        for p in parameters:
            if p.kind == p.VAR_POSITIONAL:
                raise RuntimeError("registrators should always "
                                   "specify their parameters in the signature"
                                   " of their __init__ (no varargs)."
                                   " %s with constructor %s doesn't "
                                   " follow this convention."
                                   % (cls, init_signature))
        # Extract and sort argument names excluding 'self'
        return sorted([p.name for p in parameters])

    def _get_params(self):
        """Get parameters of the session.

        Returns
        -------
        params : mapping of string to any
            Parameter names mapped to their values.
        """
        out = dict()
        for key in self._get_param_names():
            value = getattr(self, key, None)
            out[key] = value
        return out

    def _set_params(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def __repr__(self):
        class_name = self.__class__.__name__
        return '%s(%s)' % (class_name, _pprint(self._get_params(),
                           offset=len(class_name),),)

    def _set_output_dir(self):
        if self.output_dir is None:
            self.output_dir = os.getcwd()

        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)
