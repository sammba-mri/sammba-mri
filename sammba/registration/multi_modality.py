import os
import numpy as np
import nibabel
from ..externals.nipype.caching import Memory
from ..externals.nipype.interfaces import afni, fsl
from ..externals.nipype.utils.filemanip import fname_presuffix
from ..interfaces import RatsMM
from .utils import fix_obliquity


def _extract_brain(head_filename, write_dir, brain_volume, caching=False,
                   terminal_output='allatonce', environ={}):
    if caching:
        memory = Memory(write_dir)
        clip_level = memory.cache(afni.ClipLevel)
        rats = memory.cache(RatsMM)
        calc = memory.cache(afni.Calc)
        rats.interface().set_default_terminal_output(terminal_output)
        calc.interface().set_default_terminal_output(terminal_output)
    else:
        clip_level = afni.ClipLevel().run
        rats = RatsMM(terminal_output=terminal_output).run
        calc = afni.Calc(terminal_output=terminal_output).run

    out_clip_level = clip_level(in_file=head_filename)
    out_rats_func = rats(
        in_file=head_filename,
        volume_threshold=brain_volume,
        intensity_threshold=int(out_clip_level.outputs.clip_val))
    out_cacl = calc(in_file_a=head_filename,
                    in_file_b=out_rats_func.outputs.out_file,
                    expr='a*b',
                    out_file=fname_presuffix(head_filename, suffix='_brain'),
                    environ=environ)

    if not caching:
        os.remove(out_rats_func.outputs.out_file)

    return out_cacl.outputs.out_file


def _rigid_body_register(to_register_filename, to_register_brain_filename,
                         target_filename, target_brain_filename,
                         write_dir, brain_volume, caching=False,
                         terminal_output='allatonce',
                         environ={}):
    if caching:
        memory = Memory(write_dir)
        allineate = memory.cache(afni.Allineate)
        allineate2 = memory.cache(afni.Allineate)
        catmatvec = memory.cache(afni.CatMatvec)
        for step in [allineate, allineate2]:
            step.interface().set_default_terminal_output(terminal_output)
    else:
        allineate = afni.Allineate(terminal_output=terminal_output).run
        allineate2 = afni.Allineate(terminal_output=terminal_output).run  # TODO: remove after fixed bug
        catmatvec = afni.CatMatvec().run

    # Compute the transformation from functional to anatomical brain
    # XXX: why in this sense
    out_allineate = allineate2(
        in_file=target_brain_filename,
        reference=to_register_brain_filename,
        out_matrix=fname_presuffix(target_brain_filename,
                                   suffix='_shiftrot.aff12.1D',
                                   use_ext=False),
        center_of_mass='',
        warp_type='shift_rotate',
        out_file=fname_presuffix(target_brain_filename,
                                 suffix='_shiftrot'),
        environ=environ)
    rigid_transform_file = out_allineate.outputs.out_matrix
    output_files = [out_allineate.outputs.out_file]

    # apply the inverse transform to register the anatomical to the func
    catmatvec_out_file = fname_presuffix(rigid_transform_file,
                                         suffix='INV')
    if not os.path.isfile(catmatvec_out_file):
        _ = catmatvec(in_file=[(rigid_transform_file, 'I')],
                      oneline=True,
                      out_file=catmatvec_out_file)
        # XXX not cached I don't understand why
        output_files.append(catmatvec_out_file)
    out_allineate_apply = allineate(
        in_file=to_register_filename,
        master=target_filename,
        in_matrix=catmatvec_out_file,
        out_file=fname_presuffix(to_register_filename,
                                 suffix='_shiftrot'),
        environ=environ)

    # Remove intermediate output
    if not caching:
        for output_file in output_files:
            os.remove(output_file)

    return out_allineate_apply.outputs.out_file, rigid_transform_file


def _warp(to_warp_filename, target_filename, write_dir, caching=False,
          terminal_output='allatonce', environ={}):
    if caching:
        memory = Memory(write_dir)
        warp = memory.cache(afni.Warp)
        warp.interface().set_default_terminal_output(terminal_output)
    else:
        warp = afni.Warp().run

    # 3dWarp doesn't put the obliquity in the header, so do it manually
    # This step generates one file per slice and per time point, so we are
    # making sure they are removed at the end
    out_warp = warp(in_file=to_warp_filename,
                    oblique_parent=target_filename,
                    interp='quintic',
                    gridset=target_filename,
                    outputtype='NIFTI_GZ',
                    verbose=True,
                    environ=environ)
    warped_filename = out_warp.outputs.out_file
    warped_oblique_filename = fix_obliquity(
        warped_filename, target_filename,
        overwrite=not caching, verbose=(terminal_output != 'none'))

    # Concatenate all the anat to func tranforms
    mat_filename = fname_presuffix(warped_filename,
                                   suffix='_warp.mat', use_ext=False)
    output_files = []
    if not os.path.isfile(mat_filename):
        np.savetxt(mat_filename, [out_warp.runtime.stdout], fmt='%s')
        output_files.append(mat_filename)
    return warped_oblique_filename, mat_filename, output_files


def _per_slice_qwarp(to_qwarp_filename, target_filename, write_dir,
                     voxel_size_x, voxel_size_y, apply_to_filename=None,
                     caching=False, terminal_output='allatonce', environ={}):
    if caching:
        memory = Memory(write_dir)
        resample = memory.cache(afni.Resample)
        slicer = memory.cache(afni.ZCutUp)
        warp_apply = memory.cache(afni.NwarpApply)
        qwarp = memory.cache(afni.Qwarp)
        merge = memory.cache(fsl.Merge)
        for step in [resample, slicer, warp_apply, qwarp, merge]:
            step.interface().set_default_terminal_output(terminal_output)
    else:
        resample = afni.Resample(terminal_output=terminal_output).run
        slicer = afni.ZCutUp(terminal_output=terminal_output).run
        warp_apply = afni.NwarpApply(terminal_output=terminal_output).run
        qwarp = afni.Qwarp(terminal_output=terminal_output).run
        merge = fsl.Merge(terminal_output=terminal_output).run

    overwrite = not caching,
    verbose = (terminal_output != 'none')

    # Slice anatomical image
    target_img = nibabel.load(target_filename)
    target_n_slices = target_img.header.get_data_shape()[2]
    sliced_target_filenames = []
    for slice_n in range(target_n_slices):
        out_slicer = slicer(in_file=target_filename,
                            keep='{0} {0}'.format(slice_n),
                            out_file=fname_presuffix(
                                target_filename,
                                suffix='_slice%d' % slice_n),
                            environ=environ)
        _ = fix_obliquity(out_slicer.outputs.out_file,
                          target_filename,
                          overwrite=overwrite, verbose=verbose)
        sliced_target_filenames.append(out_slicer.outputs.out_file)

    # Slice mean functional
    sliced_to_qwarp_filenames = []
    img = nibabel.load(to_qwarp_filename)
    n_slices = img.header.get_data_shape()[2]
    for slice_n in range(n_slices):
        out_slicer = slicer(in_file=to_qwarp_filename,
                            keep='{0} {0}'.format(slice_n),
                            out_file=fname_presuffix(
                                to_qwarp_filename,
                                suffix='_slice%d' % slice_n),
                            environ=environ)
        _ = fix_obliquity(out_slicer.outputs.out_file,
                          to_qwarp_filename,
                          overwrite=overwrite, verbose=verbose)
        sliced_to_qwarp_filenames.append(out_slicer.outputs.out_file)

    # Below line is to deal with slices where there is no signal (for example
    # rostral end of some anatomicals)

    # The inverse warp frequently fails, Resampling can help it work better
    # XXX why specifically .1 in voxel_size ?
    voxel_size_z = target_img.header.get_zooms()[2]
    resampled_sliced_target_filenames = []
    for sliced_target_filename in sliced_target_filenames:
        out_resample = resample(in_file=sliced_target_filename,
                                voxel_size=(voxel_size_x, voxel_size_y,
                                            voxel_size_z),
                                outputtype='NIFTI_GZ',
                                environ=environ)
        resampled_sliced_target_filenames.append(out_resample.outputs.out_file)

    resampled_sliced_to_qwarp_filenames = []
    for sliced_to_qwarp_filename in sliced_to_qwarp_filenames:
        out_resample = resample(in_file=sliced_to_qwarp_filename,
                                voxel_size=(voxel_size_x, voxel_size_y,
                                            voxel_size_z),
                                outputtype='NIFTI_GZ',
                                environ=environ)
        resampled_sliced_to_qwarp_filenames.append(
            out_resample.outputs.out_file)

    # single slice non-linear functional to anatomical registration
    warped_slices = []
    warp_filenames = []
    output_files = []
    for (resampled_sliced_to_qwarp_filename,
         resampled_sliced_target_filename) in zip(
            resampled_sliced_to_qwarp_filenames,
            resampled_sliced_target_filenames):
        warped_slice = fname_presuffix(resampled_sliced_to_qwarp_filename,
                                       suffix='_qwarp')
        out_qwarp = qwarp(in_file=resampled_sliced_to_qwarp_filename,
                          base_file=resampled_sliced_target_filename,
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
            fname_presuffix(out_qwarp.outputs.warped_source,
                            suffix='_Allin.nii', use_ext=False),
            fname_presuffix(out_qwarp.outputs.warped_source,
                            suffix='_Allin.aff12.1D', use_ext=False)])

    # Resample the mean volume back to the initial resolution,
    voxel_size = nibabel.load(to_qwarp_filename).header.get_zooms()
    resampled_warped_slices = []
    for warped_slice in warped_slices:
        out_resample = resample(in_file=warped_slice,
                                voxel_size=voxel_size,
                                outputtype='NIFTI_GZ',
                                environ=environ)
        resampled_warped_slices.append(out_resample.outputs.out_file)

    # fix the obliquity
    for (sliced_target_filename, resampled_warped_slice) in zip(
            sliced_target_filenames, resampled_warped_slices):
        _ = fix_obliquity(resampled_warped_slice,
                          sliced_target_filename,
                          overwrite=overwrite, verbose=verbose)

    out_merge_func = merge(
        in_files=resampled_warped_slices, dimension='z',
        merged_file=resampled_warped_slices[0].replace('_slice0',
                                                       '_perslice'))

    # slice functional
    sliced_apply_to_filenames = []
    for slice_n in range(n_slices):
        out_slicer = slicer(in_file=apply_to_filename,
                            keep='{0} {0}'.format(slice_n),
                            out_file=fname_presuffix(
                                apply_to_filename,
                                suffix='_slice%d' % slice_n),
                            environ=environ)
        _ = fix_obliquity(out_slicer.outputs.out_file,
                          apply_to_filename,
                          overwrite=overwrite, verbose=verbose)
        sliced_apply_to_filenames.append(out_slicer.outputs.out_file)

    # Collect the outputs
    output_files.extend(sliced_target_filenames +
                        sliced_to_qwarp_filenames +
                        resampled_sliced_target_filenames +
                        resampled_sliced_to_qwarp_filenames +
                        warped_slices + resampled_warped_slices)

    # Apply the precomputed warp slice by slice
    if apply_to_filename is not None:
        warped_apply_to_slices = []
        for (sliced_apply_to_filename, warp_filename) in zip(
                sliced_apply_to_filenames, warp_filenames):
            out_warp_apply = warp_apply(in_file=sliced_apply_to_filename,
                                        master=sliced_apply_to_filename,
                                        warp=warp_filename,
                                        out_file=fname_presuffix(
                                            sliced_apply_to_filename,
                                            suffix='_qwarp'),
                                        environ=environ)
            warped_apply_to_slices.append(out_warp_apply.outputs.out_file)

        # Fix the obliquity
        for (sliced_apply_to_filename, warped_apply_to_slice) in zip(
                sliced_apply_to_filenames, warped_apply_to_slices):
            _ = fix_obliquity(warped_apply_to_slice, sliced_apply_to_filename,
                              overwrite=overwrite, verbose=verbose)

        # Finally, merge all slices !
        out_merge_apply_to = merge(
            in_files=warped_apply_to_slices,
            dimension='z',
            merged_file=warped_apply_to_slices[0].replace('_slice0',
                                                          '_perslice'))

        # Update the outputs
        output_files.extend(sliced_apply_to_filenames + warped_apply_to_slices)

        merged_apply_to_filename = out_merge_apply_to.outputs.merged_file
    else:
        merged_apply_to_filename = None

    if not caching:
        for out_file in output_files:
            os.remove(out_file)

    return (out_merge_func.outputs.merged_file, warp_filenames,
            merged_apply_to_filename)
