import os
import fnmatch
import glob
import nibabel
from nipype.caching import Memory
from nipype.interfaces import afni, fsl
from nipype.utils.filemanip import fname_presuffix
from nipype.interfaces.fsl.base import Info
from ..orientation import fix_obliquity


def _delete_orientation(in_file, write_dir=None, min_zoom=.1, caching=False,
                        verbose=True):
    if write_dir is None:
        write_dir = os.path.dirname(in_file)

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
        verbose=verbose)
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
                    save_warp=True,
                    environ=environ)

    # 3dWarp doesn't put the obliquity in the header, so do it manually
    warped_oblique_file = fix_obliquity(
        out_warp.outputs.out_file, reference_file,
        verbose=verbose, caching=caching,
        caching_dir=write_dir, environ=environ)

    return warped_oblique_file, out_warp.outputs.warp_file


def _get_fsl_slice_output_files(out_base_name, output_type):
    ext = Info.output_type_to_ext(output_type)
    suffix = '_slice_*' + ext
    exact_pattern = '_slice_[0-9][0-9][0-9][0-9]' + ext
    fname_template = os.path.abspath(
        out_base_name + suffix)
    fname_exact_pattern = os.path.abspath(
        out_base_name + exact_pattern)
    sliced_files = fnmatch.filter(sorted(glob.glob(fname_template)),
                                  fname_exact_pattern)
    return sliced_files


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
    # XXX: workaround for nipype globbing to find slicer outputs
    # Use out_slicer.outputs.out_files once fixed
    sliced_reference_files = _get_fsl_slice_output_files(
        out_slicer.inputs['out_base_name'], out_slicer.inputs['output_type'])

    # Slice mean functional
    out_slicer = slicer(in_file=to_qwarp_file,
                        out_base_name=fname_presuffix(to_qwarp_file,
                                                      newpath=per_slice_dir,
                                                      use_ext=False))
    sliced_to_qwarp_files = _get_fsl_slice_output_files(
        out_slicer.inputs['out_base_name'], out_slicer.inputs['output_type'])

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
        sliced_apply_to_files = _get_fsl_slice_output_files(
                out_slicer.inputs['out_base_name'],
                out_slicer.inputs['output_type'])
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
    sliced_apply_to_files = _get_fsl_slice_output_files(
        out_slicer.inputs['out_base_name'], out_slicer.inputs['output_type'])

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
                      interpolation=None,
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
                            'wsinc5'} or None, optional
        Interpolation type. If None, AFNI defaults are used.

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
        if interpolation is None:
            _ = allineate(
                in_file=to_register_filename,
                master=resampled_template_filename,
                in_matrix=affine_transform_filename,
                out_file=transformed_filename,
                environ=environ)
        else:
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
        if interpolation is None:
            _ = warp_apply(in_file=to_register_filename,
                           master=resampled_template_filename,
                           warp=warp,
                           inv_warp=inverse,
                           out_file=transformed_filename,
                           environ=environ)
        else:
            _ = warp_apply(in_file=to_register_filename,
                           master=resampled_template_filename,
                           warp=warp,
                           inv_warp=inverse,
                           interp=interpolation,
                           out_file=transformed_filename,
                           environ=environ)

    # XXX obliquity information is lost if resampling is done
    transformed_filename = fix_obliquity(transformed_filename,
                                         resampled_template_filename,
                                         verbose=verbose, caching=caching,
                                         caching_dir=write_dir,
                                         environ=environ)
    return transformed_filename