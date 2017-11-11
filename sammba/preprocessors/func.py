import os
import numpy as np
import nibabel
from nipype.caching import Memory
from nipype.interfaces import afni, fsl, ants
from nipype.utils.filemanip import fname_presuffix
from sammba.interfaces import RatsMM
from .utils import fix_obliquity


def register_func_to_anat(func_filename, anat_filename, tr, write_dir,
                          caching=False):
    """
    The functional volume is aligned to the anatomical, first with a rigid body
    registration and then on a per-slice basis (only a fine correction, this is
    mostly for correction of EPI distortion). This pipeline includes
    slice timing.
    """
    if caching:
        memory = Memory(write_dir, base_dir='sammba_cache')
        tshift = memory.cache(afni.TShift)
        clip_level = memory.cache(afni.ClipLevel)
        threshold = memory.cache(fsl.Threshold)
        volreg = memory.cache(afni.Volreg)
        allineate = memory.cache(afni.Allineate)
        copy_geom = memory.cache(fsl.CopyGeom)
        bias_correct = memory.cache(ants.N4BiasFieldCorrection)
        tstat = memory.cache(afni.TStat)
        rats = memory.cache(RatsMM)
        calc = memory.cache(afni.Calc)
        allineate = memory.cache(afni.Allineate)
        allineate2 = memory.cache(afni.Allineate)
        unifize = memory.cache(afni.Unifize)
        catmatvec = memory.cache(afni.CatMatvec)
        warp = memory.cache(afni.Warp)
        resample = memory.cache(afni.Resample)
        slicer = memory.cache(afni.ZCutUp)
        warp_apply = memory.cache(afni.NwarpApply)
        qwarp = memory.cache(afni.Qwarp)
        merge = memory.cache(fsl.Merge)
    else:
        tshift = afni.TShift().run
        clip_level = afni.ClipLevel().run
        threshold = fsl.Threshold().run
        volreg = afni.Volreg().run
        allineate = afni.Allineate().run
        allineate2 = afni.Allineate().run  # TODO: remove after fixed bug
        copy_geom = fsl.CopyGeom().run
        bias_correct = ants.N4BiasFieldCorrection().run
        tstat = afni.TStat().run
        rats = RatsMM().run
        calc = afni.Calc().run
        allineate = afni.Allineate().run
        unifize = afni.Unifize().run
        catmatvec = afni.CatMatvec().run
        warp = afni.Warp().run
        resample = afni.Resample().run
        slicer = afni.ZCutUp().run
        warp_apply = afni.NwarpApply().run
        qwarp = afni.Qwarp().run
        merge = fsl.Merge().run

    # Correct slice timing
    os.chdir(write_dir)
    out_tshift = tshift(in_file=func_filename,
                        outputtype='NIFTI_GZ',
                        tpattern='altplus',
                        tr=tr)
    tshifted_filename = out_tshift.outputs.out_file

    # Register to the first volume
    # XXX why do you need a thresholded image ?    
    out_clip_level = clip_level(in_file=tshifted_filename)
    out_threshold = threshold(in_file=tshifted_filename,
                              thresh=out_clip_level.outputs.clip_val)
    thresholded_filename = out_threshold.outputs.out_file

    oned_filename = fname_presuffix(thresholded_filename,
                                    suffix='Vr.1Dfile.1D',
                                    use_ext=False)
    oned_matrix_filename = fname_presuffix(thresholded_filename,
                                           suffix='Vr.aff12.1D',
                                           use_ext=False)
    out_volreg = volreg(in_file=thresholded_filename,
                        outputtype='NIFTI_GZ',
                        oned_file=oned_filename,  # XXX dfile not saved
                        oned_matrix_save=oned_matrix_filename)
    # XXX: bad output: up and down on y-axis

    # Apply the registration to the whole head
    allineated_filename = fname_presuffix(tshifted_filename, suffix='Av')
    out_allineate = allineate(in_file=tshifted_filename,
                              master=tshifted_filename,
                              in_matrix=out_volreg.outputs.oned_matrix_save,
                              out_file=allineated_filename)

    # 3dAllineate removes the obliquity. This is not a good way to readd it as
    # removes motion correction info in the header if it were an AFNI file...as
    # it happens it's NIfTI which does not store that so irrelevant!
    out_copy_geom = copy_geom(dest_file=out_allineate.outputs.out_file,
                              in_file=out_volreg.outputs.out_file)
    allineated_filename = out_copy_geom.outputs.out_file

    # XXX: bad output: up and down on y-axis

    # Create a (hopefully) nice mean image for use in the registration
    out_tstat = tstat(in_file=allineated_filename, args='-mean',
                      outputtype='NIFTI_GZ')
    averaged_filename = out_tstat.outputs.out_file

    # Correct the functional average for intensities bias
    out_bias_correct = bias_correct(input_image=averaged_filename, dimension=3)
    unbiased_func_filename = out_bias_correct.outputs.output_image

    # Bias correct the antomical image
    out_unifize = unifize(in_file=anat_filename, outputtype='NIFTI_GZ')
    unbiased_anat_filename = out_unifize.outputs.out_file

    # Mask the mean functional volume outside the brain.
    out_clip_level = clip_level(in_file=unbiased_func_filename)
    # XXX bad: brain mask cut
    out_rats = rats(in_file=unbiased_func_filename,
                    out_file='UnBm.nii.gz',
                    volume_threshold=400,
                    intensity_threshold=int(out_clip_level.outputs.clip_val))
    out_cacl = calc(in_file_a=unbiased_func_filename,
                    in_file_b=out_rats.outputs.out_file,
                    expr='a*b',
                    outputtype='NIFTI_GZ')

    # Compute the transformation from the functional image to the anatomical
    # XXX: why in this sense
    out_allineate = allineate2(
        in_file=out_cacl.outputs.out_file,
        reference=unbiased_anat_filename,
        out_matrix='BmBe_shr.aff12.1D',
        center_of_mass='',
        warp_type='shift_rotate',
        out_file='BmBe_shr.nii.gz')
    rigid_transform_file = out_allineate.outputs.out_matrix

    # apply the inverse transformation to register to the anatomical volume
    catmatvec_out_file = fname_presuffix(rigid_transform_file, suffix='INV')
    if not os.path.isfile(catmatvec_out_file):
        _ = catmatvec(in_file=[(rigid_transform_file, 'I')],
                      oneline=True,
                      out_file=catmatvec_out_file)
        # XXX not cached I don't understand why
    out_allineate = allineate(
        in_file=unbiased_anat_filename,
        master=unbiased_func_filename,
        in_matrix=catmatvec_out_file,
        out_file='BmBe_shr_in_func_space.nii.gz')

    # suppanatwarp="$base"_BmBe_shr.aff12.1D

    # Non-linear registration
    # XXX what is the difference between Warp and 3dQwarp?    
    out_warp = warp(in_file=out_allineate.outputs.out_file,
                    oblique_parent=unbiased_func_filename,
                    interp='quintic',
                    gridset=unbiased_func_filename,
                    outputtype='NIFTI_GZ',
                    verbose=True)
    registered_anat_filename = out_warp.outputs.out_file
    mat_filename = fname_presuffix(registered_anat_filename,
                                   suffix='_warp.mat', use_ext=False)
    if not os.path.isfile(mat_filename):
        np.savetxt(mat_filename, [out_warp.runtime.stdout], fmt='%s')

    # 3dWarp doesn't put the obliquity in the header, so do it manually
    # This step generates one file per slice and per time point, so we are
    # making sure they are removed at the end
    registered_anat_oblique_filename = fix_obliquity(
        registered_anat_filename, unbiased_func_filename,
        overwrite=False)

    # Slice anatomical image
    anat_img = nibabel.load(registered_anat_oblique_filename)
    anat_n_slices = anat_img.header.get_data_shape()[2]
    sliced_registered_anat_filenames = []
    for slice_n in range(anat_n_slices):
        out_slicer = slicer(in_file=registered_anat_oblique_filename,
                            keep='{0} {1}'.format(slice_n, slice_n),
                            out_file=fname_presuffix(
                                registered_anat_oblique_filename,
                                suffix='Sl%d' % slice_n))
        sliced_registered_anat_filenames.append(out_slicer.outputs.out_file)

    # Slice mean functional
    sliced_bias_corrected_filenames = []
    img = nibabel.load(func_filename)
    n_slices = img.header.get_data_shape()[2]
    for slice_n in range(n_slices):
        out_slicer = slicer(in_file=unbiased_func_filename,
                            keep='{0} {1}'.format(slice_n, slice_n),
                            out_file=fname_presuffix(unbiased_func_filename,
                                                     suffix='Sl%d' % slice_n))
        sliced_bias_corrected_filenames.append(out_slicer.outputs.out_file)

    # Below line is to deal with slices where there is no signal (for example
    # rostral end of some anatomicals)

    # The inverse warp frequently fails, Resampling can help it work better
    voxel_size_z = anat_img.header.get_zooms()[2]
    resampled_registered_anat_filenames = []
    for sliced_registered_anat_filename in sliced_registered_anat_filenames:
        out_resample = resample(in_file=sliced_registered_anat_filename,
                                voxel_size=(.1, .1, voxel_size_z),
                                outputtype='NIFTI_GZ')
        resampled_registered_anat_filenames.append(
            out_resample.outputs.out_file)

    resampled_bias_corrected_filenames = []
    for sliced_bias_corrected_filename in sliced_bias_corrected_filenames:
        out_resample = resample(in_file=sliced_bias_corrected_filename,
                                voxel_size=(.1, .1, voxel_size_z),
                                outputtype='NIFTI_GZ')
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
                          iwarp=True,
                          noneg=True,
                          blur=[0],
                          nmi=True,
                          noXdis=True,
                          allineate=True,
                          allineate_opts='-parfix 1 0 -parfix 2 0 -parfix 3 0 '
                                         '-parfix 4 0 -parfix 5 0 -parfix 6 0 '
                                         '-parfix 7 0 -parfix 9 0 '
                                         '-parfix 10 0 -parfix 12 0',
                          out_file=warped_slice)
        warped_slices.append(out_qwarp.outputs.warped_source)
        warp_filenames.append(out_qwarp.outputs.source_warp)

    # Resample the mean volume back to the initial resolution,
    voxel_size = nibabel.load(
        sliced_bias_corrected_filename).header.get_zooms()
    resampled_warped_slices = []
    for warped_slice in warped_slices:
        out_resample = resample(in_file=warped_slice,
                                voxel_size=voxel_size + (voxel_size[0],),
                                outputtype='NIFTI_GZ')
        resampled_warped_slices.append(out_resample.outputs.out_file)

    # fix the obliquity
    for (resampled_registered_anat_filename, resampled_warped_slice) in zip(
            resampled_registered_anat_filenames, resampled_warped_slices):
        _ = fix_obliquity(resampled_warped_slice,
                          resampled_registered_anat_filename)

    # Merge all slices !
#    out_merge_mean = merge(in_files=resampled_warped_slices, dimension='z')

    # slice functional
    sliced_func_filenames = []
    for slice_n in range(n_slices):
        out_slicer = slicer(in_file=allineated_filename,
                            keep='{0} {1}'.format(slice_n, slice_n),
                            out_file=fname_presuffix(allineated_filename,
                                                     suffix='Sl%d' % slice_n))
        sliced_func_filenames.append(out_slicer.outputs.out_file)

    # resample functional slices
    resampled_func_filenames = []
    for sliced_func_filename in sliced_func_filenames:
        out_resample = resample(in_file=sliced_func_filename,
                                voxel_size=(.1, .1, voxel_size_z),
                                outputtype='NIFTI_GZ')
        resampled_func_filenames.append(out_resample.outputs.out_file)

    # Apply the precomputed warp slice by slice
    warped_func_slices = []
    for (resampled_func_filename, warp_filename) in zip(
            resampled_func_filenames, warp_filenames):
        out_warp_apply = warp_apply(in_file=resampled_func_filename,
                                    warp=warp_filename,
                                    out_file=fname_presuffix(
                                        resampled_func_filename, suffix='_qw'))
        warped_func_slices.append(out_warp_apply.outputs.out_file)

    # Fix the obliquity
    # XXX why no resampling back before ?
    for (resampled_registered_anat_filename, warped_func_slice) in zip(
            resampled_registered_anat_filenames, warped_func_slices):
        _ = fix_obliquity(warped_func_slice,
                          resampled_registered_anat_filename)

    # Finally, merge all slices !
    out_merge_func = merge(in_files=warped_func_slices, dimension='z')
    out_merge_anat = merge(in_files=resampled_registered_anat_filenames,
                           dimension='z')
    return (out_merge_func.outputs.merged_file,
            out_merge_anat.outputs.merged_file)
