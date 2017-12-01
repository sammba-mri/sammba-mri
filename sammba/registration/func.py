import os
import numpy as np
import nibabel
from sammba.externals.nipype.caching import Memory
from sammba.externals.nipype.interfaces import afni, fsl, ants
from sammba.externals.nipype.utils.filemanip import fname_presuffix
from sammba.interfaces import RatsMM
from .utils import fix_obliquity
from .fmri_session import FMRISession
from .t1 import anats_to_template


def coregister_fmri_session(session_data, t_r, write_dir,
                            slice_timing=True,
                            prior_rigid_body_registration=False,
                            caching=False):
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

    caching : bool, optional
        Wether or not to use caching.

    Returns
    -------
    the same sequence with each animal_data updated: attributes
        `registered_func_` and
        `registered_anat_` are added to specify the paths to the functional
        and anatomical images registered to the template,
        and `output_dir_` is added to give output directory for each animal.
    """
    func_filename = session_data.func
    anat_filename = session_data.anat
    if caching:
        memory = Memory(write_dir)
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

    session_data._check_inputs()
    output_dir = os.path.join(os.path.abspath(write_dir),
                              session_data.animal_id)
    session_data._set_output_dir_(output_dir)
    current_dir = os.getcwd()
    os.chdir(output_dir)
    # Correct slice timing
    if slice_timing:
        out_tshift = tshift(in_file=func_filename,
                            outputtype='NIFTI_GZ',
                            tpattern='altplus',
                            tr=str(t_r))
        func_filename = out_tshift.outputs.out_file

    # Register to the first volume
    # XXX why do you need a thresholded image ?
    out_clip_level = clip_level(in_file=func_filename)
    out_threshold = threshold(in_file=func_filename,
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
    allineated_filename = fname_presuffix(func_filename, suffix='Av')
    out_allineate = allineate(in_file=func_filename,
                              master=func_filename,
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

    if prior_rigid_body_registration:
        # Mask the mean functional volume outside the brain.
        out_clip_level = clip_level(in_file=unbiased_func_filename)
        # XXX bad: brain mask cut
        out_rats = rats(in_file=unbiased_func_filename,
                        volume_threshold=400,
                        intensity_threshold=int(out_clip_level.outputs.clip_val))
        out_cacl_func = calc(in_file_a=unbiased_func_filename,
                             in_file_b=out_rats.outputs.out_file,
                             expr='a*b',
                             outputtype='NIFTI_GZ')
    
        # Mask the anatomical volume outside the brain.
        out_clip_level = clip_level(in_file=unbiased_anat_filename)
        # XXX bad: brain mask cut
        out_rats = rats(in_file=unbiased_anat_filename,
                        volume_threshold=400,
                        intensity_threshold=int(out_clip_level.outputs.clip_val))
        out_cacl_anat = calc(in_file_a=unbiased_anat_filename,
                             in_file_b=out_rats.outputs.out_file,
                             expr='a*b',
                             outputtype='NIFTI_GZ')
    
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
                                     suffix='_shr'))
        rigid_transform_file = out_allineate.outputs.out_matrix
    
        # apply the inverse transformation to register the anatomical to the func
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
            out_file=fname_presuffix(unbiased_anat_filename,
                                     suffix='_shr_in_func_space'))
        allineated_anat_filename = out_allineate.outputs.out_file
    else:
        allineated_anat_filename = unbiased_anat_filename

    # suppanatwarp="$base"_BmBe_shr.aff12.1D

    # Non-linear registration
    # XXX what is the difference between Warp and 3dQwarp?
    out_warp = warp(in_file=allineated_anat_filename,
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

    # func to anat tranform
    transform_filename = fname_presuffix(registered_anat_filename,
                                         suffix='_func_to_anat.aff12.1D',
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
        _ = fix_obliquity(out_slicer.outputs.out_file,
                          registered_anat_oblique_filename)
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
        _ = fix_obliquity(out_slicer.outputs.out_file,
                          unbiased_func_filename)
        sliced_bias_corrected_filenames.append(out_slicer.outputs.out_file)

    # Below line is to deal with slices where there is no signal (for example
    # rostral end of some anatomicals)

    # The inverse warp frequently fails, Resampling can help it work better
    # XXX why specifically .1 in voxel_size ?
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
        unbiased_func_filename).header.get_zooms()
    resampled_warped_slices = []
    for warped_slice in warped_slices:
        out_resample = resample(in_file=warped_slice,
                                voxel_size=voxel_size,
                                outputtype='NIFTI_GZ')
        resampled_warped_slices.append(out_resample.outputs.out_file)

    # fix the obliquity
    for (sliced_registered_anat_filename, resampled_warped_slice) in zip(
            sliced_registered_anat_filenames, resampled_warped_slices):
        _ = fix_obliquity(resampled_warped_slice,
                          sliced_registered_anat_filename)

    # Merge all slices !
#    out_merge_mean = merge(in_files=resampled_warped_slices, dimension='z')

    # slice functional
    sliced_func_filenames = []
    for slice_n in range(n_slices):
        out_slicer = slicer(in_file=allineated_filename,
                            keep='{0} {1}'.format(slice_n, slice_n),
                            out_file=fname_presuffix(allineated_filename,
                                                     suffix='Sl%d' % slice_n))
        _ = fix_obliquity(out_slicer.outputs.out_file,
                          allineated_filename)
        sliced_func_filenames.append(out_slicer.outputs.out_file)

    # Apply the precomputed warp slice by slice
    warped_func_slices = []
    for (sliced_func_filename, warp_filename) in zip(
            sliced_func_filenames, warp_filenames):
        out_warp_apply = warp_apply(in_file=sliced_func_filename,
                                    master=sliced_func_filename,
                                    warp=warp_filename,
                                    out_file=fname_presuffix(
                                        sliced_func_filename, suffix='_qw'))
        warped_func_slices.append(out_warp_apply.outputs.out_file)

    # Fix the obliquity
    for (sliced_func_filename, warped_func_slice) in zip(
            sliced_func_filenames, warped_func_slices):
        _ = fix_obliquity(warped_func_slice, sliced_func_filename)

    # Finally, merge all slices !
    out_merge_func = merge(in_files=warped_func_slices, dimension='z')
#    out_merge_anat = merge(in_files=resampled_registered_anat_filenames,
#                           dimension='z')

    # Update the fmri data
    setattr(session_data, "coreg_func_", out_merge_func.outputs.merged_file)
    setattr(session_data, "coreg_anat_", registered_anat_oblique_filename)
    setattr(session_data, "coreg_transform_", transform_filename)
    os.chdir(current_dir)


def _func_to_template(func_coreg_filename, template_filename, write_dir,
                      func_to_anat_oned_filename,
                      anat_to_template_oned_filename,
                      anat_to_template_warp_filename,
                      caching=False):
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
        Path to the affine 1D transform from anatomical to template space.
    """
    if caching:
        memory = Memory(write_dir)
        warp_apply = memory.cache(afni.NwarpApply)
    else:
        warp_apply = afni.NwarpApply().run

    current_dir = os.getcwd()
    os.chdir(write_dir)
    out_warp_apply = warp_apply(in_file=func_coreg_filename,
                                master=template_filename,
                                nwarp=[anat_to_template_warp_filename,
                                       anat_to_template_oned_filename,
                                       func_to_anat_oned_filename],
                                out_file=fname_presuffix(
                                    func_coreg_filename, suffix='_normalized'))

    os.chdir(current_dir)

    return out_warp_apply.outputs.out_file


def fmri_sessions_to_template(sessions_data, t_r, template_filename, write_dir,
                              prior_rigid_body_registration=False,
                              slice_timing=True,
                              caching=False):
    """ Registration of subject's functional and anatomical images to
    a given template.

    Parameters
    ----------
    sessions_data : sequence of sammba.registration.SessionData
        Animals data, giving paths to their functional and anatomical images.

    t_r : float
        Repetition time for the EPI, in seconds.

    template_filename : str
        Template to register the functional to.

    write_dir : str
        Path to the affine 1D transform from anatomical to template space.

    caching : bool, optional
        Wether or not to use caching.

    Returns
    -------
    the same sequence with each animal_data updated: attributes
        `registered_func_` and
        `registered_anat_` are added to specify the paths to the functional
        and anatomical images registered to the template,
        and `output_dir_` is added to give output directory for each animal.
    """
    if not hasattr(sessions_data, "__iter__"):
            raise ValueError(
                "'animals_data' input argument must be an iterable. You "
                "provided {0}".format(sessions_data.__class__))

    for n, data in enumerate(sessions_data):
        if not isinstance(data, FMRISession):
            raise ValueError('Each animal data must have type '
                             'sammba.registration.Animal. You provided {0} of'
                             ' type {1}'.format(data, type(data)))
        if data.animal_id is None:
            setattr(data, "animal_id", 'animal{0:03d}'.format(n + 1))

    # Check that ids are different
    animals_ids = [data.animal_id for data in sessions_data]
    if len(set(animals_ids)) != len(animals_ids):
        raise ValueError('Animals ids must be different. You'
                         ' provided {0}'.format(animals_ids))

    coreg_func_filenames = []
    func_to_anat_oned_filenames = []
    for n, animal_data in enumerate(sessions_data):
        animal_data._check_inputs()
        animal_output_dir = os.path.join(os.path.abspath(write_dir),
                                         animal_data.animal_id)
        setattr(animal_data, "output_dir_", animal_output_dir) # XXX do a function for creating new attributes ?

        sessions_data[n] = animal_data
        if not os.path.isdir(animal_data.output_dir_):
            os.makedirs(animal_data.output_dir_)

        coreg_func_filename, _, func_to_anat_oned_filename = \
            coregister_fmri_session(
                animal_data, t_r, write_dir,
                prior_rigid_body_registration=prior_rigid_body_registration,
                slice_timing=slice_timing,
                caching=caching)
        coreg_func_filenames.append(coreg_func_filename)
        func_to_anat_oned_filenames.append(func_to_anat_oned_filename)

    anat_filenames = [animal_data.anat for animal_data in sessions_data]
    (normalized_anat_filenames, anat_to_template_oned_filenames,
     anat_to_template_warp_filenames) = \
        anats_to_template(anat_filenames, template_filename,
                          animal_data.output_dir_,
                          caching=caching)
    for n, (animal_data, normalized_anat_filename) in enumerate(zip(
            sessions_data, normalized_anat_filenames)):
        setattr(animal_data, "registered_anat_", normalized_anat_filename)
        sessions_data[n] = animal_data

    for n, (animal_data, coreg_func_filename, func_to_anat_oned_filename,
        anat_to_template_oned_filename,
        anat_to_template_warp_filename) in enumerate(zip(sessions_data,
            coreg_func_filenames, func_to_anat_oned_filenames,
            anat_to_template_oned_filenames,
            anat_to_template_warp_filenames)):
        normalized_func_filename = _func_to_template(
            coreg_func_filename,
            template_filename,
            animal_data.output_dir_,
            func_to_anat_oned_filename,
            anat_to_template_oned_filename,
            anat_to_template_warp_filename,
            caching=caching)

        setattr(animal_data, "registered_func_", normalized_func_filename)
        sessions_data[n] = animal_data

    return sessions_data
