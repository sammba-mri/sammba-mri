import os
from sammba.externals.nipype.interfaces import afni, fsl
from sammba.interfaces import RatsMM
from sammba.externals.nipype.utils.filemanip import fname_presuffix
from sammba.externals.nipype.caching import Memory
from sklearn.datasets.base import Bunch


def register_to_common(t1_filenames, write_dir, brain_volume=400,
                       registration_kind='affine',
                       nonlinear_levels=[1, 2, 3],
                       nonlinear_minimal_patch=75,
                       convergence=0.005, caching=False, verbose=0):
    """ Create common template from native T1 weighted images and achieve
    their registration to it.

    Parameters
    ----------
    t1_filenames : list of str
        Paths to the T1 weighted images.

    write_dir : Str
        Path to an existant directory to save output files to.

    brain_volume : float, optional
        Volumes of the brain as passed to Rats_MM brain extraction tool.
        Default to 400 for the mouse brain.

    registration_kind : one of {'rigid', 'affine', 'nonlinear'}, optional
        The allowed transform kind.

    nonlinear_levels : list of int, optional
        Maximal levels for each nonlinear warping iteration. Passed iteratively
        to sammba.externals.nipype.interfaces.afni.Qwarp

    nonlinear_minimal_patch : int, optional
        Minimal patch for the final nonlinear warp, passed to
        sammba.externals.nipype.interfaces.afni.Qwarp

    caching : bool, optional
        If True, caching is used for all the registration steps.

    convergence : float, optional
        Convergence limit, passed to

    verbose : bool, optional
        If True, all steps are verbose.

    Returns
    -------
    data : sklearn.datasets.base.Bunch
        Dictionary-like object, the interest attributes are :

        - 'registered': string list. Paths to registered images.
        - 'transforms': string list. Paths to the transforms from the raw
                        images to the registered images.
    """
    registration_kinds = ['rigid', 'affine', 'nonlinear']
    if registration_kind not in registration_kinds:
        raise ValueError(
            'Registration kind must be one of {0}, you entered {1}'.format(
                registration_kinds, registration_kind))

    if verbose:
        terminal_output = 'allatonce'
    else:
        terminal_output = 'none'
    if caching:

        memory = Memory(write_dir)
        copy = memory.cache(afni.Copy)
        unifize = memory.cache(afni.Unifize)
        clip_level = memory.cache(afni.ClipLevel)
        rats = memory.cache(RatsMM)
        apply_mask = memory.cache(fsl.ApplyMask)
        center_mass = memory.cache(afni.CenterMass)
        refit = memory.cache(afni.Refit)
        tcat = memory.cache(afni.TCat)
        tstat = memory.cache(afni.TStat)
        undump = memory.cache(afni.Undump)
        resample = memory.cache(afni.Resample)
        allineate = memory.cache(afni.Allineate)
        allineate2 = memory.cache(afni.Allineate)
        mask_tool = memory.cache(afni.MaskTool)
        catmatvec = memory.cache(afni.CatMatvec)
        qwarp = memory.cache(afni.Qwarp)
        nwarp_cat = memory.cache(afni.NwarpCat)
        warp_apply = memory.cache(afni.NwarpApply)
        for func in [copy, unifize, rats, apply_mask, center_mass, refit,
                     tcat, tstat, undump, resample, allineate, allineate2,
                     mask_tool, catmatvec, qwarp, nwarp_cat, warp_apply]:
            func.interface().set_default_terminal_output(terminal_output)
    else:
        copy = afni.Copy(terminal_output=terminal_output).run
        unifize = afni.Unifize(terminal_output=terminal_output).run
        clip_level = afni.ClipLevel().run  # XXX fix nipype bug with 'none'
        rats = RatsMM(terminal_output=terminal_output).run
        apply_mask = fsl.ApplyMask(terminal_output=terminal_output).run
        center_mass = afni.CenterMass(terminal_output=terminal_output).run
        refit = afni.Refit(terminal_output=terminal_output).run
        tcat = afni.TCat(terminal_output=terminal_output).run
        tstat = afni.TStat(terminal_output=terminal_output).run
        undump = afni.Undump(terminal_output=terminal_output).run
        resample = afni.Resample(terminal_output=terminal_output).run
        allineate = afni.Allineate(terminal_output=terminal_output).run
        allineate2 = afni.Allineate(terminal_output=terminal_output).run
        mask_tool = afni.MaskTool(terminal_output=terminal_output).run
        catmatvec = afni.CatMatvec(terminal_output=terminal_output).run
        qwarp = afni.Qwarp(terminal_output=terminal_output).run
        nwarp_cat = afni.NwarpCat(terminal_output=terminal_output).run
        warp_apply = afni.NwarpApply(terminal_output=terminal_output).run

    current_dir = os.getcwd()
    os.chdir(write_dir)

    ###########################################################################
    # First copy anatomical files, to make sure they are never changed
    # and they have different names across individuals
    copied_t1_filenames = []
    for n, anat_file in enumerate(t1_filenames):
        suffixed_file = fname_presuffix(anat_file, suffix='_{}'.format(n))
        copied_file = os.path.join(write_dir, os.path.basename(suffixed_file))
        out_copy = copy(in_file=anat_file, out_file=copied_file)
        copied_t1_filenames.append(out_copy.outputs.out_file)

    ###########################################################################
    # Register using center of mass
    # -----------------------------
    # An initial coarse registration is done using brain centre of mass (CoM).
    #
    # First we loop through anatomical scans and correct intensities for bias.
    unifized_files = []
    for n, anat_file in enumerate(copied_t1_filenames):
        out_unifize = unifize(in_file=anat_file, outputtype='NIFTI_GZ')
        unifized_files.append(out_unifize.outputs.out_file)

    ###########################################################################
    # Second extract brains, aided by an approximate guessed brain volume,
    # and set the NIfTI image centre (as defined in the header) to the CoM
    # of the extracted brain.
    brain_files = []
    for unifized_file in unifized_files:
        out_clip_level = clip_level(in_file=unifized_file)
        out_rats = rats(
            in_file=unifized_file,
            volume_threshold=brain_volume,
            intensity_threshold=int(out_clip_level.outputs.clip_val),
            terminal_output=terminal_output)
        out_apply_mask = apply_mask(in_file=unifized_file,
                                    mask_file=out_rats.outputs.out_file)
        out_center_mass = center_mass(
            in_file=out_apply_mask.outputs.out_file,
            cm_file=fname_presuffix(unifized_file, suffix='_cm.txt',
                                    use_ext=False),
            set_cm=(0, 0, 0))
        brain_files.append(out_center_mass.outputs.out_file)

    ###########################################################################
    # Same header change, for head files.
    head_files = []
    for unifized_file, brain_file in zip(unifized_files, brain_files):
        out_refit = refit(in_file=unifized_file, duporigin_file=brain_file)
        head_files.append(out_refit.outputs.out_file)

    ###########################################################################
    # The brain files with new image center are concatenated to produce
    # a quality check video
    out_tcat = tcat(in_files=brain_files, outputtype='NIFTI_GZ',
                    terminal_output=terminal_output)

    ###########################################################################
    # and averaged
    out_tstat = tstat(in_file=out_tcat.outputs.out_file, outputtype='NIFTI_GZ')

    ###########################################################################
    # to create an empty template, with origin placed at CoM
    out_undump = undump(in_file=out_tstat.outputs.out_file,
                        outputtype='NIFTI_GZ')
    out_refit = refit(in_file=out_undump.outputs.out_file,
                      xorigin='cen', yorigin='cen', zorigin='cen')

    ###########################################################################
    # Finally, we shift heads and brains within the images to place the CoM at
    # the image center.
    shifted_head_files = []
    for head_file in head_files:
        out_resample = resample(in_file=head_file,
                                master=out_refit.outputs.out_file,
                                outputtype='NIFTI_GZ')
        shifted_head_files.append(out_resample.outputs.out_file)

    shifted_brain_files = []
    for brain_file in brain_files:
        out_resample = resample(in_file=brain_file,
                                master=out_refit.outputs.out_file,
                                outputtype='NIFTI_GZ')
        shifted_brain_files.append(out_resample.outputs.out_file)

    ###########################################################################
    # Quality check videos and average brain
    out_tcat = tcat(in_files=shifted_brain_files,
                    out_file=os.path.join(write_dir, 'shifted_brains.nii.gz'))
    out_tstat_shifted_brain = tstat(in_file=out_tcat.outputs.out_file,
                                    outputtype='NIFTI_GZ')

    ###########################################################################
    # At this point, we achieved a translation-only registration of the raw
    # anatomical images to each other's brain's (as defined by the brain
    # extractor) CoMs.

    ###########################################################################
    # Shift rotate
    # ------------
    # Now we move to rigid-body registration of CoM brains, and application of
    # this registration to CoM heads. This registration requires a target
    #  template. Here we use mean of all bias-corrected, brain-extracted,
    # mass-centered images. Other possibilities include an externally-sourced
    # image or, more biased, a nicely-aligned individual.
    shift_rotated_brain_files = []
    rigid_transform_files = []
    for shifted_brain_file in shifted_brain_files:
        out_allineate = allineate(
            in_file=shifted_brain_file,
            reference=out_tstat_shifted_brain.outputs.out_file,
            out_matrix=fname_presuffix(shifted_brain_file,
                                       suffix='_shr.aff12.1D',
                                       use_ext=False),
            convergence=convergence,
            two_blur=1,
            warp_type='shift_rotate',
            out_file=fname_presuffix(shifted_brain_file, suffix='_shr'))
        rigid_transform_files.append(out_allineate.outputs.out_matrix)
        shift_rotated_brain_files.append(out_allineate.outputs.out_file)

    ###########################################################################
    # Application to the whole head image. can also be used for a good
    # demonstration of linear vs. non-linear registration quality
    shift_rotated_head_files = []
    for shifted_head_file, rigid_transform_file in zip(shifted_head_files,
                                                       rigid_transform_files):
        out_allineate = allineate2(
            in_file=shifted_head_file,
            master=out_tstat_shifted_brain.outputs.out_file,
            in_matrix=rigid_transform_file,
            out_file=fname_presuffix(shifted_head_file, suffix='_shr'))
        shift_rotated_head_files.append(out_allineate.outputs.out_file)

    ###########################################################################
    # Note that this rigid body registration may need to be run more than once.
    # Now we produce an average of rigid body registered heads
    out_tcat = tcat(
        in_files=shift_rotated_head_files,
        out_file=os.path.join(write_dir, 'rigid_body_registered_heads.nii.gz'))
    out_tstat_shr = tstat(in_file=out_tcat.outputs.out_file,
                          outputtype='NIFTI_GZ')

    if registration_kind == 'rigid':
        os.chdir(current_dir)
        return Bunch(registered=shift_rotated_head_files,
                     transforms=rigid_transform_files)

    ###########################################################################
    # Affine transform
    # ----------------
    # We begin by achieving an affine registration on aligned heads.
    # A weighting mask is used to ...
    out_mask_tool = mask_tool(in_file=out_tcat.outputs.out_file,
                              count=True,
                              outputtype='NIFTI_GZ')

    ###########################################################################
    # The count mask is also useful for looking at brain extraction efficiency
    # and differences in brain size.
    affine_transform_files = []
    for shift_rotated_head_file, rigid_transform_file in zip(
            shift_rotated_head_files, rigid_transform_files):
        out_allineate = allineate(
            in_file=shift_rotated_head_file,
            reference=out_tstat_shr.outputs.out_file,
            out_matrix=fname_presuffix(shift_rotated_head_file,
                                       suffix='_affine.aff12.1D',
                                       use_ext=False),
            convergence=convergence,
            two_blur=1,
            one_pass=True,
            weight=out_mask_tool.outputs.out_file,
            out_file=fname_presuffix(shift_rotated_head_file,
                                     suffix='_affine'))
        catmatvec_out_file = fname_presuffix(shift_rotated_head_file,
                                             suffix='_cat.aff12.1D',
                                             use_ext=False)
        out_catmatvec = catmatvec(in_file=[(rigid_transform_file, 'ONELINE'),
                                           (out_allineate.outputs.out_matrix,
                                            'ONELINE')],
                                  out_file=catmatvec_out_file)
        affine_transform_files.append(catmatvec_out_file)

    ###########################################################################
    # Typically, convergence should be set to 0.005 but we increase it for
    # speed reason.

    ###########################################################################
    # Each resulting registration matrix is concatenated to the corresponding
    # rigid bory registration matrix then directly applied to the CoM brain
    # and head, reducing reslice errors in the final result.
    allineated_brain_files = []
    for shifted_brain_file, affine_transform_file in zip(
            shifted_brain_files, affine_transform_files):
        out_allineate = allineate2(
            in_file=shifted_brain_file,
            master=out_tstat_shr.outputs.out_file,
            in_matrix=affine_transform_file,
            out_file=fname_presuffix(shifted_brain_file, suffix='_shr_affine'))
        allineated_brain_files.append(out_allineate.outputs.out_file)

    ###########################################################################
    # The application to the whole head image can also be used for a good
    # demonstration of linear vs. non-linear registration quality.
    allineated_head_files = []
    for shifted_head_file, affine_transform_file in zip(
            shifted_head_files, affine_transform_files):
        out_allineate = allineate2(
            in_file=shifted_head_file,
            master=out_tstat_shr.outputs.out_file,
            in_matrix=affine_transform_file,
            out_file=fname_presuffix(shifted_head_file, suffix='_shr_affine'))
        allineated_head_files.append(out_allineate.outputs.out_file)

    ###########################################################################
    # Quality check videos and template
    out_tcat_head = tcat(
        in_files=allineated_head_files,
        out_file=os.path.join(write_dir, 'affine_registered_heads.nii.gz'))
    out_tstat_allineated_head = tstat(in_file=out_tcat_head.outputs.out_file,
                                      outputtype='NIFTI_GZ')

    if registration_kind == 'affine':
        os.chdir(current_dir)
        return Bunch(registered=allineated_head_files,
                     transforms=affine_transform_files)

    ###########################################################################
    # Non-linear registration
    # -----------------------
    # A weight mask that extends beyond the brain, incorporating some
    # surrounding tissue, is needed to help better define the brain head
    # boundary.
    out_mask_tool = mask_tool(in_file=out_tcat.outputs.out_file, count=True,
                              outputtype='NIFTI_GZ')
    out_mask_tool = mask_tool(in_file=out_tcat.outputs.out_file, union=True,
                              outputtype='NIFTI_GZ')
    out_mask_tool = mask_tool(in_file=out_mask_tool.outputs.out_file,
                              dilate_inputs='4',
                              outputtype='NIFTI_GZ')

    ###########################################################################
    # The input source images are initially transformed prior to registration,
    # to ensure that they are already quite well-aligned to the template.
    # To save time, we only achieve one refinement level per step
    if nonlinear_levels is None:
        nonlinear_levels = [1, 2, 3]

    warped_files = []
    warp_files = []
    for affine_transform_file, shifted_head_file in zip(affine_transform_files,
                                                        shifted_head_files):
        out_qwarp = qwarp(
            in_file=shifted_head_file,
            base_file=out_tstat_allineated_head.outputs.out_file,
            nmi=True,
            noneg=True,
            iwarp=True,
            weight=out_mask_tool.outputs.out_file,
            iniwarp=[affine_transform_file],
            inilev=0,
            maxlev=nonlinear_levels[0],
            out_file=fname_presuffix(shifted_head_file, suffix='_warped1'))
        warp_files.append(out_qwarp.outputs.source_warp)
        warped_files.append(out_qwarp.outputs.warped_source)

    out_tcat = tcat(
        in_files=warped_files,
        out_file=os.path.join(write_dir, 'warped_1iter_heads.nii.gz'))
    out_tstat_warp_head = tstat(in_file=out_tcat.outputs.out_file,
                                outputtype='NIFTI_GZ')

    ###########################################################################
    # Then iterative registration from a given level to another is achieved.
    # Note that any level below a patch size of 25 will not be done (see
    # 3dQwarp help for further detail).
    # The input transform is the former warp and needs to be concatenated to
    # IDENT initially; I forget why, I think it is to avoid some weird bug.
    if len(nonlinear_levels) > 1:
        previous_warp_files = warp_files
        warped_files = []
        warp_files = []
        for warp_file, shifted_head_file in zip(previous_warp_files,
                                                shifted_head_files):
            out_nwarp_cat = nwarp_cat(
                in_files=[('IDENT', out_tstat_warp_head.outputs.out_file),
                          warp_file], out_file='iniwarp.nii.gz')
            out_qwarp = qwarp(
                in_file=shifted_head_file,
                base_file=out_tstat_warp_head.outputs.out_file,
                nmi=True,
                noneg=True,
                iwarp=True,
                weight=out_mask_tool.outputs.out_file,
                iniwarp=[out_nwarp_cat.outputs.out_file],
                inilev=nonlinear_levels[0],
                maxlev=nonlinear_levels[1],
                out_file=fname_presuffix(shifted_head_file, suffix='_warped2'))
            warp_files.append(out_qwarp.outputs.source_warp)
            warped_files.append(out_qwarp.outputs.warped_source)

        out_tcat = tcat(in_files=warped_files,
                        out_file='warped_2iters_heads.nii.gz')
        out_tstat_warp_head = tstat(in_file=out_tcat.outputs.out_file,
                                    outputtype='NIFTI_GZ')

    ###########################################################################
    # Using previous files and concatenated transforms can be exploited to
    # avoid building up reslice errors.
    # Warp with mini-patch
    # In this particular case, minpatch=75 corresponds to a level of 4
    if len(nonlinear_levels) > 2:
        if nonlinear_minimal_patch is None:
            nonlinear_minimal_patch = 75

        for n_iter, inilev in enumerate(nonlinear_levels[2:]):
            previous_warp_files = warp_files
            warped_files = []
            warp_files = []
            for warp_file, shifted_head_file in zip(previous_warp_files,
                                                    shifted_head_files):
                out_qwarp = qwarp(
                    in_file=shifted_head_file,
                    base_file=out_tstat_warp_head.outputs.out_file,
                    nmi=True,
                    noneg=True,
                    iwarp=True,
                    weight=out_mask_tool.outputs.out_file,
                    iniwarp=[warp_file],
                    inilev=inilev,
                    minpatch=nonlinear_minimal_patch,
                    out_file=fname_presuffix(
                        shifted_head_file,
                        suffix='_warped{}'.format(n_iter + 3)))
                warped_files.append(out_qwarp.outputs.warped_source)
                warp_files.append(out_qwarp.outputs.source_warp)

            out_tcat = tcat(
                in_files=warped_files,
                out_file=os.path.join(
                    write_dir, 'warped_{0}iters_heads.nii.gz'.format(n_iter +
                                                                     3)))
            out_tstat_warp_head = tstat(in_file=out_tcat.outputs.out_file,
                                        outputtype='NIFTI_GZ')

    ###########################################################################
    # We can repeat this very last warp while using the last average until we
    # are satisfied with the template quality

    ###########################################################################
    # Register to template
    # --------------------
    # Apply non-linear registration results to uncorrected images
    warped_files = []
    for head_file, warp_file in zip(head_files, warp_files):
        out_warp_apply = warp_apply(
            in_file=head_file,
            warp=warp_file,
            master=out_tstat_warp_head.outputs.out_file,
            out_file=fname_presuffix(
                head_file, suffix='_warped{}'.format(len(nonlinear_levels))
                                     ))
        warped_files.append(out_warp_apply.outputs.out_file)

    os.chdir(current_dir)
    return Bunch(registered=warped_files,
                 transforms=warp_files)
