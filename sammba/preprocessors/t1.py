import os
from sammba.externals.nipype.interfaces import afni, fsl
from sammba.interfaces import RatsMM
from sammba.externals.nipype.utils.filemanip import fname_presuffix
from sammba.externals.nipype.caching import Memory
from sklearn.datasets.base import Bunch


def register_to_common(t1_filenames, write_dir, common_filename=None,
                       brain_volume=400, registration_kind='nonlinear',
                       nonlinear_iters=3, nonlinear_levels=[1, 2, 3],
                       nonlinear_minimal_patch=75,
                       convergence=0.005, caching=False):
    """
    The functional volume is aligned to the anatomical, first with a rigid body
    registration and then on a per-slice basis (only a fine correction, this is
    mostly for correction of EPI distortion). This pipeline includes
    slice timing.

    Parameters
    ----------
    pipeline_graph : one of {'orig', 'hierarchical', 'flat', 'exec', 'colored'}
        'orig': creates a top level graph without expanding internal
        workflow nodes;
        'flat': expands workflow nodes recursively;
        'hierarchical': expands workflow nodes recursively with a notion
        on hierarchy;
        'colored': expands workflow nodes recursively with a
        notion on hierarchy in color;
        'exec': expands workflows to depict iterables
    """
    registration_kinds = ['rigid-body', 'affine', 'nonlinear']
    if registration_kind not in registration_kinds:
        raise ValueError(
            'Registration kind must be one of {0}, you entered {1}'.format(
                registration_kinds, registration_kind))
    if caching:
        memory = Memory(write_dir)
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
    else:
        unifize = afni.Unifize().run
        clip_level = afni.ClipLevel().run
        rats = RatsMM().run
        apply_mask = fsl.ApplyMask().run
        center_mass = afni.CenterMass().run
        refit = afni.Refit().run
        tcat = afni.TCat().run
        tstat = afni.TStat().run
        undump = afni.Undump().run
        resample = afni.Resample().run
        allineate = afni.Allineate().run
        allineate2 = afni.Allineate().run
        mask_tool = afni.MaskTool().run
        catmatvec = afni.CatMatvec().run
        qwarp = afni.Qwarp().run
        nwarp_cat = afni.NwarpCat().run
        warp_apply = afni.NwarpApply().run

    current_dir = os.getcwd()
    os.chdir(write_dir)
    ###########################################################################
    # Register using center of mass
    # -----------------------------
    # An initial coarse registration is done using brain centre of mass (CoM).
    #
    # First we loop through anatomical scans and correct intensities for bias.
    unifized_files = []
    for n, anat_file in enumerate(t1_filenames):
        out_file = fname_presuffix(anat_file,
                                   suffix='_{}_corrected'.format(n))
        out_file = os.path.join(write_dir, os.path.basename(out_file))
        out_unifize = unifize(in_file=anat_file, out_file=out_file)
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
            intensity_threshold=int(out_clip_level.outputs.clip_val))
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
    out_tcat = tcat(in_files=brain_files, outputtype='NIFTI_GZ')

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

    if registration_kind == 'rigid-body':
        os.chdir(write_dir)
        return Bunch(registered_anats=shift_rotated_head_files,
                     transforms=rigid_transform_files)

    mask_node = pe.Node(afni.MaskTool(), name='generate_count_mask')
    allineate_node3 = pe.Node(afni.Allineate(), name='affine_register')
    catmatvec_node = pe.Node(afni.CatMatvec(), name='concatenate_transforms')
    allineate_node4 = pe.Node(afni.Allineate(), name='apply_transform_on_head')
    allineate_node5 = pe.Node(afni.Allineate(),
                              name='apply_transform_on_brain')
    tcat_node3 = pe.Node(
        afni.TCat(), name='concatenate_affine_registred_across_individuals')
    tstat_node3 = pe.Node(afni.TStat(),
                          name='compute_affine_registred_average')

    workflow.add_nodes([mask_node, allineate_node3, catmatvec_node,
                        allineate_node4, allineate_node5,
                        tcat_node3, tstat_node3])

    workflow.connect(tcat_node2, 'out_file', mask_node, 'in_file')
    workflow.connect(mask_node, 'out_file', allineate_node3, 'weight')
    workflow.connect(allineate_node2, 'out_file', allineate_node3, 'in_file')
    workflow.connect(allineate_node3, 'out_matrix', catmatvec_node, 'in_file')
    #XXX how can we enter multiple files ? 
    workflow.connect(catmatvec_node, 'out_file', allineate_node4, 'in_matrix')
    workflow.connect(allineate_node, 'out_file', allineate_node4, 'in_file')

    workflow.connect(unifize_node, 'out_file', catmatvec_node, 'in_file')
    workflow.connect(rats_node, 'out_file', apply_mask_node, 'mask_file')
    workflow.connect(apply_mask_node, 'out_file',
                     center_mass_node, 'in_file')
    workflow.connect(unifize_node, 'out_file', refit_node, 'in_file')
    workflow.connect(center_mass_node, 'out_file',
                     refit_node, 'duporigin_file')
    workflow.connect(center_mass_node, 'out_file', tcat_node, 'in_files')
    workflow.connect(tcat_node, 'out_file', tstat_node, 'in_file')
    workflow.connect(tstat_node, 'out_file', undump_node, 'in_file')
    workflow.connect(undump_node, 'out_file', refit_node2, 'in_file')
    workflow.connect(refit_node, 'out_file', resample_node, 'master')
    workflow.connect(refit_node2, 'out_file', resample_node, 'in_file')

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
        return Bunch(registered_anats=allineated_head_files,
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
    if nonlinear_iters > 1:
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
    if nonlinear_iters > 2:
        for n_iter in range(2, nonlinear_iters):
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
                    inilev=nonlinear_levels[n_iter],
                    minpatch=nonlinear_minimal_patch,
                    out_file=fname_presuffix(
                        shifted_head_file, suffix='_warped{}'.format(n_iter)))
                warped_files.append(out_qwarp.outputs.warped_source)
                warp_files.append(out_qwarp.outputs.source_warp)

            out_tcat = tcat(
                in_files=warped_files,
                out_file=os.path.join(write_dir,
                                      'warped_{0}iters_heads.nii.gz'.format(
                                      n_iter)))
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
            out_file=fname_presuffix(head_file,
                                     suffix='_warped{}'.format(nonlinear_iters)
                                     ))
        warped_files.append(out_warp_apply.outputs.out_file)

    os.chdir(current_dir)
    return Bunch(registered_anats=warped_files,
                 transforms=warp_files)
