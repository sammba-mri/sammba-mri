"""
Template creation
=================

Here we show how to create a template from multiple anatomical scans.
"""
import numpy as np
from nipype.interfaces import afni, fsl, ants
from nipype.caching import Memory
from nipype.utils.filemanip import fname_presuffix


def save_inverted_affine(oned_matrix_file, inverted_matrix_file):
    u11, u12, u13, v1, u21, u22, u23, v2, u31, u32, u33, v3 = np.loadtxt(
        oned_matrix_file)
    u = np.array([[u11, u12, u13], [u21, u22, u23], [u31, u32, u33]])
    v = np.array([v1, v2, v3])
    u_new = np.linalg.inv(u)
    v_new = -u_new.dot(v)
    np.savetxt(inverted_matrix_file,
               np.hstack((u_new[0], v_new[:1], u_new[1], v_new[1:2],
                          u_new[2], v_new[2:])))


def _compose_affines(in_affines):
    if len(in_affines) > 1:
        composed_affine = np.zeros((3, 4))
        previous_composed_affine = _compose_affines(in_affines[:-1])
        composed_affine[:, :3] = in_affines[1][:, :3].dot(
            previous_composed_affine[:, :3])
        composed_affine[:, 3] = in_affines[1][:, 3] +\
            in_affines[1][:, :3].dot(previous_composed_affine[:, 3])
        return composed_affine
    elif len(in_affines) == 1:
        return in_affines[0]
    else:
        raise ValueError('Can not compose empty list of affines')

import os


def save_composed_affines(in_files, out_file):
    in_affines = []
    for in_affine_file in in_files:
        u11, u12, u13, v1, u21, u22, u23, v2, u31, u32, u33, v3 = np.loadtxt(
            in_affine_file)
        u = np.array([[u11, u12, u13], [u21, u22, u23], [u31, u32, u33]])
        v = np.array([[v1], [v2], [v3]])
        in_affines.append(np.hstack((u, v)))

    composed_affine = _compose_affines(in_affines)
    out_file = os.path.abspath(out_file)
    np.savetxt(out_file, np.hstack((composed_affine[0], composed_affine[1],
                                    composed_affine[2])))
    return out_file


def save_stacked_volumes(volumes):
    """ Computes a video of stacked volumes, as well as their average.
    """

# Retrieve the data
from sammba import data_fetchers

retest = data_fetchers.fetch_zurich_test_retest(subjects=range(3))
# Initially, registration is of extracted brains. Once these are reasonably
# aligned, whole heads are registered, weighted by masks that, if parameters
# are chosen well, include some scalp. The amount of scalp is hopefully enough
# to help in differentiating the brain-scalp boundary without including so much
# head tissue that it starts to contaminate the registration with the highly
# variable head tissue.

##############################################################################
# Register using center of mass
# -----------------------------
# An initial coarse registration is done using brain centre of mass (CoM).

conv=0.005
twoblur=1
brainvol=400

memory = Memory('/tmp')

unifize = memory.cache(afni.Unifize)
clip_level = memory.cache(afni.ClipLevel)
#RATS_MM "$anat"_"$Bc".nii.gz "$anat"_"$Bc"Bm.nii.gz -v $brainvol -t $thresh
bet = memory.cache(fsl.BET)
apply_mask = memory.cache(fsl.ApplyMask)
center_mass = memory.cache(afni.CenterMass)
refit = memory.cache(afni.Refit)

unifized_files = []
unifized_masked_files = []

##############################################################################
# Loop through anatomical scans
for anat_filename in retest.anat:
    # Bias-correction. Note that rbt values might be improveable
    out_unifize = unifize(in_file=anat_filename, rbt=(18.3, 70, 80),
                          out_file='Un.nii.gz')

    # Brain extraction, aided by an approximate guessed brain volume
    out_clip_level = clip_level(in_file=out_unifize.outputs.out_file)
    out_bet = bet(in_file=out_unifize.outputs.out_file)
    out_apply_mask = apply_mask(in_file=out_unifize.outputs.out_file,
                                mask_file=out_bet.outputs.out_file,
                                out_file='UnBmBe.nii.gz')

    # Set the NIfTI image centre (as defined in the header) to the CoM
    # of the extracted brain
    out_center_mass = center_mass(in_file=out_apply_mask.outputs.out_file,
                                  cm_file='/tmp/cm.txt',
                                  set_cm=(0, 0, 0))
    unifized_masked_files.append(out_center_mass.outputs.out_file)
    out_refit = refit(in_file=out_unifize.outputs.out_file,
                      duporigin_file=out_center_mass.outputs.out_file)
    unifized_files.append(out_refit.outputs.out_file)

##############################################################################
# Quality check videos and template
tcat = memory.cache(afni.TCat)
tstat = memory.cache(afni.TStat)
out_tcat = tcat(in_files=unifized_files, out_file='Un_video.nii.gz')
_ = tstat(in_file=out_tcat.outputs.out_file, out_file='Un_mean.nii.gz')
out_tcat = tcat(in_files=unifized_masked_files, out_file='UnBmBe_video.nii.gz')
out_tstat = tstat(in_file=out_tcat.outputs.out_file,
                  out_file='UnBmBe_mean.nii.gz')

##############################################################################
# Create an empty template, with origin placed at CoM
undump = memory.cache(afni.Undump)
out_undump = undump(in_file=out_tstat.outputs.out_file,
                    out_file='emptytemplate.nii.gz')
out_refit = refit(in_file=out_undump.outputs.out_file,
                  xorigin='cen', yorigin='cen', zorigin='cen')

##############################################################################
# The head is then shifted within the image to place the CoM at the image
# center.
resample = memory.cache(afni.Resample)
resampled_files = []
for unifized_file in unifized_files:
    out_resample = resample(in_file=unifized_file,
                            master=out_refit.outputs.out_file,
                            out_file='UnCC.nii.gz')
    resampled_files.append(out_resample.outputs.out_file)

masked_resampled_files = []
for masked_file in unifized_masked_files:
    out_resample = resample(in_file=masked_file,
                            master=out_refit.outputs.out_file,
                            out_file='UnBmBeCC.nii.gz')
    masked_resampled_files.append(out_resample.outputs.out_file)

##############################################################################
# Quality check videos and template
out_tcat = tcat(in_files=resampled_files, out_file='UnCC_video.nii.gz')
_ = tstat(in_file=out_tcat.outputs.out_file, out_file='UnCC_mean.nii.gz')
out_tcat = tcat(in_files=masked_resampled_files, out_file='UnBmBeCC_video.nii.gz')
out_tstat_bm = tstat(in_file=out_tcat.outputs.out_file,
                     out_file='UnBmBeCC_mean.nii.gz')

##############################################################################
# At this point, we achieved a translation-only registration of the raw
# anatomical images to each other's brain's (as defined by the brain extractor)
# CoMs.

##############################################################################
# Shift rotation
# --------------
# Now we move to rigid-body registration of CoM brains, and application of this
# registration to CoM heads. This registration requires a target template.
# Here we use mean of all bias-corrected, brain-extracted, mass-centered
# images. Other possibilities include an externally-sourced image or, more
# biased, a nicely-aligned individual.

allineate = memory.cache(afni.Allineate)
masked_allineated_files = []
transform_matrices = []
for masked_resampled_file in masked_resampled_files:
    out_allineate = allineate(in_file=masked_resampled_file,
                              reference=out_tstat_bm.outputs.out_file,
                              out_matrix='UnBmBeCCAl3.aff12.1D',
                              convergence=0.005,
                              two_blur=1,
                              warp_type='shift_rotate',
                              out_file='UnBmBeCCAl3.nii.gz')
    transform_matrices.append(out_allineate.outputs.out_matrix)
    masked_allineated_files.append(out_allineate.outputs.out_file)

##############################################################################
# Application to the whole head image. can also be used for a good
# demonstration of linear vs. non-linear registration quality
allineated_files = []
for resampled_file, transform_matrix in zip(resampled_files,
                                            transform_matrices):
    out_allineate = allineate(in_file=resampled_file,
                              master=out_tstat_bm.outputs.out_file,
                              in_matrix=transform_matrix,
                              out_file='UnBmBeCCAa3.nii.gz')
    allineated_files.append(out_allineate.outputs.out_file)

##############################################################################
# Quality check videos and template
out_tcat = tcat(in_files=masked_allineated_files, out_file='shr_video.nii.gz')
_ = tstat(in_file=out_tcat.outputs.out_file, out_file='shr_mean.nii.gz')

##############################################################################
# Note that this rigid body registration may need to be run more than once.
out_tcat = tcat(in_files=allineated_files, out_file='shr_videohead.nii.gz')
out_tstat_shr = tstat(in_file=out_tcat.outputs.out_file,
                      out_file='shr_meanhead.nii.gz')
mask_tool = memory.cache(afni.MaskTool)
out_mask_tool = mask_tool(in_file=out_tcat.outputs.out_file, count=True,
                          out_file='shr_count.nii.gz')
# The count mask is useful for looking at brain extraction efficiency
# and differences in brain size.

##############################################################################
# Affine
# ------
# reducing reslice errors in the final result, and 
# 4) there are two n options: the first is to identify the MRIT3_shr.bash run,
# the second to give a number to the affine results to make them easier to tell
# apart from the previous linear registrations. It is important to note that
# there is not a unique weight mask for each head. Instead, since heads are
# already roughly aligned, a single generous weight mask can be used for all of
# them (typically the count mask produced by the previous MRIT3_shr.bash),
# which should also cover some scalp for most individuals. If brains vary
# massively in size this may not work well, but from experience it is fine.

# In any case, if size is that variable, the extraction (which assumes brain
#  are around a certain size) would not have worked well either, messing up
# the entire strategy. Doing instead an affine transform of extracted brains
# (or whole heads weighted by an individual-specific mask presumably based on
# brain extraction results) is worse for two reasons. Firstly, it is less
# robust if extraction is sometimes poor. The use of a mask calculated from the
# results of the whole group effectively provides safety in numbers. Secondly,
# it tends to make them too large, at least with 3dAllineate on small mammal
# data in the context of how it is being processed here. This is acceptable
# for mutual registration but gives the final template an incorrect size.
# An alternative is to skip MRIT4_aff.bash and go straight to MRIT5_Qw.bash,
# which seems robust enough to handle the rougher inputs from MRIT3_shr.bash,
# and does not inflate average brain size.

# We begin by achieving an affine registration on aligned heads.
# A weighting mask is used to ...
composed_affine_files = []
allineated_head_files = []
catmatvec = memory.cache(afni.CatMatvec)
for allineated_file, transform_matrix in zip(allineated_files,
                                             transform_matrices):
    out_allineate = allineate(in_file=allineated_file,
                              reference=out_tstat_shr.outputs.out_file,
                              out_matrix='UnBmBeAl4.aff12.1D',
                              convergence=0.005,
                              two_blur=1,
                              one_pass=True,
                              weight=out_mask_tool.outputs.out_file,
                              out_file='UnCCAl4.nii.gz')
    allineated_head_files.append(out_allineate.outputs.out_file)
    catmatvec_out_file = os.path.join(os.path.dirname(transform_matrix),
                                      'UnBmBeCCAl3UnCCAl4.aff12.1D')
    out_catmatvec = catmatvec(in_file=[(transform_matrix, 'ONELINE'),
                                       (out_allineate.outputs.out_matrix,
                                        'ONELINE')],
                              out_file=catmatvec_out_file)
    composed_affine_files.append(catmatvec_out_file)

##############################################################################
# Each resulting registration matrix is concatenated to the corresponding
# rigid bory registration matrix then directly applied to the CoM brain
# and head, reducing reslice errors in the final result.
masked_aa4 = []
for masked_allineated_file, composed_affine_file in zip(
        masked_allineated_files, composed_affine_files):
    out_allineate = allineate(in_file=masked_allineated_file,
                              master=out_tstat_shr.outputs.out_file,
                              in_matrix=composed_affine_file,
                              out_file='UnBmBeCCAa4.nii.gz')
    masked_aa4.append(out_allineate.outputs.out_file)

##############################################################################
# The application to the whole head image can also be used for a good
# demonstration of linear vs. non-linear registration quality.
unifized_aa4 = []
for allineated_head_file, composed_affine_file in zip(allineated_head_files,
                                                      composed_affine_files):
    out_allineate = allineate(in_file=allineated_head_file,
                              master=out_tstat_shr.outputs.out_file,
                              in_matrix=composed_affine_file,
                              out_file='UnCCAa.nii.gz')
    unifized_aa4.append(out_allineate.outputs.out_file)

##############################################################################
# Quality check videos and template
out_tcat = tcat(in_files=masked_aa4, out_file='aff4_video.nii.gz')
_ = tstat(in_file=out_tcat.outputs.out_file, out_file='aff4_mean.nii.gz')
out_tcat_head = tcat(in_files=unifized_aa4, out_file='aff4_videohead.nii.gz')
out_tstat_affine_head = tstat(in_file=out_tcat_head.outputs.out_file,
                              out_file='aff4_meanhead.nii.gz')

##############################################################################
# Non-linear registration
# -----------------------
# A weight mask that extends beyond the brain, incorporating some
# surrounding tissue, is needed to help better define the brain head boundary.
out_mask_tool = mask_tool(in_file=out_tcat.outputs.out_file, count=True,
                          out_file='aff4_count.nii.gz')
out_mask_tool = mask_tool(in_file=out_tcat.outputs.out_file, union=True,
                          out_file='aff4_unionmask.nii.gz')
out_mask_tool = mask_tool(in_file=out_mask_tool.outputs.out_file,
                          dilate_inputs='4',
                          out_file='aff4_unionmaskdil4.nii.gz')

##############################################################################
# The input source images are initially transformed prior to registration,
# to ensure that they are already quite well-aligned to the template.
nwarp_cat = memory.cache(afni.NwarpCat)
qwarp = memory.cache(afni.Qwarp)
warped_files = []
warp_files = []
for composed_affine_file, aa4 in zip(composed_affine_files, unifized_aa4):
    out_qwarp = qwarp(in_file=aa4,
                      base_file=out_tstat_affine_head.outputs.out_file,
                      nmi=True,
                      noneg=True,
                      iwarp=True,
                      weight=out_mask_tool.outputs.out_file,
                      iniwarp=[composed_affine_file],
                      inilev=0,
                      maxlev=4,
                      out_file='UnCCQw1.nii.gz')
    warp_files.append(out_qwarp.outputs.source_warp)
    warped_files.append(out_qwarp.outputs.warped_source)

out_tcat = tcat(in_files=warped_files, out_file='Qw1_videohead.nii.gz')
out_tstat_warp_head = tstat(in_file=out_tcat.outputs.out_file,
                            out_file='Qw1_meanhead.nii.gz')

##############################################################################
# Then iterative registration from a given level to another is achieved.
# Note that any level below a patch size of 25 will not be done (see 3dQwarp
# help for further detail).
# The input transform is the former warp and needs to be concatenated to IDENT
# initially; I forget why, I think it is to avoid some weird bug.
warped_files2 = []
warp_files2 = []
for warp_file, aa4_file in zip(warp_files, unifized_aa4):
    out_nwarp_cat = nwarp_cat(
        in_files=[('IDENT', out_tstat_warp_head.outputs.out_file),
                  warp_file], out_file='iniwarp.nii.gz')
    out_qwarp = qwarp(in_file=aa4_file,
                      base_file=out_tstat_warp_head.outputs.out_file,
                      nmi=True,
                      noneg=True,
                      iwarp=True,
                      weight=out_mask_tool.outputs.out_file,
                      iniwarp=[out_nwarp_cat.outputs.out_file],
                      inilev=5,
                      maxlev=7,
                      out_file='UnCCQw2.nii.gz')
    warp_files2.append(out_qwarp.outputs.source_warp)
    warped_files2.append(out_qwarp.outputs.warped_source)

out_tcat = tcat(in_files=warped_files2, out_file='Qw2_videohead.nii.gz')
out_tstat_warp_head2 = tstat(in_file=out_tcat.outputs.out_file,
                             out_file='Qw2_meanhead.nii.gz')

# Such possibilities can be exploited to avoid building up reslice errors. XXX check with Nad

# Warp with mini-patch
warped_files3 = []
warp_files3 = []
for warp_file, aa4_file in zip(warp_files2, unifized_aa4):
    out_qwarp = qwarp(in_file=aa4_file,
                      base_file=out_tstat_warp_head2.outputs.out_file,
                      nmi=True,
                      noneg=True,
                      iwarp=True,
                      weight=out_mask_tool.outputs.out_file,
                      iniwarp=[warp_file],
                      inilev=8,
                      minpatch=13,
                      out_file='UnCCQw3.nii.gz')
    warped_files3.append(out_qwarp.outputs.warped_source)
    warp_files3.append(out_qwarp.outputs.source_warp)

out_tcat = tcat(in_files=warped_files3, out_file='Qw3_videohead.nii.gz')
out_tstat_warp_head3 = tstat(in_file=out_tcat.outputs.out_file,
                             out_file='Qw3_meanhead.nii.gz')

# Final warp
warped_files4 = []
for warp_file, aa4_file in zip(warp_files3, unifized_aa4):
    out_qwarp = qwarp(in_file=aa4_file,
                      base_file=out_tstat_warp_head3.outputs.out_file,
                      nmi=True,
                      noneg=True,
                      iwarp=True,
                      weight=out_mask_tool.outputs.out_file,
                      iniwarp=[warp_file],
                      inilev=11,
                      minpatch=9,
                      out_file='UnCCQw4.nii.gz')
    warped_files4.append(out_qwarp.outputs.warped_source)

out_tcat = tcat(in_files=warped_files4, out_file='Qw4_videohead.nii.gz')
out_tstat_warp_head4 = tstat(in_file=out_tcat.outputs.out_file,
                             out_file='Qw4_meanhead.nii.gz')

##############################################################################
# Registration to template
# ------------------------
# Apply warp to anat and atlas
# Apply non-linear registration results to uncorrected images
warp_apply = memory.cache(afni.NwarpApply)
for resampled_file, warp_file in zip(resampled_files, warp_files):
    out_warp_apply = warp_apply(in_file=resampled_file,
                                warp=warp_file,
                                master=out_tstat.outputs.out_file,
                                out_file='Na.nii.gz')

# Transform the atlas from template space to individual spaces
# This can then be useful for morphological measurements in each individual.
for resampled_file, warp_file in zip(resampled_files, warp_files):
    out_nwarp_cat = nwarp_cat(in_files=[('INV', warp_file),
                                        atlas_warp_file],  # check with Nad
                              out_file='atlas_WARPINV.nii.gz')
    out_warp_apply = warp_apply(in_file=resampled_file,
                                warp=warp_file,
                                master=out_tstat.outputs.out_file,
                                ainterp='NN',
                                out_file='atlas_Na.nii.gz')