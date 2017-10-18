"""
Template creation
=================

Here we show how to create a template from multiple anatomical scans.
Initially, registration is of extracted brains. Once these are reasonably
aligned, whole heads are registered, weighted by masks that, if parameters
are chosen well, include some scalp. The amount of scalp is hopefully enough
to help in differentiating the brain-scalp boundary without including so much
head tissue that it starts to contaminate the registration with the highly
variable head tissue.
"""
##############################################################################
# Retrieve the data
from sammba import data_fetchers

retest = data_fetchers.fetch_zurich_test_retest(subjects=range(2))

##############################################################################
# and define the caching repository
from nipype.caching import Memory

memory = Memory('/tmp')

##############################################################################
# Register using center of mass
# -----------------------------
# An initial coarse registration is done using brain centre of mass (CoM).
from nipype.interfaces import afni, fsl

unifize = memory.cache(afni.Unifize)

##############################################################################
# First we loop through anatomical scans and correct intensities for bias.
unifized_files = []
for anat_filename in retest.anat:
    # rbt values might be improveable
    out_unifize = unifize(in_file=anat_filename, rbt=(18.3, 70, 80),
                          out_file='Un.nii.gz')
    unifized_files.append(out_unifize.outputs.out_file)

##############################################################################
# Second extract brains, aided by an approximate guessed brain volume,
# and set the NIfTI image centre (as defined in the header) to the CoM
# of the extracted brain.
from sammba.interfaces import RatsMM

clip_level = memory.cache(afni.ClipLevel)
rats = memory.cache(RatsMM)
apply_mask = memory.cache(fsl.ApplyMask)
center_mass = memory.cache(afni.CenterMass)
brain_files = []
for unifized_file in unifized_files:
    out_clip_level = clip_level(in_file=unifized_file)
    out_rats = rats(in_file=unifized_file,
                    out_file='UnBm.nii.gz',
                    volume_threshold=400,
                    intensity_threshold=int(out_clip_level.outputs.clip_val))
    out_apply_mask = apply_mask(in_file=unifized_file,
                                mask_file=out_rats.outputs.out_file,
                                out_file='UnBmBe.nii.gz')
    out_center_mass = center_mass(in_file=out_apply_mask.outputs.out_file,
                                  cm_file='/tmp/cm.txt',
                                  set_cm=(0, 0, 0))
    brain_files.append(out_center_mass.outputs.out_file)

##############################################################################
# Same header change, for head files.
refit = memory.cache(afni.Refit)
head_files = []
for brain_file in brain_files:
    out_refit = refit(in_file=unifized_file, duporigin_file=brain_file)
    head_files.append(out_refit.outputs.out_file)

##############################################################################
# Quality check videos and average head and brain.
tcat = memory.cache(afni.TCat)
tstat = memory.cache(afni.TStat)
out_tcat = tcat(in_files=head_files, out_file='Un_video.nii.gz')
_ = tstat(in_file=out_tcat.outputs.out_file, out_file='Un_mean.nii.gz')
out_tcat = tcat(in_files=brain_files, out_file='UnBmBe_video.nii.gz')
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
# Finally, shift heads and brains within the images to place the CoM at the
# image center.
resample = memory.cache(afni.Resample)
shifted_head_files = []
for head_file in head_files:
    out_resample = resample(in_file=head_file,
                            master=out_refit.outputs.out_file,
                            out_file='UnCC.nii.gz')
    shifted_head_files.append(out_resample.outputs.out_file)

shifted_brain_files = []
for brain_file in brain_files:
    out_resample = resample(in_file=brain_file,
                            master=out_refit.outputs.out_file,
                            out_file='UnBmBeCC.nii.gz')
    shifted_brain_files.append(out_resample.outputs.out_file)

##############################################################################
# Quality check videos and template
out_tcat = tcat(in_files=brain_files, out_file='UnCC_video.nii.gz')
_ = tstat(in_file=out_tcat.outputs.out_file, out_file='UnCC_mean.nii.gz')
out_tcat = tcat(in_files=shifted_brain_files, out_file='UnBmBeCC_video.nii.gz')
out_tstat_shifted_brain = tstat(in_file=out_tcat.outputs.out_file,
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
shift_rotated_brain_files = []
rigid_transform_files = []
for shifted_brain_file in shifted_brain_files:
    out_allineate = allineate(in_file=shifted_brain_file,
                              reference=out_tstat_shifted_brain.outputs.out_file,
                              out_matrix='UnBmBeCCAl3.aff12.1D',
                              convergence=0.005,
                              two_blur=1,
                              warp_type='shift_rotate',
                              out_file='UnBmBeCCAl3.nii.gz')
    rigid_transform_files.append(out_allineate.outputs.out_matrix)
    shift_rotated_brain_files.append(out_allineate.outputs.out_file)

##############################################################################
# Application to the whole head image. can also be used for a good
# demonstration of linear vs. non-linear registration quality
shift_rotated_head_files = []
for shifted_head_file, transform_matrix in zip(shifted_head_files,
                                               rigid_transform_files):
    out_allineate = allineate(in_file=shifted_head_file,
                              master=out_tstat_shifted_brain.outputs.out_file,
                              in_matrix=transform_matrix,
                              out_file='UnBmBeCCAa3.nii.gz')
    shift_rotated_head_files.append(out_allineate.outputs.out_file)

##############################################################################
# Quality check videos and template
out_tcat = tcat(in_files=shift_rotated_brain_files,
                out_file='shr_video.nii.gz')
_ = tstat(in_file=out_tcat.outputs.out_file, out_file='shr_mean.nii.gz')

##############################################################################
# Note that this rigid body registration may need to be run more than once.
out_tcat = tcat(in_files=shift_rotated_head_files,
                out_file='shr_videohead.nii.gz')
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
# We begin by achieving an affine registration on aligned heads.
# A weighting mask is used to ...
import os

affine_transform_files = []
catmatvec = memory.cache(afni.CatMatvec)
for shift_rotated_head_file, rigid_transform_file in zip(shift_rotated_head_files,
                                                         rigid_transform_files):
    out_allineate = allineate(in_file=shift_rotated_head_file,
                              reference=out_tstat_shr.outputs.out_file,
                              out_matrix='UnBmBeAl4.aff12.1D',
                              convergence=0.005,
                              two_blur=1,
                              one_pass=True,
                              weight=out_mask_tool.outputs.out_file,
                              out_file='UnCCAl4.nii.gz')
    catmatvec_out_file = os.path.join(os.path.dirname(rigid_transform_file),
                                      'UnBmBeCCAl3UnCCAl4.aff12.1D')
    out_catmatvec = catmatvec(in_file=[(rigid_transform_file, 'ONELINE'),
                                       (out_allineate.outputs.out_matrix,
                                        'ONELINE')],
                              out_file=catmatvec_out_file)
    affine_transform_files.append(catmatvec_out_file)

##############################################################################
# Each resulting registration matrix is concatenated to the corresponding
# rigid bory registration matrix then directly applied to the CoM brain
# and head, reducing reslice errors in the final result.
allineated_brain_files = []
for shifted_brain_file, affine_transform_file in zip(
        shifted_brain_files, affine_transform_files):
    out_allineate = allineate(in_file=shifted_brain_file,
                              master=out_tstat_shr.outputs.out_file,
                              in_matrix=affine_transform_file,
                              out_file='UnBmBeCCAa4.nii.gz')
    allineated_brain_files.append(out_allineate.outputs.out_file)

##############################################################################
# The application to the whole head image can also be used for a good
# demonstration of linear vs. non-linear registration quality.
allineated_head_files = []
for shifted_head_file, affine_transform_file in zip(shifted_brain_files,
                                                    affine_transform_files):
    out_allineate = allineate(in_file=shifted_head_file,
                              master=out_tstat_shr.outputs.out_file,
                              in_matrix=affine_transform_file,
                              out_file='UnCCAa.nii.gz')
    allineated_head_files.append(out_allineate.outputs.out_file)

##############################################################################
# Quality check videos and template
out_tcat = tcat(in_files=allineated_brain_files, out_file='aff4_video.nii.gz')
_ = tstat(in_file=out_tcat.outputs.out_file, out_file='aff4_mean.nii.gz')
out_tcat_head = tcat(in_files=allineated_head_files,
                     out_file='aff4_videohead.nii.gz')
out_tstat_allineated_head = tstat(in_file=out_tcat_head.outputs.out_file,
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
for affine_transform_file, shifted_head_file in zip(affine_transform_files,
                                                    shifted_head_files):
    out_qwarp = qwarp(in_file=shifted_head_file,
                      base_file=out_tstat_allineated_head.outputs.out_file,
                      nmi=True,
                      noneg=True,
                      iwarp=True,
                      weight=out_mask_tool.outputs.out_file,
                      iniwarp=[affine_transform_file],
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
for warp_file, shifted_head_file in zip(warp_files, shifted_head_files):
    out_nwarp_cat = nwarp_cat(
        in_files=[('IDENT', out_tstat_warp_head.outputs.out_file),
                  warp_file], out_file='iniwarp.nii.gz')
    out_qwarp = qwarp(in_file=shifted_head_file,
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

##############################################################################
# Using previous files and concatenated transforms can be exploited to avoid
# building up reslice errors.
# Warp with mini-patch
warped_files3 = []
warp_files3 = []
for warp_file, shifted_head_file in zip(warp_files2, shifted_head_files):
    out_qwarp = qwarp(in_file=shifted_head_file,
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

##############################################################################
# In this particular case, minpathc=13 corresponds to a level of 10

##############################################################################
# Final warp
warped_files4 = []
for warp_file, shifted_head_file in zip(warp_files3, shifted_head_files):
    out_qwarp = qwarp(in_file=shifted_head_file,
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
stop
##############################################################################
# Registration to template
# ------------------------
# Apply warp to anat and atlas
# Apply non-linear registration results to uncorrected images
warp_apply = memory.cache(afni.NwarpApply)
for head_file, warp_file in zip(head_files, warp_files):
    out_warp_apply = warp_apply(in_file=head_file,
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