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


def save_composed_affines(in_files, out_file):
    in_affines = []
    for in_affine_file in in_files:
        u11, u12, u13, v1, u21, u22, u23, v2, u31, u32, u33, v3 = np.loadtxt(
            in_affine_file)
        u = np.array([[u11, u12, u13], [u21, u22, u23], [u31, u32, u33]])
        v = np.array([[v1], [v2], [v3]])
        in_affines.append(np.hstack((u, v)))

    composed_affine = _compose_affines(in_affines)
    np.savetxt(out_file, np.hstack((composed_affine[0], composed_affine[1],
                                    composed_affine[2])))
    return out_file

# Retrieve the data
from sammba import data_fetchers

retest = data_fetchers.fetch_zurich_test_retest(subjects=range(3))

##############################################################################
# Register using center of mass
# -----------------------------
# An initial coarse registration is done using brain centre of mass (CoM).

conv=0.005
twoblur=1
brainvol=400
#-rbt values might be improveable
Urad=18.3
b=70
t=80

memory = Memory('/tmp')
template_file = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz'
head_file = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm.nii.gz'
weight_file = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain_mask_dil.nii.gz'
atlas_warp_file = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz'


unifize = memory.cache(afni.Unifize)
clip_level = memory.cache(afni.ClipLevel)
#RATS_MM "$anat"_"$Bc".nii.gz "$anat"_"$Bc"Bm.nii.gz -v $brainvol -t $thresh
bet = memory.cache(fsl.BET)
apply_mask = memory.cache(fsl.ApplyMask)
center_mass = memory.cache(afni.CenterMass)
refit = memory.cache(afni.Refit)

unifized_files = []
masked_files = []
# Loop through anatomical scans
for anat_filename in retest.anat:
    # Bias-correction. Note that rbt values might be improveable
    out_unifize = unifize(in_file=anat_filename,
                          rbt=(18.3, 70, 80),
                          out_file=fname_presuffix(anat_filename, suffix='Un'))

    # Brain extraction, aided by an approximate guessed brain volume
    out_clip_level = clip_level(in_file=out_unifize.outputs.out_file)
    out_bet = bet(in_file=out_unifize.outputs.out_file)
    out_apply_mask = apply_mask(in_file=out_unifize.outputs.out_file,
                                mask_file=out_bet.outputs.out_file)

    # Set the NIfTI image centre (as defined in the header) to the CoM
    # of the extracted brain
    out_center_mass = center_mass(in_file=out_apply_mask.outputs.out_file,
                                  cm_file='/tmp/cm.txt',
                                  set_cm=(0, 0, 0))
    masked_files.append(out_center_mass.outputs.out_file)
    out_refit = refit(in_file=out_unifize.outputs.out_file,
                      duporigin_file=out_center_mass.outputs.out_file)
    unifized_files.append(out_refit.outputs.out_file)

# Quality check videos
tcat = memory.cache(afni.TCat)
tstat = memory.cache(afni.TStat)
out_tcat = tcat(in_files=unifized_files, out_file='Un_video.nii.gz')
out_tstat = tstat(in_file=out_tcat.outputs.out_file, out_file='Un_mean.nii.gz')
out_tcat = tcat(in_files=masked_files, out_file='Un_masked_video.nii.gz')
out_tstat = tstat(in_file=out_tcat.outputs.out_file,
                  out_file='Un_masked_mean.nii.gz')

undump = memory.cache(afni.Undump)
out_undump = undump(in_file=out_tstat.outputs.out_file,
                    out_file='emptytemplate.nii.gz')
out_refit = refit(in_file=out_undump.outputs.out_file,
                  xorigin='cen', yorigin='cen', zorigin='cen')

#The head is then shifted within the image to place the CoM at the image center.

resample = memory.cache(afni.Resample)
resampled_files = []
for anat_file in retest.anat:
    out_resample = resample(in_file=anat_file,
                            master=out_refit.outputs.out_file,
                            out_file=fname_presuffix(anat_file, suffix='CC'))
    resampled_files.append(out_resample.outputs.out_file)

unifized_resampled_files = []
for unifized_file in unifized_files:
    out_resample = resample(in_file=unifized_file,
                            master=out_refit.outputs.out_file,
                            out_file=fname_presuffix(unifized_file,
                                                     suffix='CC'))
    unifized_resampled_files.append(out_resample.outputs.out_file)

for brain_mask_file in brain_mask_files:
    out_resample = resample(in_file=brain_mask_file,
                            master=out_refit.outputs.out_file,
                            out_file=fname_presuffix(brain_mask_file,
                                                     suffix='CC'))
for brain_mask_file in brain_mask_files:
    out_resample = resample(in_file=brain_mask_file,
                            master=out_refit.outputs.out_file,
                            out_file=fname_presuffix(brain_mask_file,
                                                     suffix='CC'))
masked_resampled_files = []
for masked_file in masked_files:
    out_resample = resample(in_file=masked_file,
                            master=out_refit.outputs.out_file,
                            out_file=fname_presuffix(masked_file, suffix='CC'))
    masked_resampled_files.append(out_resample.outputs.out_file)

out_tcat = tcat(in_files=unifized_files, out_file='UnCC_video.nii.gz')
out_tstat = tstat(in_file=out_tcat.outputs.out_file,
                  out_file='UnCC_mean.nii.gz')
out_tcat = tcat(in_files=masked_files, out_file='UnBmBeCC_video.nii.gz')
out_tstat = tstat(in_file=out_tcat.outputs.out_file,
                  out_file='UnBmBeCC_mean.nii.gz')

# At this point, we achieved a translation-only registration of the raw
# anatomical images to each other's brain's (as defined by the brain extractor)
# CoMs.

##############################################################################
# Shift rotation
allineate = memory.cache(afni.Allineate)
masked_allineated_files = []
unifized_allineated_files = []
out_matrices = []
for unifized_resampled_file, masked_resampled_file in zip(
        unifized_resampled_files, masked_resampled_files):
    out_allineate = allineate(in_file=masked_resampled_file,
                              reference=template_file,
                              out_matrix=basename + '_UnBmBeAl3.aff12.1D',
                              convergence=4,
                              two_blur=5,
                              warp_type='shift_rotate',
                              out_file=fname_presuffix(masked_resampled_file,
                                                       suffix='Al3'))
    out_matrices.append(out_allineate.outputs.out_matrix)
    masked_allineated_files.append(out_allineate.outputs.out_file)
    # Application to the whole head image. can also be used for a good
    # demonstration of linear vs. non-linear registration quality
    out_allineate2 = allineate(in_file=unifized_resampled_file,
                               master=template_file,
                               in_matrix=out_allineate.outputs.out_matrix,
                               out_file=fname_presuffix(masked_resampled_file,
                                                        suffix='Aa3'))
    unifized_allineated_files.append(out_allineate2.outputs.out_file)

out_tcat = tcat(in_files=masked_allineated_files, out_file='shr3_video.nii.gz')
out_tstat = tstat(in_file=out_tcat.outputs.out_file,
                  out_file='shr_mean.nii.gz')
out_tcat2 = tcat(in_files=unifized_allineated_files,
                 out_file='shr3_videohead.nii.gz')
out_tstat = tstat(in_file=out_tcat2.outputs.out_file,
                  out_file='shr3_meanhead.nii.gz')
mask_tool = memory.cache(afni.MaskTool)
out_mask_tool = mask_tool(in_file=out_tcat.outputs.out_file, count=True,
                          out_file='shr3_count.nii.gz')

##############################################################################
# Affine
masked_aa4 = []
unifized_aa4 = []
composed_files = []
for masked_allineated_file, unifized_allineated_file, out_matrix in zip(
        masked_allineated_files, unifized_allineated_files, out_matrices):
    out_allineate = allineate(in_file=unifized_allineated_file,
                              reference=template_file,
                              out_matrix=basename + '_UnBmBeAl4.aff12.1D',
                              convergence=5,
                              two_blur=6,
                              one_pass=True,
                              weight_file=weight_file,
                              out_file=fname_presuffix(unifized_resampled_file,
                                                       suffix='Al4'))
    composed_file = save_composed_affines(
        in_files=[out_matrix, out_allineate.outputs.out_matrix],
        out_file=fname_presuffix('baseline',
                                 suffix='UnBmBeCCAl3UnCCAl4.aff12.1D'))
    composed_files.append(composed_file)
    # Application to the whole head image. can also be used for a good
    # demonstration of linear vs. non-linear registration quality
    out_allineate = allineate(in_file=unifized_resampled_file,
                              master=template_file,
                              in_matrix=composed_file,
                              out_file=fname_presuffix(unifized_resampled_file,
                                                       suffix='Aa4'))
    unifized_aa4.append(out_allineate2.outputs.out_file)
    out_allineate2 = allineate(in_file=masked_resampled_file,
                               master=template_file,
                               in_matrix=composed_file,
                               out_file=fname_presuffix(masked_resampled_file,
                                                        suffix='Aa4'))
    masked_aa4.append(out_allineate2.outputs.out_file)


out_tcat = tcat(in_files=masked_aa4, out_file='aff4_video.nii.gz')
out_tstat = tstat(in_file=out_tcat.outputs.out_file,
                  out_file='aff4_mean.nii.gz')
out_tcat2 = tcat(in_files=unifized_aa4,
                 out_file='aff4_videohead.nii.gz')
out_tstat = tstat(in_file=out_tcat2.outputs.out_file,
                  out_file='aff4_meanhead.nii.gz')
out_mask_tool = mask_tool(in_file=out_tcat.outputs.out_file, count=True,
                          out_file='aff4_count.nii.gz')

##############################################################################
# Warp
nwarp_cat = memory.cache(afni.NwarpCat)
qwarp = memory.cache(afni.Qwarp)
warped_files = []
for composed_file, masked_aa4_file in zip(composed_files, masked_aa4):
    out_nwarp_cat = nwarp_cat(in_files=[('IDENT', template_file),
                                        composed_file],  # check with Nad
                              out_file='iniwarp.nii.gz')
    out_qwarp = qwarp(in_file=out_allineate2.outputs.out_file,
                      base_file=template_file,
                      nmi=True,
                      noneg=True,
                      iwarp=True,
                      weight=weight_file,
                      iniwarp=[out_nwarp_cat.outputs.out_file],
                      inilev=6,
                      maxlev=7,
                      out_file=basename + '_UnCCQw8.nii.gz')
    warped_files.append(out_qwarp.inputs['out_file'])

out_tcat = tcat(in_files=warped_files, out_file='Qw8_videohead.nii.gz')
out_tstat = tstat(in_file=out_tcat.outputs.out_file,
                  out_file='Qw8_meanhead.nii.gz')

# Warp with mini-patch
warped_files = []
warp_files = []
for composed_file, masked_aa4_file in zip(composed_files, masked_aa4):
    out_nwarp_cat = nwarp_cat(in_files=[('IDENT', template_file),
                                        composed_file],
                              out_file='iniwarp.nii.gz')
    out_qwarp = qwarp(in_file=out_allineate2.outputs.out_file,
                      base_file=template_file,
                      nmi=True,
                      noneg=True,
                      iwarp=True,
                      weight=weight_file,
                      iniwarp=[out_nwarp_cat.outputs.out_file],
                      inilev=6,
                      maxlev=7,
#                      minpatch=8,  # take s a lot of time
                      out_file=basename + '_UnCCQw8.nii.gz')
    warped_files.append(out_qwarp.outputs.warped_source)
    warp_files.append(out_qwarp.outputs.source_warp)

out_tcat = tcat(in_files=warped_files, out_file='Qw8_videohead.nii.gz')
out_tstat = tstat(in_file=out_tcat.outputs.out_file,
                  out_file='Qw8_meanhead.nii.gz')

##############################################################################
# Apply warp to anat and atlas
warp_apply = memory.cache(afni.NwarpApply)
for resampled_file, warp_file in zip(resampled_files, warp_files):
    out_warp_apply = warp_apply(in_file=resampled_file,
                                warp=warp_file,
                                master=out_tstat.outputs.out_file,
                                out_file='Na.nii.gz')
    out_nwarp_cat = nwarp_cat(in_files=[('INV', warp_file),
                                        atlas_warp_file],  # check with Nad
                              out_file='atlas_WARPINV.nii.gz')
    out_warp_apply = warp_apply(in_file=resampled_file,
                                warp=warp_file,
                                master=out_tstat.outputs.out_file,
                                ainterp='NN',
                                out_file='atlas_Na.nii.gz')