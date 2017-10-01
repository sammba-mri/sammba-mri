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

memory = Memory('/tmp')
template_file = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz'
head_file = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm.nii.gz'

basenames = ['/tmp/test_anat/anat']
unifized_files = []
masked_files = []
brain_mask_files = []
# Loop on subjects
for basename in basenames:
    bias_correct = memory.cache(ants.N4BiasFieldCorrection)
    out_bias_correct = bias_correct(input_image=basename + '.nii.gz',
                                    output_image=basename + '_N4.nii.gz')

    unifize = memory.cache(afni.Unifize)
    out_unifize = unifize(in_file=out_bias_correct.outputs.output_image,
                          urad=18.3,
                          out_file=basename + '_Un.nii.gz')

    clip_level = memory.cache(afni.ClipLevel)
    out_clip_level = clip_level(in_file=out_unifize.outputs.out_file)

    #RATS_MM "$anat"_"$Bc".nii.gz "$anat"_"$Bc"Bm.nii.gz -v $brainvol -t $thresh
    bet = memory.cache(fsl.BET)
    out_bet = bet(in_file=out_unifize.outputs.out_file,
                  out_file=basename + '_UnBm.nii.gz')

    apply_mask = memory.cache(fsl.ApplyMask)
    out_apply_mask = apply_mask(in_file=out_unifize.outputs.out_file,
                                mask_file=out_bet.outputs.out_file,
                                out_file=basename + '_UnBmBe.nii.gz')

    center_mass = memory.cache(afni.CenterMass)
    out_center_mass = center_mass(in_file=out_apply_mask.outputs.out_file,
                                  cm_file='/tmp/cm.txt',
                                  set_cm=(0, 0, 0))

    refit = memory.cache(afni.Refit)
    out_refit = refit(in_file=basename + '.nii.gz',
                      duporigin_file=out_center_mass.outputs.out_file)
    unifized_files.append(out_refit.outputs.out_file)
    out_refit = refit(in_file=basename + '_Un.nii.gz',
                      duporigin_file=basename + '_UnBmBe.nii.gz')
    brain_mask_files.append(out_refit.outputs.out_file)
    out_refit = refit(in_file=basename + '_UnBm.nii.gz',
                      duporigin_file=basename + '_UnBmBe.nii.gz')
    masked_files.append(out_refit.outputs.out_file)

tcat = memory.cache(afni.TCat)
out_tcat = tcat(in_files=unifized_files, out_file='Un_video.nii.gz')
tstat = memory.cache(afni.TStat)
out_tstat = tstat(in_file=out_tcat.outputs.out_file, out_file='Un_mean.nii.gz')
out_tcat = tcat(in_files=masked_files, out_file='UnBmBe_video.nii.gz')
out_tstat = tstat(in_file=out_tcat.outputs.out_file,
                  out_file='UnBmBe_mean.nii.gz')

undump = memory.cache(afni.Undump)
out_undump = undump(in_file=out_tstat.outputs.out_file,
                    out_file='emptytemplate.nii.gz')
out_refit = refit(in_file=out_undump.outputs.out_file,
                  xorigin='cen', yorigin='cen', zorigin='cen')

resample = memory.cache(afni.Resample)
for anat_file in anat_files:
    out_resample = resample(in_file=anat_file, master='emptytemplate.nii.gz',
                            out_file=fname_presuffix(anat_file, suffix='CC'))

unifized_resampled = []
for unifized_file in unifized_files:
    out_resample = resample(in_file=unifized_file,
                            master='emptytemplate.nii.gz',
                            out_file=fname_presuffix(unifized_file,
                                                     suffix='CC'))
    unifized_resampled_files.append(out_resample.outputs.out_file)

for brain_mask_file in brain_mask_files:
    out_resample = resample(in_file=brain_mask_file,
                            master='emptytemplate.nii.gz',
                            out_file=fname_presuffix(brain_mask_file,
                                                     suffix='CC'))
for brain_mask_file in brain_mask_files:
    out_resample = resample(in_file=brain_mask_file,
                            master='emptytemplate.nii.gz',
                            out_file=fname_presuffix(brain_mask_file,
                                                     suffix='CC'))
masked_resampled_files = []
for masked_file in masked_files:
    out_resample = resample(in_file=masked_file, master='emptytemplate.nii.gz',
                            out_file=fname_presuffix(masked_file, suffix='CC'))
    masked_resampled_files.append(out_resample.outputs.out_file)

out_tcat = tcat(in_files=unifized_files, out_file='UnCC_video.nii.gz')
out_tstat = tstat(in_file=out_tcat.outputs.out_file,
                  out_file='UnCC_mean.nii.gz')
out_tcat = tcat(in_files=masked_files, out_file='UnBmBeCC_video.nii.gz')
out_tstat = tstat(in_file=out_tcat.outputs.out_file,
                  out_file='UnBmBeCC_mean.nii.gz')

for masked_resampled_file, unifized_resampled_file in zip(
        unifized_resampled_file, masked_resampled_files):
    allineate = memory.cache(afni.Allineate)
    out_allineate = allineate(in_file=masked_resampled_file,
                              reference=template_file,
                              master=template_file,
                              out_matrix=basename + '_UnBmBeAl3',
                              cost='nmi',
                              convergence=.05,
                              two_pass=True,
                              two_blur=6.,
                              center_of_mass='',
                              maxrot=90,
                              out_file=fname_presuffix(masked_resampled_file,
                                                       suffix='Al3'))
    # Application to the whole head image. can also be used for a good
    # demonstration of linear vs. non-linear registration quality
    out_allineate2 = allineate(in_file=unifized_resampled_file,
                               master=head_file,
                               in_matrix=out_allineate.outputs.matrix,
                               out_file=basename + '_UnAa3.nii.gz')

# Application to the whole head image. can also be used for a good
# demonstration of linear vs. non-linear registration quality
out_allineate2 = allineate(in_file=out_unifize.outputs.out_file,
                           master=head_file,
                           in_matrix=out_allineate.outputs.matrix,
                           out_file=basename + '_UnAa3.nii.gz')

for dir in $(find -mindepth 1 -maxdepth 1 -type d); do
	
	3dAllineate -base $template -source $dir/UnBmBeCC.nii.gz -prefix $dir/UnBmBeCCAl$n.nii.gz -1Dmatrix_save $dir/UnBmBeCCAl$n.aff12.1D -conv $conv -twoblur $twoblur -warp shr
	3dAllineate -input $dir/UnCC.nii.gz -master $template -prefix $dir/UnCCAa$n.nii.gz -1Dmatrix_apply $dir/UnBmBeCCAl$n.aff12.1D
	
done

3dTcat -prefix $savedir/shr"$n"_video.nii.gz $(find $savedir -name UnBmBeCCAl"$n".nii.gz)
3dTstat -prefix $savedir/shr"$n"_mean.nii.gz $savedir/shr"$n"_video.nii.gz
3dTcat -prefix $savedir/shr"$n"_videohead.nii.gz $(find $savedir -name UnCCAa"$n".nii.gz)
3dTstat -prefix $savedir/shr"$n"_meanhead.nii.gz $savedir/shr"$n"_videohead.nii.gz
3dmask_tool -count -inputs $savedir/shr"$n"_video.nii.gz -prefix $savedir/shr"$n"_count.nii.gz



# The actual T1anat to template registration using the brain extracted image
# could do in one 3dQwarp step using allineate flags but will separate as
# 3dAllineate performs well on brain image, and 3dQwarp well on whole head
allineate = memory.cache(afni.Allineate)
out_allineate = allineate(in_file=out_apply_mask.outputs.out_file,
                          reference=template_file,
                          master=template_file,
                          out_matrix=basename + '_UnBmBeAl3',
                          cost='nmi',
                          convergence=.05,
                          two_pass=True,
                          two_blur=6.,
                          center_of_mass='',
                          maxrot=90,
                          out_file=basename + '_UnBmBeAl3.nii.gz')

save_inverted_affine(out_allineate.outputs.matrix,
                     basename + '_UnBmBeAl3_INV.aff12.1D')

# Application to the whole head image. can also be used for a good
# demonstration of linear vs. non-linear registration quality
out_allineate2 = allineate(in_file=out_unifize.outputs.out_file,
                           master=head_file,
                           in_matrix=out_allineate.outputs.matrix,
                           out_file=basename + '_UnAa3.nii.gz')

# Skipping all video related commands
out_allineate3 = allineate(in_file=out_allineate2.outputs.out_file,
                           master=head_file,
                           out_matrix=basename + '_UnAa4',
                           convergence=5,
                           two_pass=True,
                           two_blur=6.,
#                          weight=,  # Fix nipype bug: weight can be a str, not only an existing file
                           out_file=basename + '_UnAa4.nii.gz')

cat_mat_vec = memory.cache(afni.CatMatvec)
out_cat = cat_mat_vec(in_file=[(out_allineate.outputs.matrix,
                                out_allineate3.outputs.matrix)],
                      oneline=True,
                      out_file=basename + 'UnBmBeAl3UnCCAl4.aff12.1D',
                      fourxfour=False)


	cat_matvec -ONELINE $dir/UnBmBeCCAl"$n1".aff12.1D $dir/UnCCAl"$n2".aff12.1D > $dir/UnBmBeCCAl"$n1"UnCCAl"$n2".aff12.1D
	3dAllineate -input $dir/UnCC.nii.gz -master $template -prefix $dir/UnCCAa"$n2".nii.gz -1Dmatrix_apply $dir/UnBmBeCCAl"$n1"UnCCAl"$n2".aff12.1D
	3dAllineate -input $dir/UnBmBeCC.nii.gz -master $template -prefix $dir/UnBmBeCCAa"$n2".nii.gz -1Dmatrix_apply $dir/UnBmBeCCAl"$n1"UnCCAl"$n2".aff12.1D

# Non-linear registration of affine pre-registered whole head image to template
# Don't initiate straight from the original with an iniwarp due to weird errors
# (like it creating an Allin it then can't find)
qwarp = memory.cache(afni.Qwarp)
out_qwarp = qwarp(in_file=out_allineate2.outputs.out_file,
                  base_file=template_file,
                  nmi=True,
                  noneg=True,
                  iwarp=True,
                  blur=[0],
                  out_file=basename + '_UnAaQw.nii.gz')
