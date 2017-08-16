import numpy as np
from nipype.interfaces import afni, fsl, ants
from nipype.caching import Memory


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
basename = '/tmp/test_anat/anat'
template_file = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz'
head_file = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm.nii.gz'

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

# The actual T1anat to template registration using the brain extracted image
# could do in one 3dQwarp step using allineate flags but will separate as
# 3dAllineate performs well on brain image, and 3dQwarp well on whole head
allineate = memory.cache(afni.Allineate)
out_allineate = allineate(in_file=out_apply_mask.outputs.out_file,
                          reference=template_file,
                          master=template_file,
                          out_matrix=basename + '_UnBmBeAl',
                          cost='nmi',
                          convergence=.05,
                          two_pass=True,
                          two_blur=6.,
                          center_of_mass='',
                          maxrot=90,
                          out_file=basename + '_UnBmBeAl.nii.gz')

save_inverted_affine(out_allineate.outputs.matrix,
                     basename + '_UnBmBeAl_INV.aff12.1D')

# Application to the whole head image. can also be used for a good
# demonstration of linear vs. non-linear registration quality
out_allineate2 = allineate(in_file=out_unifize.outputs.out_file,
                           master=head_file,
                           in_matrix=out_allineate.outputs.matrix,
                           out_file=basename + '_UnAa.nii.gz')


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
