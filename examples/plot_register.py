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

import nibabel
from nilearn.image.resampling import coord_transform
from scipy.ndimage import center_of_mass

out_apply_mask_img = nibabel.load(out_apply_mask.outputs.out_file)
i, j, k = center_of_mass(out_apply_mask_img.get_data() > 0)
x, y, z = np.array(coord_transform(i, j, k, out_apply_mask_img.affine)).T


ambmc_data = left_ambmc_img.get_data()
roi_label = ambmc_atlas.labels.id[ambmc_atlas.labels.name == region_name][0]
roi_mask_array = ambmc_data == ambmc_atlas.labels.id[roi_label]
roi_mask_img = image.new_img_like(left_ambmc_img, roi_mask_array)
i, j, k = center_of_mass(roi_mask_array)
x, y, z = np.array(coord_transform(i, j, k, left_ambmc_img.affine)).T


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
