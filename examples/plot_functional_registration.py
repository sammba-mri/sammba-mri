"""
fMRI registration
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

retest = data_fetchers.fetch_zurich_test_retest(subjects=[1])

##############################################################################
# and define the caching repository
from nipype.caching import Memory

memory = Memory('/tmp')

##############################################################################
# Register using center of mass
# -----------------------------
# An initial coarse registration is done using brain centre of mass (CoM).
from nipype.interfaces import afni, fsl

##############################################################################
# First we loop through anatomical scans and correct intensities for bias.
import os
import numpy as np
import nibabel
from nipype.utils.filemanip import fname_presuffix

# TODO: t_r must be set by user
tshift = memory.cache(afni.TShift)
func_filename = retest.func[0]
tr = '1'
img = nibabel.load(func_filename)
n_slices = img.header.get_data_shape()[2]
out_tshift = tshift(in_file=func_filename,
                    outputtype='NIFTI_GZ',
                    tpattern='altplus',
                    tr=tr)
tshifted_filename = out_tshift.outputs.out_file

clip_level = memory.cache(afni.ClipLevel)
threshold = memory.cache(fsl.Threshold)
out_clip_level = clip_level(in_file=tshifted_filename)
out_threshold = threshold(in_file=tshifted_filename,
                          thresh=out_clip_level.outputs.clip_val)
thresholded_filename = out_threshold.outputs.out_file

# Register to the first time point
volreg = memory.cache(afni.Volreg)
oned_filename = fname_presuffix(thresholded_filename,
                                suffix='Vr.1Dfile.1D',
                                use_ext=False)
oned_matrix_filename = fname_presuffix(thresholded_filename,
                                       suffix='Vr.aff12.1D',
                                       use_ext=False)
out_volreg = volreg(in_file=thresholded_filename,
                    outputtype='NIFTI_GZ',
                    oned_file=oned_filename,
                    oned_matrix_save=oned_matrix_filename)  # XXX dfile not saved
registered_thresholded_filename = out_volreg.outputs.out_file
matrix_filename = out_volreg.outputs.oned_matrix_save

# Apply the registration to the whole head
allineate = memory.cache(afni.Allineate)
allineated_filename = fname_presuffix(tshifted_filename, suffix='Av')
out_allineate = allineate(in_file=tshifted_filename,
                          master=tshifted_filename,
                          in_matrix=matrix_filename,
                          out_file=allineated_filename)
allineated_filename = out_allineate.outputs.out_file

# 3dAllineate removes the obliquity. This is not a good way to readd it as
# removes motion correction info in the header if it were an AFNI file...as it
# happens it's NIfTI which does not store that so irrelevant!
copy_geom = memory.cache(fsl.CopyGeom)
out_copy_geom = copy_geom(dest_file=allineated_filename,
                          in_file=registered_thresholded_filename)
allineated_oblique_filename = out_copy_geom.outputs.out_file

tstat = memory.cache(afni.TStat)
# Create a (hopefully) nice mean image for use in the registration
#averaged_filename = fname_presuffix(allineated_filename, suffix='Av')
out_tstat = tstat(in_file=allineated_filename, args='-mean',
                  outputtype='NIFTI_GZ')
averaged_filename = out_tstat.outputs.out_file

from nipype.interfaces import ants

bias_correct = memory.cache(ants.N4BiasFieldCorrection)
out_bias_correct = bias_correct(input_image=averaged_filename, dimension=3)
bias_corrected_filename = out_bias_correct.outputs.output_image

# optional prior whole brain rigid body registration.
rats = memory.cache(RatsMM)
out_clip_level = clip_level(in_file=unifized_file)
out_rats = rats(in_file=unifized_file,
                out_file='UnBm.nii.gz',
                volume_threshold=400,
                intensity_threshold=int(out_clip_level.outputs.clip_val))

	3dcalc -a "$base".nii.gz -b "$base"_Bm.nii.gz -expr "a*b" -prefix "$base"_BmBe.nii.gz

apply_mask = memory.cache(fsl.ApplyMask)
center_mass = memory.cache(afni.CenterMass)
    out_apply_mask = apply_mask(in_file=unifized_file,
                                mask_file=out_rats.outputs.out_file)
    out_center_mass = center_mass(
        in_file=out_apply_mask.outputs.out_file,
        cm_file=os.path.join(cache_directory, 'cm.txt'),
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
# The brain files with new image center are concatenated to produce
# a quality check video
tcat = memory.cache(afni.TCat)
out_tcat = tcat(in_files=brain_files, outputtype='NIFTI_GZ')

##############################################################################
# and averaged
tstat = memory.cache(afni.TStat)
out_tstat = tstat(in_file=out_tcat.outputs.out_file, outputtype='NIFTI_GZ')

##############################################################################
# to create an empty template, with origin placed at CoM
undump = memory.cache(afni.Undump)
out_undump = undump(in_file=out_tstat.outputs.out_file, outputtype='NIFTI_GZ')
out_refit = refit(in_file=out_undump.outputs.out_file,
                  xorigin='cen', yorigin='cen', zorigin='cen')

# Finally, we shift heads and brains within the images to place the CoM at the
# image center.
resample = memory.cache(afni.Resample)
out_resample = resample(in_file=head_file,
                        master=out_refit.outputs.out_file,
                        outputtype='NIFTI_GZ')
shifted_head_files = out_resample.outputs.out_file

out_resample = resample(in_file=brain_file,
                        master=out_refit.outputs.out_file,
                        outputtype='NIFTI_GZ')
shifted_brain_files = out_resample.outputs.out_file

# Quality check videos and average brain
out_tcat = tcat(in_files=shifted_brain_files, outputtype='NIFTI_GZ')
out_tstat_shifted_brain = tstat(in_file=out_tcat.outputs.out_file,
                                outputtype='NIFTI_GZ')

# Shift rotate
allineate = memory.cache(afni.Allineate)
out_allineate = allineate(
    in_file=shifted_brain_file,
    reference=out_tstat_shifted_brain.outputs.out_file,
    out_matrix='UnBmBeCCAl3.aff12.1D',
    convergence=0.005,
    two_blur=1,
    warp_type='shift_rotate',
    out_file='UnBmBeCCAl3.nii.gz')
rigid_transform_files = out_allineate.outputs.out_matrix
shift_rotated_brain_files = out_allineate.outputs.out_file

# Application to the whole head image. can also be used for a good
out_allineate = allineate(in_file=shifted_head_file,
                          master=out_tstat_shifted_brain.outputs.out_file,
                          in_matrix=rigid_transform_file,
                          out_file='UnBmBeCCAa3.nii.gz')
shift_rotated_head_files = out_allineate.outputs.out_file


	3dAllineate -base "$anatlab"BmBe.nii.gz -source "$base"_BmBe.nii.gz -prefix "$base"_BmBe_shr.nii.gz -1Dmatrix_save "$base"_BmBe_shr.aff12.1D -cmass -warp shr
	cat_matvec -ONELINE "$base"_BmBe_shr.aff12.1D -I > "$base"_BmBe_shr_INV.aff12.1D
	3dAllineate -input "$anatlab".nii.gz -master $base.nii.gz -prefix "$anatlab"_shr_"$base".nii.gz -1Dmatrix_apply "$base"_BmBe_shr_INV.aff12.1D
	3dAllineate -input "$atlaslab".nii.gz -master $base.nii.gz -prefix "$atlaslab"_shr_"$base".nii.gz -1Dmatrix_apply "$base"_BmBe_shr_INV.aff12.1D #not used later (yet)
	anatlab="$anatlab"_shr_"$base"
suppanatwarp="$base"_BmBe_shr.aff12.1D

# transform anatomical and atlas to functional space. atlas is already in
# anatomical space, so only need to record matrix once, from the anatomical
# XXX atlas not done
unifize = memory.cache(afni.Unifize)
out_unifize = unifize(in_file=retest.anat[0], outputtype='NIFTI_GZ')
anat_filename = out_unifize.outputs.out_file

warp = memory.cache(afni.Warp)
out_warp = warp(in_file=anat_filename,
                oblique_parent=bias_corrected_filename,
                interp='quintic',
                gridset=bias_corrected_filename,
                outputtype='NIFTI_GZ',
                verbose=True)
allineated_anat_filename = out_warp.outputs.out_file
mat_filename = '/tmp/test.mat'
np.savetxt(mat_filename, [out_warp.runtime.stdout], fmt='%s')

# Reformat transform matrix
catmatvec = memory.cache(afni.CatMatvec)
out_catmatvec = catmatvec(in_file=[(mat_filename, 'ONELINE')],
                          out_file=fname_presuffix(mat_filename,
                                                   suffix='.aff12.1D ',
                                                   use_ext=False))  # XXX not used

# 3dWarp doesn't put the obliquity in the header, so do it manually
# This step generates one file per slice and per time point, so we are making
# sure they are removed at the end
import shutil

copy = memory.cache(afni.Copy)
refit = memory.cache(afni.Refit)
tmp_folder = os.path.join(
    os.path.dirname(bias_corrected_filename), 'tmp')
if not os.path.isdir(tmp_folder):
    os.makedirs(tmp_folder)

bias_corrected_basename = os.path.basename(bias_corrected_filename)
orig_bias_corrected_filename = fname_presuffix(os.path.join(
    tmp_folder, bias_corrected_basename), suffix='+orig.BRIK',
    use_ext=False)
out_copy_oblique = copy(in_file=bias_corrected_filename,
                        out_file=orig_bias_corrected_filename,
                        environ={'AFNI_DECONFLICT': 'OVERWRITE'})

allineated_anat_basename = os.path.basename(allineated_anat_filename)
orig_allineated_anat_filename = fname_presuffix(os.path.join(
    tmp_folder, allineated_anat_basename), suffix='+orig.BRIK',
    use_ext=False)
out_copy = copy(in_file=allineated_anat_filename,
                out_file=orig_allineated_anat_filename,
                environ={'AFNI_DECONFLICT': 'OVERWRITE'})

out_refit = refit(in_file=out_copy.outputs.out_file,
                  atrcopy=(out_copy_oblique.outputs.out_file,
                           'IJK_TO_DICOM_REAL'))
out_copy = copy(in_file=out_refit.outputs.out_file,
                environ={'AFNI_DECONFLICT': 'OVERWRITE'},
                out_file=allineated_anat_filename)
# shutil.rmtree(tmp_folder) # XXX to do later on

# slice anatomical
anat_img = nibabel.load(allineated_anat_filename)
anat_n_slices = anat_img.header.get_data_shape()[2]
slicer = memory.cache(afni.ZCutUp)
sliced_allineated_anat_filenames = []
for slice_n in range(anat_n_slices):
    out_slicer = slicer(in_file=allineated_anat_filename,
                        keep='{0} {1}'.format(slice_n, slice_n),
                        out_file=fname_presuffix(out_copy.outputs.out_file,
                                                 suffix='Sl%d' % slice_n))
    sliced_allineated_anat_filenames.append(out_slicer.outputs.out_file)

sliced_bias_corrected_filenames = []
for slice_n in range(n_slices):
    out_slicer = slicer(in_file=bias_corrected_filename,
                        keep='{0} {1}'.format(slice_n, slice_n),
                        out_file=fname_presuffix(bias_corrected_filename,
                                                 suffix='Sl%d' % slice_n))
    sliced_bias_corrected_filenames.append(out_slicer.outputs.out_file)

# XXX have the impression that you register to the average file


# Below line is to deal with slices where there is no signal (for example
# rostral end of some anatomicals)
empty_slices = []
for n, (sliced_bias_corrected_filename, sliced_allineated_anat_filename) in enumerate(
        zip(sliced_bias_corrected_filenames, sliced_allineated_anat_filenames)):
    out_clip_anat = clip_level(in_file=sliced_allineated_anat_filename)
    out_clip_func = clip_level(in_file=sliced_bias_corrected_filename)
    if out_clip_anat.outputs.clip_val < 1e-7 or out_clip_func.outputs.clip_val < 1e-7:
#        sliced_bias_corrected_filenames.pop(n)
#        sliced_allineated_anat_filenames.pop(n)
        empty_slices.append(n)
    print(out_clip_anat.outputs.clip_val)
    print(out_clip_func.outputs.clip_val)


# Single slice non-linear functional to anatomical registration; the iwarp
# frequently fails. Resampling can help it work better

# Read voxel size along third dimension
voxel_size_z = anat_img.header.get_zooms()[2]
resample = memory.cache(afni.Resample)
resampled_allineated_anat_filenames = []
for sliced_allineated_anat_filename in sliced_allineated_anat_filenames:
    out_resample = resample(in_file=sliced_allineated_anat_filename,
                            voxel_size=(.1, .1, voxel_size_z),
                            outputtype='NIFTI_GZ')
    print(nibabel.load(sliced_allineated_anat_filename).get_data().max())
    print(nibabel.load(out_resample.outputs.out_file).get_data().max())
    resampled_allineated_anat_filenames.append(out_resample.outputs.out_file)

resampled_bias_corrected_filenames = []
for sliced_bias_corrected_filename in sliced_bias_corrected_filenames:
    out_resample = resample(in_file=sliced_bias_corrected_filename,
                            voxel_size=(.1, .1, voxel_size_z),
                            outputtype='NIFTI_GZ')
    print(nibabel.load(out_resample.outputs.out_file).get_data().max())
    resampled_bias_corrected_filenames.append(out_resample.outputs.out_file)

qwarp = memory.cache(afni.Qwarp)
warped_slices = []
for (resampled_bias_corrected_filename,
     resampled_allineated_anat_filename) in zip(
    resampled_bias_corrected_filenames, resampled_allineated_anat_filenames):
    print(nibabel.load(resampled_bias_corrected_filename).get_data().max())
    print(nibabel.load(resampled_allineated_anat_filename).get_data().max())
    try:
        out_qwarp = qwarp(in_file=resampled_bias_corrected_filename,
                      base_file=resampled_allineated_anat_filename,
                      iwarp=True,
                      noneg=True,
                      blur=[0],
                      nmi=True,
                      noXdis=True,
                      allineate=True,
                      allineate_opts='-parfix 1 0 -parfix 2 0 -parfix 3 0 '
                                     '-parfix 4 0 -parfix 5 0 -parfix 6 0 -parfix '
                                      '7 0 -parfix 9 0 -parfix 10 0 -parfix 12 0',
                      out_file=warped_slice)
        warped_slices.append(out_qwarp.out_put.warped_source)
    except RuntimeError:
        pass
voxel_size = nibabel(warped_slice).header.get_zooms()[:3]
out_resample = resample(in_file=sliced_bias_corrected_filename,
                        voxel_size=voxel_size_z,
                        out_file=sliced_bias_corrected_filename,
                        environ={'AFNI_DECONFLICT': 'OVERWRITE'})
