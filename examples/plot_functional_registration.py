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


def _fix_obliquity(to_fix_filename, reference_filename, caching_dir=None,
                   clear_memory=False):
    if caching_dir is None:
        caching_dir = os.getcwd()

    memory = Memory(caching_dir)
    copy = memory.cache(afni.Copy)
    tmp_folder = os.path.join(caching_dir, 'tmp')
    if not os.path.isdir(tmp_folder):
        os.makedirs(tmp_folder)

    reference_basename = os.path.basename(reference_filename)
    orig_reference_filename = fname_presuffix(os.path.join(
        tmp_folder, reference_basename), suffix='+orig.BRIK',
        use_ext=False)
    out_copy_oblique = copy(in_file=reference_filename,
                            out_file=orig_reference_filename,
                            environ={'AFNI_DECONFLICT': 'OVERWRITE'})

    to_fix_basename = os.path.basename(to_fix_filename)
    orig_to_fix_filename = fname_presuffix(os.path.join(
        tmp_folder, to_fix_basename), suffix='+orig.BRIK',
        use_ext=False)
    out_copy = copy(in_file=to_fix_filename,
                    out_file=orig_to_fix_filename,
                    environ={'AFNI_DECONFLICT': 'OVERWRITE'})

    out_refit = refit(in_file=out_copy.outputs.out_file,
                      atrcopy=(out_copy_oblique.outputs.out_file,
                               'IJK_TO_DICOM_REAL'))
    out_copy = copy(in_file=out_refit.outputs.out_file,
                    environ={'AFNI_DECONFLICT': 'OVERWRITE'},
                    out_file=to_fix_filename)
    if clear_memory:
        shutil.rmtree(tmp_folder)  # XXX to do later on
        memory.clear_previous_run()

##############################################################################
# Retrieve the data
from sammba import data_fetchers

retest = data_fetchers.fetch_zurich_test_retest(subjects=[1])

##############################################################################
# and define the caching repository
from nipype.caching import Memory

cache_directory = '/tmp'
memory = Memory(cache_directory)

from nipype.interfaces import afni, fsl

##############################################################################
# Correct orientation
# -------------------
# Before all, we need to correct the header for missing/wron sforms and qforms
import shutil
import os
import numpy as np
import nibabel

refit = memory.cache(afni.Refit)
center_mass = memory.cache(afni.CenterMass)
catmatvec = memory.cache(afni.CatMatvec)

func_filename = '/tmp/func.nii.gz'  # TODO: change to subject specific names
anat_filename = '/tmp/anat.nii.gz'
if not os.path.isfile(func_filename) or not os.path.isfile(anat_filename):
    out_refit = refit(in_file=retest.func[0], xyzscale=.1)
    out_center_mass = center_mass(
        in_file=out_refit.outputs.out_file,
        cm_file=os.path.join(cache_directory, 'cm.txt'),
        set_cm=(0, 0, 0))

    # orientation correctors for use later
    # https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions
    # Euler_angles_.28_z-y.E2.80.99-x.E2.80.B3_intrinsic.29_.E2.86.92_Rotation_matrix
    for raw_filename, filename in zip([retest.anat[0],
                                       out_center_mass.outputs.out_file],
                                      [anat_filename, func_filename]):
        shutil.copy(raw_filename, filename)
        img = nibabel.load(filename)
        header = img.header
        sform, code = header.get_sform(coded=True)
        sform[0] *= -1
        sform[[1, 2]] = sform[[2, 1]]
        header.set_sform(sform)
        header.set_qform(sform, int(code), strip_shears=False)
        nibabel.Nifti1Image(img.get_data(), sform, header).to_filename(filename)

#XXX Flip left right not done
##############################################################################
# Correct slice timing
# --------------------
tshift = memory.cache(afni.TShift)
out_tshift = tshift(in_file=func_filename,
                    outputtype='NIFTI_GZ',
                    tpattern='altplus',
                    tr='1')
tshifted_filename = out_tshift.outputs.out_file

##############################################################################
# Register to the first volume
# ----------------------------
# XXX why do you need a thresholded image ?
from nipype.utils.filemanip import fname_presuffix

clip_level = memory.cache(afni.ClipLevel)
threshold = memory.cache(fsl.Threshold)
out_clip_level = clip_level(in_file=tshifted_filename)
out_threshold = threshold(in_file=tshifted_filename,
                          thresh=out_clip_level.outputs.clip_val)
thresholded_filename = out_threshold.outputs.out_file

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
# XXX: bad output: up and down on y-axis
##############################################################################
# Apply the registration to the whole head
allineate = memory.cache(afni.Allineate)
allineated_filename = fname_presuffix(tshifted_filename, suffix='Av')
out_allineate = allineate(in_file=tshifted_filename,
                          master=tshifted_filename,
                          in_matrix=out_volreg.outputs.oned_matrix_save,
                          out_file=allineated_filename)

##############################################################################
# 3dAllineate removes the obliquity. This is not a good way to readd it as
# removes motion correction info in the header if it were an AFNI file...as it
# happens it's NIfTI which does not store that so irrelevant!
copy_geom = memory.cache(fsl.CopyGeom)
out_copy_geom = copy_geom(dest_file=out_allineate.outputs.out_file,
                          in_file=out_volreg.outputs.out_file)
allineated_filename = out_copy_geom.outputs.out_file

# XXX: bad output: up and down on y-axis
##############################################################################
# Create a (hopefully) nice mean image for use in the registration
tstat = memory.cache(afni.TStat)
out_tstat = tstat(in_file=allineated_filename, args='-mean',
                  outputtype='NIFTI_GZ')
averaged_filename = out_tstat.outputs.out_file

##############################################################################
# Rigid body registration
# -----------------------
# In this step, we register the anatomical volume to the mean functional
# volume.
#
# First we correct the functional average for intensities bias
from nipype.interfaces import ants

bias_correct = memory.cache(ants.N4BiasFieldCorrection)
out_bias_correct = bias_correct(input_image=averaged_filename, dimension=3)
unbiased_func_filename = out_bias_correct.outputs.output_image

##############################################################################
# as well as antomical image
unifize = memory.cache(afni.Unifize)
out_unifize = unifize(in_file=anat_filename, outputtype='NIFTI_GZ')
unbiased_anat_filename = out_unifize.outputs.out_file

##############################################################################
# Then we mask the mean functional volume outside the brain.
from sammba.interfaces import RatsMM

rats = memory.cache(RatsMM)
calc = memory.cache(afni.Calc)
out_clip_level = clip_level(in_file=unbiased_func_filename)
# XXX bad: brain mask cut
out_rats = rats(in_file=unbiased_func_filename,
                out_file='UnBm.nii.gz',
                volume_threshold=400,
                intensity_threshold=int(out_clip_level.outputs.clip_val))
out_cacl = calc(in_file_a=unbiased_func_filename,
                in_file_b=out_rats.outputs.out_file,
                expr='a*b',
                outputtype='NIFTI_GZ')

##############################################################################
# We compute the transformation from the functional image to the anatomical
# XXX: why in this sense
allineate = memory.cache(afni.Allineate)
out_allineate = allineate(
    in_file=out_cacl.outputs.out_file,
    reference=unbiased_anat_filename,
    out_matrix='BmBe_shr.aff12.1D',
    center_of_mass='',
    warp_type='shift_rotate',
    out_file='BmBe_shr.nii.gz')
rigid_transform_file = out_allineate.outputs.out_matrix

##############################################################################
# and apply its inverse to register to the anatomical volume
catmatvec = memory.cache(afni.CatMatvec)
catmatvec_out_file = fname_presuffix(rigid_transform_file, suffix='INV')
if not os.path.isfile(catmatvec_out_file):
    out_catmatvec = catmatvec(in_file=[(rigid_transform_file, 'I')],
                              oneline=True,
                              out_file=catmatvec_out_file)
    # XXX not cached I don't understand why
out_allineate = allineate(
    in_file=unbiased_anat_filename,
    master=unbiased_func_filename,
    in_matrix=catmatvec_out_file,
    out_file='BmBe_shr_in_func_space.nii.gz')

#suppanatwarp="$base"_BmBe_shr.aff12.1D

##############################################################################
# Non-linear registration
# -----------------------
# XXX atlas not done
# XXX what is the difference between Warp and 3dQwarp?
warp = memory.cache(afni.Warp)
out_warp = warp(in_file=out_allineate.outputs.out_file,
                oblique_parent=unbiased_func_filename,
                interp='quintic',
                gridset=unbiased_func_filename,
                outputtype='NIFTI_GZ',
                verbose=True)
registered_anat_filename = out_warp.outputs.out_file
mat_filename = fname_presuffix(registered_anat_filename, suffix='_warp.mat',
                               use_ext=False)
if not os.path.isfile(mat_filename):
    np.savetxt(mat_filename, [out_warp.runtime.stdout], fmt='%s')

# Reformat transform matrix
catmatvec = memory.cache(afni.CatMatvec)
out_catmatvec = catmatvec(in_file=[(mat_filename, 'ONELINE')],
                          out_file=fname_presuffix(mat_filename,
                                                   suffix='.aff12.1D ',
                                                   use_ext=False))  # XXX not used

##############################################################################
# 3dWarp doesn't put the obliquity in the header, so do it manually
# This step generates one file per slice and per time point, so we are making
# sure they are removed at the end
_fix_obliquity(registered_anat_filename, unbiased_func_filename,
               caching_dir=os.path.dirname(unbiased_func_filename))

##############################################################################
# Per-slice registration
# -----------------------
# Here we achieve single slice non-linear functional to anatomical
# registration.
#
# Slice anatomical image
anat_img = nibabel.load(registered_anat_filename)
anat_n_slices = anat_img.header.get_data_shape()[2]
slicer = memory.cache(afni.ZCutUp)
sliced_registered_anat_filenames = []
for slice_n in range(anat_n_slices):
    out_slicer = slicer(in_file=registered_anat_filename,
                        keep='{0} {1}'.format(slice_n, slice_n),
                        out_file=fname_presuffix(registered_anat_filename,
                                                 suffix='Sl%d' % slice_n))
    sliced_registered_anat_filenames.append(out_slicer.outputs.out_file)

##############################################################################
# and mean functional
sliced_bias_corrected_filenames = []
img = nibabel.load(func_filename)
n_slices = img.header.get_data_shape()[2]
for slice_n in range(n_slices):
    out_slicer = slicer(in_file=unbiased_func_filename,
                        keep='{0} {1}'.format(slice_n, slice_n),
                        out_file=fname_presuffix(unbiased_func_filename,
                                                 suffix='Sl%d' % slice_n))
    sliced_bias_corrected_filenames.append(out_slicer.outputs.out_file)

# XXX Fix caching of slicer

##############################################################################
# Below line is to deal with slices where there is no signal (for example
# rostral end of some anatomicals)

##############################################################################
# The inverse warp frequently fails, Resampling can help it work better
voxel_size_z = anat_img.header.get_zooms()[2]
resample = memory.cache(afni.Resample)
resampled_registered_anat_filenames = []
for sliced_registered_anat_filename in sliced_registered_anat_filenames:
    out_resample = resample(in_file=sliced_registered_anat_filename,
                            voxel_size=(.1, .1, voxel_size_z),
                            outputtype='NIFTI_GZ')
    resampled_registered_anat_filenames.append(out_resample.outputs.out_file)

resampled_bias_corrected_filenames = []
for sliced_bias_corrected_filename in sliced_bias_corrected_filenames:
    out_resample = resample(in_file=sliced_bias_corrected_filename,
                            voxel_size=(.1, .1, voxel_size_z),
                            outputtype='NIFTI_GZ')
    resampled_bias_corrected_filenames.append(out_resample.outputs.out_file)

##############################################################################
# We are now ready for the non-linear warp
qwarp = memory.cache(afni.Qwarp)
warped_slices = []
warp_filenames = []
for (resampled_bias_corrected_filename,
     resampled_registered_anat_filename) in zip(
     resampled_bias_corrected_filenames, resampled_registered_anat_filenames):
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
                                     '-parfix 7 0 -parfix 9 0 -parfix 10 0 '
                                     '-parfix 12 0',
                      out_file=warped_slice)
    warped_slices.append(out_qwarp.outputs.warped_source)
    warp_filenames.append(out_qwarp.outputs.source_warp)

##############################################################################
# Finally, resample the mean volume back to the initial resolution,
voxel_size = nibabel.load(sliced_bias_corrected_filename).header.get_zooms()
for warped_slice in warped_slices:
    out_resample = resample(in_file=warped_slice,
                            voxel_size=voxel_size + (voxel_size[0],),
                            out_file=warped_slice,
                            environ={'AFNI_DECONFLICT': 'OVERWRITE'})

##############################################################################
# fix the obliquity,
for (resampled_registered_anat_filename, warped_slice) in zip(
        resampled_registered_anat_filenames, warped_slices):
    _fix_obliquity(warped_slice, resampled_registered_anat_filename)

###############################################################################
# and merge all slices !
merge = memory.cache(fsl.Merge)
out_merge = merge(in_files=warped_slices, dimension='z')

###############################################################################
# Apply warp
# ----------
# The precomputed warp can be applied to all functional volumes, after
# repeating the previous steps of slicing,
sliced_func_filenames = []
for slice_n in range(n_slices):
    out_slicer = slicer(in_file=allineated_filename,
                        keep='{0} {1}'.format(slice_n, slice_n),
                        out_file=fname_presuffix(allineated_filename,
                                                 suffix='Sl%d' % slice_n))
    sliced_func_filenames.append(out_slicer.outputs.out_file)

###############################################################################
# and resampling.
resampled_func_filenames = []
for sliced_func_filename in sliced_func_filenames:
    out_resample = resample(in_file=sliced_func_filename,
                            voxel_size=(.1, .1, voxel_size_z),
                            outputtype='NIFTI_GZ')
    resampled_func_filenames.append(out_resample.outputs.out_file)

###############################################################################
# The warp is applied slice by slice
warp_apply = memory.cache(afni.NwarpApply)
warped_func_slices = []
for (resampled_func_filename, warp_filename) in zip(
        resampled_func_filenames, warp_filenames):
    out_warp_apply = warp_apply(in_file=resampled_func_filename,
                                warp=warp_filename,
                                out_file=fname_presuffix(
                                    resampled_func_filename, suffix='_qw'))
    warped_func_slices.append(out_warp_apply.outputs.out_file)

###############################################################################
# and fix the obliquity
# XXX why no resampling back before ?
for (resampled_registered_anat_filename, warped_func_slice) in zip(
        resampled_registered_anat_filenames, warped_func_slices):
    _fix_obliquity(warped_func_slice, resampled_registered_anat_filename)

###############################################################################
# Finally, merge all slices !
out_merge_func = merge(in_files=warped_func_slices, dimension='z')
out_merge_anat = merge(in_files=resampled_registered_anat_filenames,
                       dimension='z')

###############################################################################
# Check out the result
from nilearn import plotting

display = plotting.plot_anat(out_merge_func.outputs.merged_file)
display.add_edges(out_merge_anat.outputs.merged_file)
plotting.show()
