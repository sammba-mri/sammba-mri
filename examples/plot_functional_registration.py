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
t_r = 1.
for func_filename in :
    tpattern_filename = fname_presuffix(func_filename, suffix='_tpattern.txt',
                                        use_ext=False)
    if not os.path.isfile(tpattern_filename):
        img = nibabel.load(func_filename)
        n_slices = img.header.get_data_shape()[2]
        interslicetime = t_r / n_slices
        evenslices = np.arange(0, n_slices - 1, 2) * interslicetime
        oddslices = np.arange(1, n_slices, 2) * interslicetime
        slices = np.hstack((evenslices, oddslices))
        np.savetxt(tpattern_filename, np.atleast_2d(slices), delimiter=' ',
                   fmt='%f')

    tshifted_filename = fname_presuffix(func_filename, suffix='_Ts')
    out_tshift = tshift(in_file=func_filename,
                        out_file=tshifted_filename,
                        tpattern='@' + tpattern_filename)
    tshifted_filenames.append(out_tshift.outputs.out_file)

clip_level = memory.cache(afni.ClipLevel)
threshold = memory.cache(fsl.Threshold)
thresholded_filenames = []
for tshifted_filename in tshifted_filenames:
    out_clip_level = clip_level(in_file=tshifted_filename)
    out_threshold = threshold(in_file=tshifted_filename,
                              thresh=out_clip_level.outputs.clip_val)
    thresholded_filenames.append(out_threshold.outputs.out_file)

# Register to the first time point
volreg = memory.cache(afni.Volreg)
matrices_filenames = []
registered_thresholded_filenames = []
for thresholded_filename in thresholded_filenames:
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
    registered_thresholded_filenames.append(out_volreg.outputs.out_file)
    matrices_filenames.append(out_volreg.outputs.oned_matrix_save)

# Apply the registration to the whole head
allineate = memory.cache(afni.Allineate)
allineated_filenames = []
for tshifted_filename, matrix_filename in zip(tshifted_filenames,
                                              matrices_filenames):
    allineated_filename = fname_presuffix(tshifted_filename, suffix='Av')
    out_allineate = allineate(in_file=tshifted_filename,
                              master=tshifted_filename,
                              in_matrix=matrix_filename,
                              out_file=allineated_filename)
    allineated_filenames.append(out_allineate.outputs.out_file)

# 3dAllineate removes the obliquity. This is not a good way to readd it as
# removes motion correction info in the header if it were an AFNI file...as it
# happens it's NIfTI which does not store that so irrelevant!
copy_geom = memory.cache(fsl.CopyGeom)
allineated_oblique_filenames = []
for allineated_filename, registered_thresholded_filename in zip(
        allineated_filenames, registered_thresholded_filenames):
    out_copy_geom = copy_geom(dest_file=allineated_filename,
                              in_file=registered_thresholded_filename)
    allineated_oblique_filenames.append(out_copy_geom.outputs.out_file)

tstat = memory.cache(afni.TStat)
# Create a (hopefully) nice mean image for use in the registration
averaged_filenames = []
for allineated_oblique_filename in allineated_oblique_filenames:
    averaged_filename = fname_presuffix(allineated_oblique_filename,
                                        suffix='Av')
    out_tstat = tstat(in_file=allineated_filename, args='-mean',
                      out_file=averaged_filename)
    averaged_filenames.append(averaged_filename)

# N4 fails for some reason. Not tried 3dUnifize yet.
# XXX says 0 spacing not allowed. Replaced by Unifize
bias_correct = memory.cache(afni.Unifize)
bias_corrected_filenames = []
for averaged_filename in averaged_filenames:
    out_bias_correct = bias_correct(in_file=averaged_filename,
                                    outputtype='NIFTI_GZ')
    bias_corrected_filenames.append(out_bias_correct.outputs.out_file)

# optional prior whole brain rigid body registration.
# transform anatomical and atlas to functional space. atlas is already in
# anatomical space, so only need to record matrix once, from the anatomical
# XXX atlas not done
anat_filenames = [fname_presuffix(a, suffix='Un') for a in retest.anat]
warp = memory.cache(afni.Warp)
catmatvec = memory.cache(afni.CatMatvec)
allineated_anat_filenames = []
for anat_filename, bias_corrected_filename in zip(anat_filenames,
                                                  bias_corrected_filenames):
    out_warp = warp(in_file=anat_filename,
                    oblique_parent=bias_corrected_filename,
                    interp='quintic',
                    gridset=bias_corrected_filename,
                    outputtype='NIFTI_GZ',
                    verbose=True)  # 'allatonce'
    allineated_anat_filenames.append(out_warp.outputs.out_file)
    mat_filename = '/tmp/test.mat'
    np.savetxt(mat_filename, [out_warp.runtime.stdout], fmt='%s')

    # Reformat transform matrix
    out_catmatvec = catmatvec(in_file=[(mat_filename, 'ONELINE')],
                              out_file=fname_presuffix(mat_filename,
                                                       suffix='.aff12.1D ',
                                                       use_ext=False))

# 3dWarp doesn't put the obliquity in the header, so do it manually
# This step generates one file per slice and per time point, so we are making
# sure they are removed at the end
import shutil

copy = memory.cache(afni.Copy)
refit = memory.cache(afni.Refit)
for allineated_anat_filename, bias_corrected_filename in zip(
        allineated_anat_filenames, bias_corrected_filenames):
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
    shutil.rmtree(tmp_folder)

# slice functional and  anatomical
fsl_slice = memory.cache(fsl.slice)

fsl5.0-fslslice "$anatlab"_Op_"$base".nii.gz "$anatlab"_Op_"$base"
fsl5.0-fslslice "$atlaslab"_Op_"$base".nii.gz "$atlaslab"_Op_"$base" #not used later (yet)


goodoblique=bias_corrected_filename
badoblique=$out_warp
goodobliquename=$(basename -s .nii.gz $1)
badobliquename=$(basename -s .nii.gz $2)

3dcopy $goodoblique $workdir/"$goodobliquename"+orig
3dcopy $badoblique $workdir/"$badobliquename"+orig
3drefit -atrcopy $workdir/"$goodobliquename"+orig IJK_TO_DICOM_REAL $workdir/"$badobliquename"+orig
3dcopy $workdir/"$badobliquename"+orig $badoblique


fixobliquity.bash $base.nii.gz "$anatlab"_Op_"$base".nii.gz


br=$(basename -s .nii.gz $file)

			imfilearray=("$br"_TsAvAvN3 "$br"_TsAv)
			echo ${imfilearray[@]} > "$p"_imfilearray.txt
			bash $subpipeline $procdir $(basename -s AaQw.nii.gz ${a_Qwanats[0]}) "$p"_imfilearray.txt atlas_Na1 no $template $rminterregfiles $tmpdir 0.1 $registerfunctional $brainvol


