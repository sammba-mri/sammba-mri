import os
from nose.tools import assert_true, assert_equal, assert_less
from nose import with_setup
import nibabel
from nilearn.datasets.tests import test_utils as tst
from nilearn._utils.testing import assert_raises_regex
from nilearn._utils.niimg_conversions import _check_same_fov
from nilearn.image import mean_img
from sammba.registration import FMRISession, func
from sammba.registration.utils import _check_same_obliquity
from sammba import testing_data


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_coregister():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')
    mean_func_file = os.path.join(tst.tmpdir, 'mean_func.nii.gz')
    mean_img(func_file).to_filename(mean_func_file)

    bunch = func.coregister(anat_file, mean_func_file, tst.tmpdir,
                            slice_timing=False, verbose=False)
    assert_true(_check_same_fov(nibabel.load(bunch.coreg_func_),
                                nibabel.load(bunch.coreg_anat_)))
    assert_true(_check_same_obliquity(bunch.coreg_anat_,
                                      bunch.coreg_func_))
    assert_true(os.path.isfile(bunch.coreg_transform_))
    assert_less(0, len(bunch.coreg_warps_))
    assert_true(bunch.coreg_warps_[-1] is None)  # Last slice in functional
                                                 # is without signal
    for warp_file in bunch.coreg_warps_[:-1]:
        assert_true(os.path.isfile(warp_file))

    # Check environement variables setting
    assert_raises_regex(RuntimeError,
                        "3dcopy",
                        func.coregister, anat_file, mean_func_file, tst.tmpdir,
                        slice_timing=False, verbose=False,
                        AFNI_DECONFLICT='NO')

    # Check caching does not change the paths
    bunch2 = func.coregister(anat_file, mean_func_file, tst.tmpdir,
                             slice_timing=False, verbose=False, caching=True,
                             AFNI_DECONFLICT='OVERWRITE')
    assert_equal(bunch.coreg_func_, bunch2.coreg_func_)
    assert_equal(bunch.coreg_anat_, bunch2.coreg_anat_)
    assert_equal(bunch.coreg_transform_, bunch2.coreg_transform_)
    for warp_file, warp_file2 in zip(bunch.coreg_warps_, bunch2.coreg_warps_):
        assert_equal(warp_file, warp_file2)


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_fmri_sessions_to_template():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')
    mammal_data = FMRISession(anat=anat_file, func=func_file)

    template_file = anat_file
    t_r = 1.
    brain_volume = 400

    registered_data = func.fmri_sessions_to_template([mammal_data], t_r,
                                                     template_file,
                                                     tst.tmpdir,
                                                     brain_volume,
                                                     slice_timing=False,
                                                     verbose=True,
                                                     use_rats_tool=False)
    assert_true(os.path.isdir(registered_data[0].output_dir_))
    assert_true(os.path.isfile(registered_data[0].registered_func_))
    assert_true(os.path.isfile(registered_data[0].registered_anat_))

    assert_raises_regex(ValueError,
                        "'animals_data' input argument must be an iterable",
                        func.fmri_sessions_to_template, mammal_data, t_r,
                        template_file, tst.tmpdir, brain_volume)

    assert_raises_regex(ValueError,
                        "Each animal data must have type",
                        func.fmri_sessions_to_template, [mammal_data, ''], t_r,
                        template_file, tst.tmpdir, brain_volume)

    assert_raises_regex(ValueError,
                        "Animals ids must be different",
                        func.fmri_sessions_to_template,
                        [mammal_data, mammal_data],
                        t_r, template_file, tst.tmpdir, brain_volume)
