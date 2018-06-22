import os
from nose import with_setup
from nose.tools import assert_true, assert_equal
import numpy as np
import nibabel
from nilearn.datasets.tests import test_utils as tst
from nilearn._utils.testing import assert_raises_regex
from nilearn._utils.niimg_conversions import _check_same_fov
from sammba import testing_data
from sammba.registration.coregistrator import Coregistrator


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_segment():
    registrator = Coregistrator(brain_volume=400, output_dir=tst.tmpdir,
                                use_rats_tool=False, verbose=False)
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    unifized_file, brain_file = registrator.segment(anat_file)
    assert_true(os.path.isfile(brain_file))


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_coregistrator():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')

    registrator = Coregistrator(output_dir=tst.tmpdir, use_rats_tool=False, verbose=False)
    assert_raises_regex(
        ValueError, 'has not been fitted. ', registrator.fit_modality,
        func_file, 'func')

    registrator.fit_anat(anat_file)
    assert_equal(registrator.anat_, anat_file)
    assert_true(os.path.isfile(registrator.anat_brain_))

    assert_raises_regex(
        ValueError, "Only 'func' and 'perf' ", registrator.fit_modality,
        func_file, 'diffusion')
    assert_raises_regex(
        ValueError, "'t_r' is needed for slice ", registrator.fit_modality,
        func_file, 'func')
    assert_raises_regex(
        ValueError, '`brain_volume` must be ', registrator.fit_modality,
        func_file, 'func', slice_timing=False,
        prior_rigid_body_registration=True)

    # without slice timing
    registrator.fit_modality(func_file, 'func', slice_timing=False)
    func_img = nibabel.load(func_file)
    assert_true(_check_same_fov(nibabel.load(registrator.undistorted_func_), func_file))
    assert_true(_check_same_fov(nibabel.load(registrator.anat_in_func_space_), func_file))

    # with slice timing
    registrator.fit_modality(func_file, 'func', t_r=1.)
    assert_raises_regex(
        ValueError, 'has not been perf fitted',
        registrator.transform_modality_like, func_file, 'perf')
    func_img = nibabel.load(func_file)
    _check_same_fov(nibabel.load(registrator.undistorted_func_), func_img)
    _check_same_fov(nibabel.load(registrator.anat_in_func_space_), func_img)


    m0_img = nibabel.Nifti1Image(func_img.get_data()[..., :-1], func_img.affine)
    m0_file = os.path.join(tst.tmpdir, 'm0.nii.gz')
    m0_img.to_filename(m0_file)
    registrator.fit_modality(m0_file, 'perf')
    assert_true(_check_same_fov(nibabel.load(registrator.undistorted_perf_), m0_img))
    assert_true(_check_same_fov(nibabel.load(registrator.anat_in_perf_space_), m0_img))

    # test transform_modality_like on an image with same orientation as the functional
    func_like_img = nibabel.Nifti1Image(np.zeros(func_img.get_data().shape[:-1]),
                                        func_img.affine)
    func_like_file = os.path.join(tst.tmpdir, 'func_like.nii.gz')
    func_like_img.to_filename(func_like_file)
    transformed_file = registrator.transform_modality_like(func_like_file, 'func')
    registered_func_img = nibabel.load(registrator.undistorted_func_)
    transformed_img = nibabel.load(transformed_file)
    np.testing.assert_array_equal(registered_func_img.affine,
                                  transformed_img.affine)
    np.testing.assert_array_equal(registered_func_img.get_data().shape[:-1],
                                  transformed_img.get_data().shape)
