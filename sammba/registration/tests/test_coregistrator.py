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

    registrator = Coregistrator(output_dir=tst.tmpdir, use_rats_tool=False,
                                verbose=False)
    assert_raises_regex(
        ValueError, 'has not been fitted. ', registrator.fit_modality,
        func_file, 'func')

    registrator.fit_anat(anat_file)
    assert_equal(registrator.anat_, anat_file)
    assert_true(registrator._anat_brain_mask is None)

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

    # without rigid body registration
    registrator.fit_modality(func_file, 'func', slice_timing=False)
    func_img = nibabel.load(func_file)
    assert_true(_check_same_fov(nibabel.load(registrator.undistorted_func_),
                                func_img))
    np.testing.assert_array_almost_equal(
        nibabel.load(registrator.anat_in_func_space_).affine, func_img.affine)
    np.testing.assert_array_equal(
        nibabel.load(registrator.anat_in_func_space_).shape,
        func_img.shape[:-1])

    # test transform_modality_like on an image with oriented as the functional
    func_like_img = nibabel.Nifti1Image(np.zeros(func_img.shape[:-1]),
                                        func_img.affine)
    func_like_file = os.path.join(tst.tmpdir, 'func_like.nii.gz')
    func_like_img.to_filename(func_like_file)
    transformed_file = registrator.transform_modality_like(func_like_file,
                                                           'func')
    transformed_img = nibabel.load(transformed_file)
    np.testing.assert_array_almost_equal(transformed_img.affine,
                                         func_img.affine)
    np.testing.assert_array_equal(transformed_img.shape, func_img.shape[:-1])

    # Similarly with rigid body registration
    registrator.fit_modality(func_file, 'func', slice_timing=False,
                             prior_rigid_body_registration=True)
    func_img = nibabel.load(func_file)
    assert_true(_check_same_fov(nibabel.load(registrator.undistorted_func_),
                                func_img))
    np.testing.assert_array_almost_equal(
        nibabel.load(registrator.anat_in_func_space_).affine, func_img.affine)
    np.testing.assert_array_equal(
        nibabel.load(registrator.anat_in_func_space_).shape,
        func_img.shape[:-1])

    func_like_img = nibabel.Nifti1Image(np.zeros(func_img.shape[:-1]),
                                        func_img.affine)
    func_like_file = os.path.join(tst.tmpdir, 'func_like.nii.gz')
    func_like_img.to_filename(func_like_file)
    transformed_file = registrator.transform_modality_like(func_like_file,
                                                           'func')
    transformed_img = nibabel.load(transformed_file)
    np.testing.assert_array_almost_equal(transformed_img.affine,
                                         func_img.affine)
    np.testing.assert_array_equal(transformed_img.shape, func_img.shape[:-1])

    # Similarly with perf
    m0_img = nibabel.Nifti1Image(func_img.get_data()[..., 0], func_img.affine)
    m0_file = os.path.join(tst.tmpdir, 'm0.nii.gz')
    m0_img.to_filename(m0_file)
    registrator.fit_modality(m0_file, 'perf')
    assert_true(_check_same_fov(nibabel.load(registrator.undistorted_perf_),
                                m0_img))
    assert_true(_check_same_fov(nibabel.load(registrator.anat_in_perf_space_),
                                m0_img))

    m0_like_img = nibabel.Nifti1Image(np.zeros(m0_img.shape), m0_img.affine)
    m0_like_file = os.path.join(tst.tmpdir, 'm0_like.nii.gz')
    m0_like_img.to_filename(m0_like_file)
    transformed_file = registrator.transform_modality_like(m0_like_file,
                                                           'perf')
    transformed_img = nibabel.load(transformed_file)
    assert_true(_check_same_fov(transformed_img, m0_img))
