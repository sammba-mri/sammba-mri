import os
from nose import with_setup
from nose.tools import assert_true
import numpy as np
import nibabel
from nilearn.datasets.tests import test_utils as tst
from nilearn._utils.testing import assert_raises_regex
from nilearn._utils.niimg_conversions import _check_same_fov
from sammba import testing_data
from sammba.registration.template_registrator import TemplateRegistrator


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_segment():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    registrator = TemplateRegistrator(anat_file, 400, output_dir=tst.tmpdir,
                                      use_rats_tool=False, verbose=False)
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    _, brain_file = registrator.segment(anat_file)
    assert_true(os.path.isfile(brain_file))


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_fit_transform_anat():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    anat_img = nibabel.load(anat_file)
    template_file = os.path.join(tst.tmpdir, 'template.nii.gz')
    template_affine = .2 * np.eye(4)
    template_affine[0, 1] = .01
    template_affine[1, 0] = .01
    template_affine[3, 3] = 1
    template_data = anat_img.get_data()[1:]
    template_img = nibabel.Nifti1Image(template_data, template_affine)
    template_img.to_filename(template_file)
    registrator = TemplateRegistrator(template_file, 400, output_dir=tst.tmpdir,
                                      use_rats_tool=False, verbose=False,
                                      registration_kind='affine')
    assert_raises_regex(
        ValueError, 'has not been fitted',
        registrator.transform_anat_like, anat_file)

    registrator.fit_anat(anat_file)
    registered_anat_img = nibabel.load(registrator.registered_anat_)
    assert_true(_check_same_fov(registered_anat_img, template_img))

    # test transform_anat_like on an image with same orientation as the template
    anat_like_img = nibabel.Nifti1Image(np.zeros(anat_img.get_data().shape), anat_img.affine)
    anat_like_file = os.path.join(tst.tmpdir, 'anat_like.nii.gz')
    anat_like_img.to_filename(anat_like_file)
    transformed_file = registrator.transform_anat_like(anat_like_file)
    transformed_img = nibabel.load(transformed_file)
    assert_true(_check_same_fov(transformed_img, template_img))
    
    # Same tests, with nonlinear registration
    registrator = TemplateRegistrator(template_file, 400, output_dir=tst.tmpdir,
                                      use_rats_tool=False, verbose=False)
    assert_raises_regex(
        ValueError, 'has not been fitted',
        registrator.transform_anat_like, anat_file)

    registrator.fit_anat(anat_file)
    registered_anat_img = nibabel.load(registrator.registered_anat_)
    assert_true(_check_same_fov(registered_anat_img, template_img))

    anat_like_img = nibabel.Nifti1Image(np.zeros(anat_img.get_data().shape), anat_img.affine)
    anat_like_file = os.path.join(tst.tmpdir, 'anat_like.nii.gz')
    anat_like_img.to_filename(anat_like_file)
    transformed_file = registrator.transform_anat_like(anat_like_file)
    transformed_img = nibabel.load(transformed_file)
    assert_true(_check_same_fov(transformed_img, template_img))


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_fit_transform_modality():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')
    anat_img = nibabel.load(anat_file)
    template_file = os.path.join(tst.tmpdir, 'template.nii.gz')
    template_affine = .2 * np.eye(4)
    template_affine[0, 1] = .01
    template_affine[1, 0] = .01
    template_affine[3, 3] = 1
    template_data = anat_img.get_data()[1:]
    template_img = nibabel.Nifti1Image(template_data, template_affine)
    template_img.to_filename(template_file)

    registrator = TemplateRegistrator(template_file, 400, output_dir=tst.tmpdir,
                                      use_rats_tool=False, verbose=False,
                                      registration_kind='affine')
    registrator.fit_anat(anat_file)

    assert_raises_regex(
        ValueError, "Only 'func' and 'perf' ", registrator.fit_modality,
        func_file, 'diffusion')
    assert_raises_regex(
        ValueError, "'t_r' is needed for slice ", registrator.fit_modality,
        func_file, 'func')
    assert_raises_regex(
        ValueError, 'has not been func fitted',
        registrator.transform_modality_like, func_file, 'func')

    # test fit_modality with func
    registrator.fit_modality(func_file, 'func', slice_timing=False)
    registered_func_img = nibabel.load(registrator.registered_func_)
    np.testing.assert_array_almost_equal(registered_func_img.affine, template_img.affine)
    np.testing.assert_array_equal(registered_func_img.shape[:-1], template_img.shape)

    # test transform_modality
    func_img = nibabel.load(func_file)
    func_like_img = nibabel.Nifti1Image(np.zeros(func_img.get_data().shape[:-1]),
                                        func_img.affine)
    func_like_file = os.path.join(tst.tmpdir, 'func_like.nii.gz')
    func_like_img.to_filename(func_like_file)
    transformed_file = registrator.transform_modality_like(func_like_file, 'func')
    transformed_img = nibabel.load(transformed_file)
    assert_true(_check_same_fov(transformed_img, template_img))

    # test transform then inverse transform brings back to the original image
    inverse_transformed_file = registrator.inverse_transform_towards_modality(
        transformed_file, 'func')
    inverse_transformed_img = nibabel.load(inverse_transformed_file)
    assert_true(_check_same_fov(inverse_transformed_img, func_img))
    np.testing.assert_array_equal(inverse_transformed_img.get_data(),
                                  func_like_img.get_data())

    # test inverse transform then transform brings back to the original image 
    transformed_file2 = registrator.transform_modality_like(inverse_transformed_file, 'func')
    assert_true(_check_same_fov(nibabel.load(transformed_file2), transformed_img))

    # similar tests with perf
    assert_raises_regex(
        ValueError, 'has not been perf fitted',
        registrator.transform_modality_like, func_file, 'perf')
    m0_img = nibabel.Nifti1Image(func_img.get_data()[..., 0], func_img.affine)
    m0_file = os.path.join(tst.tmpdir, 'm0.nii.gz')
    m0_img.to_filename(m0_file)
    registrator.fit_modality(m0_file, 'perf')
    _check_same_fov(nibabel.load(registrator.registered_perf_), template_img)

    # test transform_modality
    m0_like_img = nibabel.Nifti1Image(np.zeros(m0_img.shape), m0_img.affine)
    m0_like_file = os.path.join(tst.tmpdir, 'm0_like.nii.gz')
    m0_like_img.to_filename(m0_like_file)
    transformed_file = registrator.transform_modality_like(m0_like_file, 'perf')
    _check_same_fov(nibabel.load(transformed_file), template_img)

    # Similarly with nonlinear registration
    registrator = TemplateRegistrator(template_file, 400, output_dir=tst.tmpdir,
                                      use_rats_tool=False, verbose=False)
    registrator.fit_anat(anat_file)
    # test fit_modality with func
    registrator.fit_modality(func_file, 'func', slice_timing=False)
    registered_func_img = nibabel.load(registrator.registered_func_)
    np.testing.assert_array_almost_equal(registered_func_img.affine, template_img.affine)
    np.testing.assert_array_equal(registered_func_img.shape[:-1], template_img.shape)

    # test transform_modality
    func_img = nibabel.load(func_file)
    func_like_img = nibabel.Nifti1Image(np.zeros(func_img.get_data().shape[:-1]),
                                        func_img.affine)
    func_like_file = os.path.join(tst.tmpdir, 'func_like.nii.gz')
    func_like_img.to_filename(func_like_file)
    transformed_file = registrator.transform_modality_like(func_like_file, 'func')
    transformed_img = nibabel.load(transformed_file)
    assert_true(_check_same_fov(nibabel.load(transformed_file), template_img))

    # test transform then inverse transform brings back to the original image
    inverse_transformed_file = registrator.inverse_transform_towards_modality(
        transformed_file, 'func')
    inverse_transformed_img = nibabel.load(inverse_transformed_file)
    assert_true(_check_same_fov(inverse_transformed_img, func_img))
    np.testing.assert_array_equal(inverse_transformed_img.get_data(),
                                  func_like_img.get_data())

    # test inverse transform then transform brings back to the original image 
    transformed_file2 = registrator.transform_modality_like(inverse_transformed_file, 'func')
    transformed_img2 = nibabel.load(transformed_file2)
    np.testing.assert_array_equal(transformed_img2.affine,
                                  transformed_img.affine)
    np.testing.assert_array_equal(transformed_img2.get_data(),
                                  transformed_img.get_data())

    # similar tests with perf
    assert_raises_regex(
        ValueError, 'has not been perf fitted',
        registrator.transform_modality_like, func_file, 'perf')
    m0_img = nibabel.Nifti1Image(func_img.get_data()[..., 0], func_img.affine)
    m0_file = os.path.join(tst.tmpdir, 'm0.nii.gz')
    m0_img.to_filename(m0_file)
    registrator.fit_modality(m0_file, 'perf')
    _check_same_fov(nibabel.load(registrator.registered_perf_), template_img)

    # test transform_modality
    m0_like_img = nibabel.Nifti1Image(np.zeros(m0_img.shape), m0_img.affine)
    m0_like_file = os.path.join(tst.tmpdir, 'm0_like.nii.gz')
    m0_like_img.to_filename(m0_like_file)
    transformed_file = registrator.transform_modality_like(m0_like_file, 'perf')
    _check_same_fov(nibabel.load(transformed_file), template_img)

    # TODO: tests when resampling
