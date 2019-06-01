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


def crop_and_oblique(in_file, out_file):
    img = nibabel.load(in_file)
    oblique_affine = .2 * np.eye(4)
    oblique_affine[0, 1] = .01
    oblique_affine[1, 0] = .01
    oblique_affine[3, 3] = 1
    oblique_data = img.get_data()[1:]
    oblique_img = nibabel.Nifti1Image(oblique_data, oblique_affine)
    oblique_img.to_filename(out_file)


def empty_img_like(in_file, out_file):
    img = nibabel.load(in_file)
    new_img = nibabel.Nifti1Image(np.zeros(img.get_data().shape),
                                        img.affine)
    new_img.to_filename(out_file)


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
def test_fit_anat_and_transform_anat_like():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    template_file = os.path.join(tst.tmpdir, 'template.nii.gz')
    # Create template
    crop_and_oblique(anat_file, template_file)
    registrator = TemplateRegistrator(template_file, 400,
                                      output_dir=tst.tmpdir,
                                      use_rats_tool=False, verbose=False,
                                      registration_kind='affine')
    assert_raises_regex(
        ValueError, 'has not been anat fitted',
        registrator.transform_anat_like, anat_file)

    # test fit_anat
    registrator.fit_anat(anat_file)
    assert_true(_check_same_fov(nibabel.load(registrator.registered_anat_),
                                nibabel.load(template_file)))

    # test transform_anat_like
    anat_like_file = os.path.join(tst.tmpdir, 'anat_like.nii.gz')
    empty_img_like(anat_file, anat_like_file)
    registrator.fit_anat(anat_file)
    
    transformed_file = registrator.transform_anat_like(anat_like_file)
    assert_true(_check_same_fov(nibabel.load(transformed_file),
                                nibabel.load(template_file)))


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_fit_transform_and_inverse_modality_with_func():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')
    template_file = os.path.join(tst.tmpdir, 'template.nii.gz')
    crop_and_oblique(anat_file, template_file)

    registrator = TemplateRegistrator(template_file, 400,
                                      output_dir=tst.tmpdir,
                                      use_rats_tool=False, verbose=False,
                                      registration_kind='affine')
    registrator.fit_anat(anat_file)

    assert_raises_regex(
        ValueError, "Only 'func' and 'perf' ", registrator.fit_modality,
        func_file, 'diffusion')
    assert_raises_regex(
        ValueError, "'t_r' is needed for slice ", registrator.fit_modality,
        func_file, 'func', reorient_only=True)
    assert_raises_regex(
        ValueError, 'has not been func fitted',
        registrator.transform_modality_like, func_file, 'func')

    # test fit_modality for func
    registrator.fit_modality(func_file, 'func', slice_timing=False,
                             reorient_only=True)
    registered_func_img = nibabel.load(registrator.registered_func_)
    template_img = nibabel.load(template_file)
    np.testing.assert_array_almost_equal(registered_func_img.affine,
                                         template_img.affine)
    np.testing.assert_array_equal(registered_func_img.shape[:-1],
                                  template_img.shape)


    # test transform_modality for func
    func_like_file = os.path.join(tst.tmpdir, 'func_like.nii.gz')
    empty_img_like(func_file, func_like_file)
    transformed_file = registrator.transform_modality_like(func_like_file,
                                                           'func')
    transformed_img = nibabel.load(transformed_file)
    assert_true(_check_same_fov(transformed_img, nibabel.load(template_file)))

    # test transform then inverse transform brings back to the original image
    inverse_transformed_file = registrator.inverse_transform_towards_modality(
        transformed_file, 'func')
    inverse_transformed_img = nibabel.load(inverse_transformed_file)
    func_like_img = nibabel.load(func_like_file)
    assert_true(_check_same_fov(inverse_transformed_img, func_like_img))
    np.testing.assert_array_equal(inverse_transformed_img.get_data(),
                                  func_like_img.get_data())

    # test inverse transform then transform brings back to the original image 
    transformed_file2 = registrator.transform_modality_like(
        inverse_transformed_file, 'func')
    transformed_img2 = nibabel.load(transformed_file2)
    assert_true(_check_same_fov(transformed_img2,
                                transformed_img))
    np.testing.assert_array_equal(transformed_img2.get_data(),
                                  transformed_img.get_data())


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_fit_and_transform_modality_with_perf():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')
    template_file = os.path.join(tst.tmpdir, 'template.nii.gz')
    crop_and_oblique(anat_file, template_file)
    registrator = TemplateRegistrator(template_file, 400,
                                      output_dir=tst.tmpdir,
                                      use_rats_tool=False, verbose=False,
                                      registration_kind='affine')
    registrator.fit_anat(anat_file)
    assert_raises_regex(
        ValueError, 'has not been perf fitted',
        registrator.transform_modality_like, func_file, 'perf')
    func_img = nibabel.load(func_file)
    m0_img = nibabel.Nifti1Image(func_img.get_data()[..., 0], func_img.affine)
    m0_file = os.path.join(tst.tmpdir, 'm0.nii.gz')
    m0_img.to_filename(m0_file)

    # test fit_modality for perf
    registrator.fit_modality(m0_file, 'perf', reorient_only=True)
    assert_true(_check_same_fov(nibabel.load(registrator.registered_perf_),
                                nibabel.load(template_file)))

    # test transform_modality for perf    
    m0_like_file = os.path.join(tst.tmpdir, 'm0_like.nii.gz')
    empty_img_like(m0_file, m0_like_file)
    transformed_file = registrator.transform_modality_like(m0_like_file,
                                                           'perf')
    assert_true(_check_same_fov(nibabel.load(transformed_file),
                                nibabel.load(template_file)))
