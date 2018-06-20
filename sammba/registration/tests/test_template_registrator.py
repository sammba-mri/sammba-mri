import os
from nose import with_setup
from nose.tools import assert_true
import numpy as np
import nibabel
from nilearn.datasets.tests import test_utils as tst
from nilearn._utils.testing import assert_raises_regex
from sammba import testing_data
from sammba.registration.template_registrator import TemplateRegistrator


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_segment():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    registrator = TemplateRegistrator(template=anat_file,
                                      brain_volume=400, write_dir=tst.tmpdir,
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

    registrator = TemplateRegistrator(template=anat_file, write_dir=tst.tmpdir,
                                      use_rats_tool=False, verbose=False)
    assert_raises_regex(
        ValueError, 'has not been fitted. ', registrator.fit_modality,
        func_file, 'func')
    assert_raises_regex(
        ValueError, 'has not been perf fitted',
        registrator.transform_anat_like, anat_file)

    registrator.fit_anat(anat_file)
    transformed_anat_file = registrator.transform_modality(anat_file)
    anat_img = nibabel.load(anat_file)
    transformed_anat_img = nibabel.load(transformed_anat_file)
    np.testing.assert_array_equal(anat_img.affine, transformed_anat_img.affine)
    np.testing.assert_array_equal(anat_img.get_data(),
                                  transformed_anat_img.get_data())

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

    registrator.fit_modality(func_file, 'func', slice_timing=False)
    assert_raises_regex(
        ValueError, 'has not been perf fitted',
        registrator.transform_modality_like, func_file, 'perf')

    transformed_func_file = registrator.transform_modality_like(func_file, 'func')
    registered_func_img = nibabel.load(registrator.registered_func_)
    transformed_func_img = nibabel.load(transformed_func_file)
    np.testing.assert_array_equal(registered_func_img.affine,
                                  transformed_func_img.affine)
    np.testing.assert_array_equal(registered_func_img.get_data(),
                                  transformed_func_img.get_data())
    func_file2 = registrator.inverse_transform_towards_modality(
        transformed_func_file, 'func')
    inverse_transformed_func_img = nibabel.load(func_file2)
    func_img = nibabel.load(registrator.registered_func_)
    inverse_transformed_func_img = nibabel.load(registrator.registered_func_)
    np.testing.assert_array_equal(func_img.affine,
                                  inverse_transformed_func_img.affine)
    np.testing.assert_array_equal(func_img.get_data(),
                                  inverse_transformed_func_img.get_data())
