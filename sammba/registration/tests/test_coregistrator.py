import os
from nose import with_setup
from nose.tools import assert_true
import numpy as np
import nibabel
from nilearn.datasets.tests import test_utils as tst
from nilearn._utils.testing import assert_raises_regex
from sammba import testing_data
from sammba.registration.base_registrator import Coregistrator


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_segment():
    registrator = Coregistrator(brain_volume=400, write_dir=tst.tmpdir)
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

    registrator = Coregistrator(write_dir=tst.tmpdir)
    assert_raises_regex(
        ValueError, 'has not been fitted. ', registrator.fit_modality,
        func_file, 'func')

    registrator.fit_anat(anat_file)
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

    transformed_func_file = registrator.transform_modality(func_file, 'func')
    registered_func_img = nibabel.load(registrator.registered_func_)
    transformed_func_img = nibabel.load(transformed_func_file)
    np.testing.assert_array_equal(registered_func_img.affine,
                                  transformed_func_img.affine)
    np.testing.assert_array_equal(registered_func_img.get_data(),
                                  transformed_func_img.get_data())
