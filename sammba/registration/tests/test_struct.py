import os
import numpy as np
from nose.tools import assert_true
from nose import with_setup
from numpy.testing import assert_array_almost_equal
from nilearn._utils.testing import assert_raises_regex
from nilearn.image import math_img
from nilearn.datasets.tests import test_utils as tst
from sammba.registration import struct
from nilearn._utils.niimg_conversions import _check_same_fov
from sammba.registration.utils import _check_same_obliquity
from sammba import testing_data


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_anats_to_common():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    brain_mask_file = os.path.join(os.path.dirname(testing_data.__file__),
                                   'mask.nii.gz')

    # Check error is raised if wrong registration kind
    assert_raises_regex(ValueError, "Registration kind must be one of ",
                        struct.anats_to_common, [anat_file],
                        [brain_mask_file], tst.tmpdir,
                        registration_kind='rigidd')

    # test common space of one image is itself
    rigid = struct.anats_to_common(
        [anat_file], [brain_mask_file], tst.tmpdir, registration_kind='rigid',
        verbose=0)
    transform = np.loadtxt(rigid.transforms[0])
    assert_array_almost_equal(transform,
                              [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0])


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_anat_to_template():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    brain_mask_file = os.path.join(os.path.dirname(testing_data.__file__),
                              'mask.nii.gz')
    brain_file = os.path.join( tst.tmpdir, 'brain.nii.gz')
    math_img('img1 * img2', img1=anat_file, img2=brain_mask_file).to_filename(
        brain_file)

    # Check error is raised if wrong registration kind
    assert_raises_regex(ValueError, "Registration kind must be one of ",
                        struct.anat_to_template, anat_file, brain_file,
                        anat_file, brain_file,
                        tst.tmpdir, registration_kind='rigidd')

    # test common space of one image is itself
    register_result = struct.anat_to_template(anat_file,  brain_file,
                                              anat_file, brain_file, anat_file,
                                              tst.tmpdir,
                                              registration_kind='rigid',
                                              use_rats_tool=False,
                                              verbose=0)
    transform = np.loadtxt(register_result.pre_transforms[0])
    assert_array_almost_equal(transform,
                              [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0], decimal=1)
    assert_true(os.path.isfile(register_result.registered[0]))
    assert_true(register_result.transforms[0] is None)
