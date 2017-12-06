import os
import shutil
import tempfile
import numpy as np
from nose.tools import assert_true
from numpy.testing import assert_array_almost_equal
from nilearn._utils.testing import assert_raises_regex
import nibabel
from sammba.externals.nipype.interfaces import afni
from sammba.registration import struct
from nilearn._utils.niimg_conversions import _check_same_fov
from sammba.registration.utils import _check_same_obliquity
from sammba import testing_data


def test_anats_to_common():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    tempdir = tempfile.mkdtemp()

    # Check error is raised if wrong registration kind
    assert_raises_regex(ValueError, "Registration kind must be one of ",
                        struct.anats_to_common, [anat_file], tempdir, 400,
                        registration_kind='rigidd')

    # test common space of one image is itself
    if afni.Info().version():
        rigid = struct.anats_to_common([
            anat_file], tempdir, 400, registration_kind='rigid', verbose=0)
        transform = np.loadtxt(rigid.transforms[0])
        assert_array_almost_equal(transform,
                                  [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0])

    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)


def test_anat_to_template():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    tempdir = tempfile.mkdtemp()

    # test common space of one image is itself
    if afni.Info().version():
        register_result = struct.anats_to_template(anat_file, anat_file,
                                                   tempdir, 400)
        transform = np.loadtxt(register_result.affine_transform)
        assert_array_almost_equal(transform,
                                  [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0])
        assert_true(_check_same_fov(
            nibabel.load(register_result.registered_anat_),
            nibabel.load(anat_file)))
        assert_true(_check_same_obliquity(
            nibabel.load(register_result.registered_anat_,
            nibabel.load(anat_file))))
        assert_true(os.path.isfile(register_result.warp_transform_))

        os.remove(register_result.coreg_anat_)
        os.remove(register_result.affine_transform_)
        os.remove(register_result.warp_transform_)

    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)
