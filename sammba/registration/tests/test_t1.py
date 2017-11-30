import os
import shutil
import tempfile
import numpy as np
from numpy.testing import assert_array_almost_equal
from nilearn._utils.testing import assert_raises_regex
from sammba.externals.nipype.interfaces import afni
from sammba.registration import t1
from sammba import testing_data


def test_anats_to_common():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    tempdir = tempfile.mkdtemp()

    # Check error is raised if wrong registration kind
    assert_raises_regex(ValueError, "Registration kind must be one of ",
                        t1.anats_to_common, [anat_file], tempdir,
                        registration_kind='rigidd')

    # test common space of one image is itself
    if afni.Info().version():
        rigid = t1.anats_to_common([
            anat_file], tempdir, registration_kind='rigid', verbose=0)
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
        register_result = t1.anats_to_template(anat_file, anat_file, tempdir)
        transform = np.loadtxt(register_result.affine_transform)
        assert_array_almost_equal(transform,
                                  [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0])

    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)
