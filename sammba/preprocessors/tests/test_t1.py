import os
import shutil
import tempfile
import numpy as np
from numpy.testing import assert_array_almost_equal
from nilearn._utils.testing import assert_raises_regex
from sammba.preprocessors import t1
from sammba import testing_data


def test_register_to_common():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    tempdir = tempfile.mkdtemp()

    # Check error is raised if wrong registration kind
    assert_raises_regex(ValueError, "Registration kind must be one of ",
                        t1.register_to_common, [anat_file], tempdir,
                        registration_kind='rigidd')

    # test common space of one image is itself
    registration = t1.register_to_common([anat_file], tempdir,
                                         registration_kind='rigid', verbose=0)
    transform = np.loadtxt(registration.transforms[0])
    assert_array_almost_equal(transform, [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0])

    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)
