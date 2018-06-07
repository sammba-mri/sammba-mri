import os
from nose import with_setup
from nose.tools import assert_true
from nilearn.datasets.tests import test_utils as tst
from sammba.registration import base
from sammba import testing_data


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_compute_brain_mask():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    brain_mask_file = base.compute_brain_mask(
        anat_file, 400, tst.tmpdir, use_rats_tool=False)
    assert_true(os.path.isfile(brain_mask_file))
