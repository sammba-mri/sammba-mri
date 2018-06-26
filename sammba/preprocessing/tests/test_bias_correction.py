import os
from nose import with_setup
from nose.tools import assert_true
import nibabel
from nilearn.datasets.tests import test_utils as tst
from sammba.preprocessing import bias_correction
from sammba import testing_data
from nilearn._utils.niimg_conversions import _check_same_fov


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_ants_n4():
    in_file = os.path.join(os.path.dirname(testing_data.__file__),
                           'func.nii.gz')
    unbiased_file = bias_correction.ants_n4(in_file, write_dir=tst.tmpdir,
                                            verbose=False)
    assert_true(_check_same_fov(nibabel.load(unbiased_file),
                                nibabel.load(in_file)))


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_afni_unifize():
    in_file = os.path.join(os.path.dirname(testing_data.__file__),
                           'func.nii.gz')
    unbiased_file = bias_correction.afni_unifize(in_file,
                                                 write_dir=tst.tmpdir,
                                                 verbose=False)
    assert_true(_check_same_fov(nibabel.load(unbiased_file),
                                nibabel.load(in_file)))
