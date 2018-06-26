import os
from nose import with_setup
from nose.tools import assert_true
import nibabel
from nilearn.datasets.tests import test_utils as tst
from nilearn.image import index_img
from sammba.registration import base
from sammba import testing_data
from nilearn._utils.niimg_conversions import _check_same_fov


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_warp():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')
    func_file0 = os.path.join(tst.tmpdir, 'mean_func.nii.gz')
    func_img0 = index_img(func_file, 0)
    func_img0.to_filename(func_file0)

    registered_anat_oblique_file, mat_file =\
        base._warp(anat_file, func_file0, write_dir=tst.tmpdir,
                   caching=False, verbose=False)
    assert_true(_check_same_fov(nibabel.load(registered_anat_oblique_file),
                                func_img0))
    assert_true(os.path.isfile(mat_file))
