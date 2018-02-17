import os
from nose.tools import assert_true, assert_equal
from nose import with_setup
from nilearn.datasets.tests import test_utils as tst
from nilearn._utils.testing import assert_raises_regex
from nilearn._utils.niimg_conversions import _check_same_fov
from sammba.registration import DWISession
from sammba.registration.utils import _check_same_obliquity
from sammba import testing_data
import nibabel


def test_dwi_session_init():
    data_dir = os.path.dirname(testing_data.__file__)
    anat_file = os.path.join(data_dir, 'lemur_anat.nii.gz')
    dwi_file = os.path.join(data_dir, 'lemur_dwi.nii.gz')
    bvals_file = os.path.join(data_dir, 'bvals.txt')

    data = DWISession(anat=anat_file, dwi=dwi_file, bvals=bvals_file)
    assert_equal(data.anat, anat_file)
    assert_equal(data.dwi, dwi_file)
    assert_equal(data.bvals, bvals_file)


def test_dwi_session_check_inputs():
    data_dir = os.path.dirname(testing_data.__file__)
    anat_file = os.path.join(data_dir, 'lemur_anat.nii.gz')
    dwi_file = os.path.join(data_dir, 'lemur_dwi.nii.gz')
    bvals_file = os.path.join(data_dir, 'bvals.txt')

    dwi_session = DWISession(anat='', dwi=dwi_file, bvals=bvals_file)
    assert_raises_regex(IOError, "anat must be an existing image file",
                        dwi_session._check_inputs)

    dwi_session = DWISession(anat=anat_file, dwi='', bvals=bvals_file)
    assert_raises_regex(IOError, "dwi must be an existing image file",
                        dwi_session._check_inputs)

    dwi_session = DWISession(anat=anat_file, dwi=dwi_file, bvals='',)
    assert_raises_regex(IOError, "bvals must be an existing b-values file",
                        dwi_session._check_inputs)


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_coregister():
    data_dir = os.path.dirname(testing_data.__file__)
    anat_file = os.path.join(data_dir, 'lemur_anat.nii.gz')
    dwi_file = os.path.join(data_dir, 'lemur_dwi.nii.gz')
    bvals_file = os.path.join(data_dir, 'bvals.txt')
    dwi_session = DWISession(anat=anat_file, dwi=dwi_file, bvals=bvals_file,
                             output_dir=tst.tmpdir)

    assert_raises_regex(ValueError,
                        "'brain_volume' value is needed to perform rigid",
                        dwi_session.coregister,
                        prior_rigid_body_registration=True)

    dwi_session.coregister(verbose=False)
    assert_true(_check_same_fov(nibabel.load(dwi_session.coreg_dwi_),
                                nibabel.load(dwi_session.coreg_anat_)))
    assert_true(_check_same_obliquity(dwi_session.coreg_anat_,
                                      dwi_session.coreg_dwi_))
    assert_true(os.path.isfile(dwi_session.coreg_transform_))
    assert_equal(tst.tmpdir, dwi_session.output_dir)

    # Check environement variables setting
    assert_raises_regex(RuntimeError, "already exists", dwi_session.coregister,
                        verbose=False, AFNI_DECONFLICT='NO')
