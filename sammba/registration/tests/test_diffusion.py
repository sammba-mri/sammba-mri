import os
from nose.tools import assert_true, assert_equal
from nose import with_setup
from nilearn.datasets.tests import test_utils as tst
from nilearn._utils.testing import assert_raises_regex
from nilearn._utils.niimg_conversions import _check_same_fov
from sammba.registration import DWISession
from sammba.registration.utils import _check_same_obliquity
from sammba.externals.nipype.interfaces import afni
from sammba import testing_data
import nibabel


def test_dwi_session_init():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'lemur_anat.nii.gz')
    dwi_file = os.path.join(os.path.dirname(testing_data.__file__),
                            'lemur_dwi.nii.gz')
    bvals_file = os.path.join(os.path.dirname(testing_data.__file__),
                              'bvals.txt')
    data = DWISession(anat=anat_file, dwi=dwi_file, bvals=bvals_file)
    assert_equal(data.anat, anat_file)
    assert_equal(data.dwi, dwi_file)
    assert_equal(data.bvals, bvals_file)


def test_dwi_session_check_inputs():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'lemur_anat.nii.gz')
    dwi_file = os.path.join(os.path.dirname(testing_data.__file__),
                            'lemur_dwi.nii.gz')
    bvals_file = os.path.join(os.path.dirname(testing_data.__file__),
                              'bvals.txt')
    data = DWISession(anat='', dwi=dwi_file, bvals=bvals_file,
                      animal_id='animal001')
    assert_raises_regex(IOError, "anat must be an existing image file",
                        data._check_inputs)
    data = DWISession(anat=anat_file, dwi='', bvals=bvals_file,
                      animal_id='animal001')
    assert_raises_regex(IOError, "dwi must be an existing image file",
                        data._check_inputs)
    data = DWISession(anat=anat_file, dwi=dwi_file, bvals='',
                      animal_id='animal001')
    assert_raises_regex(IOError, "bvals must be an existing b-values file",
                        data._check_inputs)
    data = DWISession(anat=anat_file, dwi=dwi_file, bvals=bvals_file)
    assert_raises_regex(ValueError, "animal_id must be a string",
                        data._check_inputs)


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_coregister():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'lemur_anat.nii.gz')
    dwi_file = os.path.join(os.path.dirname(testing_data.__file__),
                            'lemur_dwi.nii.gz')
    bvals_file = os.path.join(os.path.dirname(testing_data.__file__),
                              'bvals.txt')
    dwi_session = DWISession(anat=anat_file, func=dwi_file,
                             animal_id='test_coreg_dwi', brain_volume=1800,
                             tst.tmpdir)
    dwi_session.coregister(verbose=False, use_rats_tool=False)
    assert_true(_check_same_fov(nibabel.load(dwi_session.coreg_dwi_),
                                nibabel.load(dwi_session.coreg_anat_)))
    assert_true(_check_same_obliquity(dwi_session.coreg_anat_,
                                      dwi_session.coreg_dwi_))
    assert_true(os.path.isfile(dwi_session.coreg_transform_))
    assert_equal(os.path.join(tst.tmpdir, 'test_coreg_dwi'),
                 dwi_session.output_dir)

    # Check environement variables setting
    current_dir = os.getcwd()  # coregister_fmri_session changes the directory
    assert_raises_regex(RuntimeError,
                        "already exists",
                        dwi_session.coregister, tst.tmpdir,
                        use_rats_tool=False, AFNI_DECONFLICT='NO')
    os.chdir(current_dir)