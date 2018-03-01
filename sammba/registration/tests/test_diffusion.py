import os
import glob
from nose.tools import assert_true, assert_equal, assert_not_equal
from nose import with_setup
import numpy as np
from numpy.testing import assert_array_equal
import nibabel
from nilearn.datasets.tests import test_utils as tst
from nilearn._utils.testing import assert_raises_regex
from nilearn._utils.niimg_conversions import _check_same_fov
from sammba.registration import DWISession
from sammba.registration.utils import _check_same_obliquity
from sammba import testing_data


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
    assert_equal(dwi_session.output_dir, tst.tmpdir)

    # Check output files
    output_dir_contents = glob.glob(os.path.join(dwi_session.output_dir, '*'))
    outputs = [dwi_session.coreg_transform_, dwi_session.coreg_anat_,
               dwi_session.coreg_dwi_, 'per_slice']
    per_slice_files = glob.glob(os.path.join(dwi_session.output_dir,
                                'per_slice', '*'))
    assert_array_equal(sorted(output_dir_contents).astype(str),
                       sorted(outputs))
    assert_equal(len(per_slice_files), 28)

    # Test same results with caching
    coregistered_dwi = dwi_session.coreg_dwi_
    coregistered_anat = dwi_session.coreg_anat_
    coregister_transform = dwi_session.coreg_transform_
    dwi_session.coregister(verbose=False, caching=True)

    assert_equal(dwi_session.coreg_dwi_, coregistered_dwi)
    assert_equal(dwi_session.coreg_anat_, coregistered_anat)
    assert_equal(dwi_session.coreg_transform_, coregister_transform)
    assert_equal(dwi_session.output_dir, tst.tmpdir)

    assert_array_equal(nibabel.load(dwi_session.coreg_dwi).get_data(),
                       nibabel.load(coregistered_dwi).get_data())
    assert_true(_check_same_obliquity(dwi_session.coreg_dwi_,
                                      coregistered_dwi))
    assert_array_equal(nibabel.load(dwi_session.coreg_anat).get_data(),
                       nibabel.load(coregistered_anat).get_data())
    assert_true(_check_same_obliquity(dwi_session.coreg_anat_,
                                      coregistered_anat))
    assert_array_equal(np.loadtxt(dwi_session.coreg_transform_),
                       np.loadtxt(coregister_transform))

    # Check environement variables setting
    assert_raises_regex(RuntimeError, "already exists", dwi_session.coregister,
                        verbose=False, AFNI_DECONFLICT='NO')
