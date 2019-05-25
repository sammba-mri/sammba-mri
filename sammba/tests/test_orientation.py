import os
import shutil
from numpy.testing import assert_array_equal, assert_array_almost_equal
from nose.tools import assert_true, assert_false
from nose import with_setup
import nibabel
from nilearn.datasets.tests import test_utils as tst
from nilearn._utils.testing import assert_raises_regex
from nipype.interfaces import afni
from sammba import orientation
from sammba import testing_data


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_reset_affines():
    # test qform is replaced by sfrom with default parameters
    in_file = os.path.join(os.path.dirname(testing_data.__file__),
                           'anat.nii.gz')
    source_file = os.path.join(os.path.dirname(testing_data.__file__),
                               'anat.nii.gz')
    in_file = os.path.join(tst.tmpdir, 'in_anat.nii.gz')
    shutil.copy(source_file, in_file)
    in_header = nibabel.load(in_file).header
    in_sform, in_scode = in_header.get_sform(coded=True)
    out_file = os.path.join(tst.tmpdir, 'out_anat.nii.gz')
    orientation._reset_affines(in_file, out_file)
    out_header = nibabel.load(out_file).header
    out_sform, out_scode = out_header.get_sform(coded=True)
    out_qform, out_qcode = out_header.get_qform(coded=True)

    assert_array_equal(out_sform, in_sform)
    assert_array_almost_equal(out_qform, out_sform)
    assert_array_equal(out_qcode, out_scode)

    # same test, overwriting
    orientation._reset_affines(in_file, out_file, overwrite=True)
    out_header = nibabel.load(out_file).header
    assert_array_equal(out_sform, in_sform)
    assert_array_almost_equal(out_qform, out_sform)
    assert_array_equal(out_qcode, out_scode)

    # previous test, with verbose=0
    orientation._reset_affines(in_file, out_file, overwrite=True, verbose=0)
    out_header = nibabel.load(out_file).header
    assert_array_equal(out_sform, in_sform)
    assert_array_almost_equal(out_qform, out_sform)
    assert_array_equal(out_qcode, out_scode)


def test_check_same_obliquity():
    img_filename1 = os.path.join(os.path.dirname(testing_data.__file__),
                                 'anat.nii.gz')
    img_filename2 = os.path.join(os.path.dirname(testing_data.__file__),
                                 'func.nii.gz')
    assert_true(orientation._check_same_obliquity(img_filename1, img_filename1))
    assert_false(orientation._check_same_obliquity(img_filename1, img_filename2))


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_fix_obliquity():
    target_filename = os.path.join(os.path.dirname(testing_data.__file__),
                                   'anat.nii.gz')
    # ZCutUp removes the obliquity
    tmp_filename = os.path.join(tst.tmpdir, 'img_test_obliquity.nii.gz')
    slicer = afni.ZCutUp().run
    _ = slicer(in_file=target_filename,
               keep='0 27',
               out_file=tmp_filename)
    assert_false(
        orientation._check_same_obliquity(tmp_filename, target_filename))
    tmp_filename_oblique = orientation.fix_obliquity(tmp_filename,
                                                     target_filename)
    assert_true(
        orientation._check_same_obliquity(tmp_filename_oblique,
                                          target_filename))


def test_check_same_geometry():
    img_filename1 = os.path.join(os.path.dirname(testing_data.__file__),
                                 'anat.nii.gz')
    img_filename2 = os.path.join(os.path.dirname(testing_data.__file__),
                                 'func.nii.gz')
    assert_true(orientation._check_same_geometry(img_filename1, img_filename1))
    assert_false(orientation._check_same_geometry(img_filename1, img_filename2))


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_copy_geometry():
    func_filename = os.path.join(os.path.dirname(testing_data.__file__),
                                 'func.nii.gz')
    anat_filename = os.path.join(os.path.dirname(testing_data.__file__),
                                 'anat.nii.gz')

    # Check error is raised if resampling is not allowed
    assert_raises_regex(
        ValueError, 'images have different shapes', orientation.copy_geometry,
        anat_filename, func_filename,
        out_filename=os.path.join(tst.tmpdir, 'geometry_test.nii.gz'),
        in_place=False)

    changed_filename = orientation.copy_geometry(
        anat_filename, func_filename,
        out_filename=os.path.join(tst.tmpdir, 'geometry_test.nii.gz'),
        in_place=False, allow_resampling=True)
    assert_false(orientation._check_same_geometry(func_filename, anat_filename))
    assert_true(orientation._check_same_geometry(changed_filename, anat_filename))