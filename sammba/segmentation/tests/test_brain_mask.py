import os
from nose import with_setup
from nose.tools import assert_true, assert_equal
import nibabel
from nilearn.datasets.tests import test_utils as tst
from sammba.segmentation import brain_mask
from sammba import testing_data
from nilearn._utils.niimg_conversions import _check_same_fov


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_compute_histo_brain_mask():
    head_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    brain_mask_file = brain_mask.compute_brain_mask(
        head_file, 400, write_dir=tst.tmpdir, use_rats_tool=False)
    assert_true(os.path.isfile(brain_mask_file))
    assert_true(_check_same_fov(nibabel.load(brain_mask_file),
                                nibabel.load(head_file)))


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_brain_extraction_report():
    head_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    expected_report = u"""\
             x extent    y extent  z extent   volume

     0.2       0.83      0.79      0.81        24
     None      0.33      0.10      0.15        31
 """
    report = brain_mask.brain_extraction_report(head_file, 400,
                                                write_dir=tst.tmpdir,
                                                use_rats_tool=False)
    assert_equal(report, expected_report)
