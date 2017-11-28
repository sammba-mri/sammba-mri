import os
from nose.tools import assert_equal
from nilearn._utils.testing import assert_raises_regex
from sammba.registration import Subject
from sammba import testing_data


def test_subject_init():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')
    subject = Subject(anat=anat_file, func=func_file)
    assert_equal(subject.anat, anat_file)
    assert_equal(subject.func, func_file)


def test_subject_check_inputs():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')

    subject = Subject(anat='', func=func_file)
    assert_raises_regex(IOError, "anat must be an existing image file",
                        subject._check_inputs)
    subject = Subject(anat=anat_file, func='')
    assert_raises_regex(IOError, "func must be an existing image file",
                        subject._check_inputs)
