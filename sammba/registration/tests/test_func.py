import os
import shutil
import tempfile
from nose.tools import assert_equal
from nilearn._utils.testing import assert_raises_regex
from sammba.registration import Subject, subjects_to_template
from sammba import testing_data


def test_subjects_to_template():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')
    subject = Subject(anat=anat_file, func=func_file)
    assert_equal(subject.anat, anat_file)
    assert_equal(subject.func, func_file)

    tempdir = tempfile.mkdtemp()
    write_dir = os.path.join(tempdir, 'test_func_dir')
    template_file = anat_file
    t_r = 1.
    assert_raises_regex(ValueError,
                        "'subjects' input argument must be an iterable",
                        subjects_to_template, subject, t_r,
                        template_file,
                        write_dir=write_dir)

    assert_raises_regex(ValueError,
                        "Each subject must have type",
                        subjects_to_template, [subject, 'subject'], t_r,
                        template_file,
                        write_dir=write_dir)

    assert_raises_regex(ValueError,
                        "Subjects ids must be different",
                        subjects_to_template, [subject, subject], t_r,
                        template_file,
                        write_dir=write_dir)

    if os.path.exists(write_dir):
        shutil.rmtree(write_dir)

    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)
