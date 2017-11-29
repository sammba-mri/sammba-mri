import os
from nose.tools import assert_equal
from nilearn._utils.testing import assert_raises_regex
from sammba.registration import FMRIData
from sammba import testing_data


def test_subject_init():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')
    data = FMRIData(anat=anat_file, func=func_file)
    assert_equal(data.anat, anat_file)
    assert_equal(data.func, func_file)


def test_subject_check_inputs():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')

    data = FMRIData(anat='', func=func_file, animal_id='animal001')
    assert_raises_regex(IOError, "anat must be an existing image file",
                        data._check_inputs)

    data = FMRIData(anat=anat_file, func='', animal_id='animal001')
    assert_raises_regex(IOError, "func must be an existing image file",
                        data._check_inputs)

    data = FMRIData(anat=anat_file, func=func_file)
    assert_raises_regex(ValueError, "animal_id must be a string",
                        data._check_inputs)
