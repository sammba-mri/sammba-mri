import os
import tempfile
from nose.tools import assert_true, assert_equal
from nilearn._utils.testing import assert_raises_regex
from nilearn._utils.niimg_conversions import _check_same_fov
from sammba.registration import FMRISession, func
from sammba.registration.utils import _check_same_obliquity
from sammba import testing_data
import nibabel


def test_coregister_fmri_session():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')

    animal_session = FMRISession(anat=anat_file, func=func_file,
                                 animal_id='test_coreg_dir')

    tempdir = tempfile.mkdtemp()
    func.coregister_fmri_session(animal_session, 1., tempdir, 400,
                                 slice_timing=False, verbose=False,
                                 use_rats_tool=False)
    assert_true(_check_same_fov(nibabel.load(animal_session.coreg_func_),
                                nibabel.load(animal_session.coreg_anat_)))
    assert_true(_check_same_obliquity(animal_session.coreg_anat_,
                                      animal_session.coreg_func_))
    assert_true(os.path.isfile(animal_session.coreg_transform_))
    assert_equal(os.path.join(tempdir, 'test_coreg_dir'),
                 animal_session.output_dir_)

    if os.path.exists(tempdir):
        for out_file in os.listdir(animal_session.output_dir_):
            os.remove(os.path.join(animal_session.output_dir_, out_file))
        os.removedirs(animal_session.output_dir_)

    # Check environement variables setting
    tempdir = tempfile.mkdtemp(suffix='?')
    assert_raises_regex(RuntimeError,
                        "badly formed filename",
                        func.coregister_fmri_session, animal_session, 1.,
                        tempdir, 400, slice_timing=False, use_rats_tool=False)
    assert_raises_regex(RuntimeError,
                        "Program Death",
                        func.coregister_fmri_session, animal_session, 1.,
                        tempdir, 400, slice_timing=False, use_rats_tool=False,
                        AFNI_ALLOW_ARBITRARY_FILENAMES='YES')

    if os.path.exists(tempdir):
        output_dir = os.path.join(tempdir, 'animal001')
        for out_file in os.listdir(output_dir):
            os.remove(os.path.join(output_dir, out_file))
        os.removedirs(output_dir)


def test_fmri_sessions_to_template():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')
    mammal_data = FMRISession(anat=anat_file, func=func_file)

    tempdir = tempfile.mkdtemp()
    template_file = anat_file
    t_r = 1.
    brain_volume = 400
    assert_raises_regex(ValueError,
                        "'animals_data' input argument must be an iterable",
                        func.fmri_sessions_to_template, mammal_data, t_r,
                        template_file, tempdir, brain_volume)

    assert_raises_regex(ValueError,
                        "Each animal data must have type",
                        func.fmri_sessions_to_template, [mammal_data, ''], t_r,
                        template_file, tempdir, brain_volume)

    assert_raises_regex(ValueError,
                        "Animals ids must be different",
                        func.fmri_sessions_to_template,
                        [mammal_data, mammal_data],
                        t_r, template_file, tempdir, brain_volume)

    registered_data = func.fmri_sessions_to_template([mammal_data], t_r,
                                                     template_file,
                                                     tempdir,
                                                     brain_volume,
                                                     slice_timing=False,
                                                     verbose=False,
                                                     use_rats_tool=False)
    assert_true(os.path.isdir(registered_data[0].output_dir_))
    assert_true(os.path.isfile(registered_data[0].registered_func_))
    assert_true(os.path.isfile(registered_data[0].registered_anat_))

    if os.path.exists(tempdir):
        for out_file in os.listdir(registered_data[0].output_dir_):
            os.remove(os.path.join(registered_data[0].output_dir_, out_file))
        os.removedirs(registered_data[0].output_dir_)
