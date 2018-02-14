import os
from nose.tools import assert_true, assert_equal
from nose import with_setup
from nilearn.datasets.tests import test_utils as tst
from nilearn._utils.testing import assert_raises_regex
from nilearn._utils.niimg_conversions import _check_same_fov
from sammba.registration import FMRISession, func
from sammba.registration.utils import _check_same_obliquity
from sammba import testing_data
import nibabel


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_coregister_fmri_session():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')

    animal_session = FMRISession(anat=anat_file, func=func_file,
                                 animal_id='test_coreg_dir')

    func.coregister_fmri_session(animal_session, 1., tst.tmpdir, 600,
                                 slice_timing=False, verbose=False,
                                 use_rats_tool=False)
    assert_true(_check_same_fov(nibabel.load(animal_session.coreg_func_),
                                nibabel.load(animal_session.coreg_anat_)))
    assert_true(_check_same_obliquity(animal_session.coreg_anat_,
                                      animal_session.coreg_func_))
    assert_true(os.path.isfile(animal_session.coreg_transform_))
    assert_equal(os.path.join(tst.tmpdir, 'test_coreg_dir'),
                 animal_session.output_dir_)

    # Check environement variables setting
    current_dir = os.getcwd()  # coregister_fmri_session changes the directory
    assert_raises_regex(RuntimeError,
                        "already exists",
                        func.coregister_fmri_session, animal_session, 1.,
                        tst.tmpdir, 400, slice_timing=False,
                        use_rats_tool=False, AFNI_DECONFLICT='NO')
    os.chdir(current_dir)



@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_fmri_sessions_to_template():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')
    mammal_data = FMRISession(anat=anat_file, func=func_file)

    template_file = anat_file
    t_r = 1.
    brain_volume = 400

    registered_data = func.fmri_sessions_to_template([mammal_data], t_r,
                                                     template_file,
                                                     tst.tmpdir,
                                                     brain_volume,
                                                     slice_timing=False,
                                                     verbose=True,
                                                     use_rats_tool=False)
    assert_true(os.path.isdir(registered_data[0].output_dir_))
    assert_true(os.path.isfile(registered_data[0].registered_func_))
    assert_true(os.path.isfile(registered_data[0].registered_anat_))

    assert_raises_regex(ValueError,
                        "'animals_data' input argument must be an iterable",
                        func.fmri_sessions_to_template, mammal_data, t_r,
                        template_file, tst.tmpdir, brain_volume)

    assert_raises_regex(ValueError,
                        "Each animal data must have type",
                        func.fmri_sessions_to_template, [mammal_data, ''], t_r,
                        template_file, tst.tmpdir, brain_volume)

    assert_raises_regex(ValueError,
                        "Animals ids must be different",
                        func.fmri_sessions_to_template,
                        [mammal_data, mammal_data],
                        t_r, template_file, tst.tmpdir, brain_volume)
