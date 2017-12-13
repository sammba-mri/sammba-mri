import os
import tempfile
from nose.tools import assert_true, assert_equal
import nibabel
from nilearn._utils.testing import assert_raises_regex
from nilearn._utils.niimg_conversions import _check_same_fov
from sammba.registration import fmri_session
from sammba import testing_data
from sammba.registration.utils import _check_same_obliquity
from sammba.externals.nipype.interfaces import afni


def test_fmri_session_init():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')
    data = fmri_session.FMRISession(anat=anat_file, func=func_file)
    assert_equal(data.anat, anat_file)
    assert_equal(data.func, func_file)


def test_fmri_session_check_inputs():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')

    data = fmri_session.FMRISession(anat='', func=func_file, animal_id='animal001')
    assert_raises_regex(IOError, "anat must be an existing image file",
                        data._check_inputs)

    data = fmri_session.FMRISession(anat=anat_file, func='', animal_id='animal001')
    assert_raises_regex(IOError, "func must be an existing image file",
                        data._check_inputs)

    data = fmri_session.FMRISession(anat=anat_file, func=func_file)
    assert_raises_regex(ValueError, "animal_id must be a string",
                        data._check_inputs)


def test_coregister():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')

    animal_session = fmri_session.FMRISession(anat=anat_file, func=func_file,
                                 animal_id='test_coreg_dir')

    if afni.Info().version():
        tempdir = tempfile.mkdtemp()
        func.coregister_fmri_session(animal_session, 1., tempdir, 400,
                                     slice_timing=False, verbose=False)
        assert_true(_check_same_fov(nibabel.load(animal_session.coreg_func_),
                                    nibabel.load(animal_session.coreg_anat_)))
        assert_true(_check_same_obliquity(animal_session.coreg_anat_,
                                          animal_session.coreg_func_))
        assert_true(os.path.isfile(animal_session.coreg_transform_))
        assert_equal(tempdir, animal_session.output_dir_)

        os.remove(animal_session.coreg_transform_)
        os.remove(animal_session.coreg_anat_)
        os.remove(animal_session.coreg_func_)
        os.removedirs(tempdir)

        # Check environement variables setting
        tempdir = tempfile.mkdtemp(suffix='$')
        assert_raises_regex(RuntimeError,
                            "Illegal dataset name",
                            func.coregister_fmri_session, animal_session, 1.,
                            tempdir, 400, slice_timing=False)
        func.coregister_fmri_session(
            animal_session, 1., tempdir, 400, slice_timing=False,
            AFNI_ALLOW_ARBITRARY_FILENAMES='YES',
            verbose=False)
        assert_true(_check_same_fov(nibabel.load(animal_session.coreg_func_),
                                    nibabel.load(animal_session.coreg_anat_)))
        assert_true(_check_same_obliquity(animal_session.coreg_anat_,
                                          animal_session.coreg_func_))
        assert_true(os.path.isfile(animal_session.coreg_transform_))
        assert_equal(tempdir, animal_session.output_dir_)

        os.remove(animal_session.coreg_transform_)
        os.remove(animal_session.coreg_anat_)
        os.remove(animal_session.coreg_func_)
        os.removedirs(tempdir)


def test_register_to_template():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')
    mammal_data = fmri_session.FMRISession(anat=anat_file, func=func_file)

    tempdir = tempfile.mkdtemp()
    template_file = anat_file
    t_r = 1.
    brain_volume = 400
    assert_raises_regex(ValueError,
                        "'animals_data' input argument must be an iterable",
                        fmri_session.register_to_template, t_r, brain_volume,
                        template_file, tempdir)

    assert_raises_regex(ValueError,
                        "Each animal data must have type",
                        func.fmri_sessions_to_template, [mammal_data, ''], t_r,
                        template_file, tempdir, brain_volume)

    assert_raises_regex(ValueError,
                        "Animals ids must be different",
                        func.fmri_sessions_to_template,
                        [mammal_data, mammal_data],
                        t_r, template_file, tempdir, brain_volume)

    if afni.Info().version():
        registered_data = func.fmri_sessions_to_template([mammal_data], t_r,
                                                         template_file,
                                                         tempdir,
                                                         brain_volume,
                                                         slice_timing=False,
                                                         verbose=False)
        assert_true(os.path.isdir(registered_data.output_dir_))
        assert_true(os.path.isfile(registered_data.registered_func_))
        assert_true(os.path.isfile(registered_data.registered_anat_))

    if os.path.exists(tempdir):
        os.removedirs(tempdir)
