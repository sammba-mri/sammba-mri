import os
import shutil
import tempfile
from nose.tools import assert_true
from nilearn._utils.testing import assert_raises_regex
from nilearn._utils.niimg_conversions import _check_same_fov
from sammba.registration import FMRISession, func
from sammba.registration.utils import _have_same_obliquity
from sammba.externals.nipype.interfaces import afni
from sammba import testing_data
import nibabel


def test_coregister_fmri_session():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')
    tempdir = tempfile.mkdtemp()

    animal_session = FMRISession(anat=anat_file, func=func_file,
                                 animal_id='test_coreg_dir')
    if afni.Info().version():
        func.coregister_fmri_session(animal_session, 1., tempdir,
                                     slice_timing=False)
        assert_true(_check_same_fov(nibabel.load(animal_session.coreg_func_),
                                    nibabel.load(animal_session.coreg_anat_)))
        assert_true(_have_same_obliquity(animal_session.coreg_anat_,
                                         animal_session.coreg_func_))

    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)


def test_fmri_sessions_to_template():
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    func_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'func.nii.gz')
    mammal_data = FMRISession(anat=anat_file, func=func_file)

    tempdir = tempfile.mkdtemp()
    write_dir = os.path.join(tempdir, 'test_func_to_template_dir')
    template_file = anat_file
    t_r = 1.
    assert_raises_regex(ValueError,
                        "'animals_data' input argument must be an iterable",
                        func.fmri_sessions_to_template, mammal_data, t_r,
                        template_file,
                        write_dir=write_dir)

    assert_raises_regex(ValueError,
                        "Each animal data must have type",
                        func.fmri_sessions_to_template, [mammal_data, ''], t_r,
                        template_file,
                        write_dir=write_dir)

    assert_raises_regex(ValueError,
                        "Animals ids must be different",
                        func.fmri_sessions_to_template,
                        [mammal_data, mammal_data],
                        t_r, template_file, write_dir=write_dir)

    if afni.Info().version():
        registered_data = func.fmri_sessions_to_template([mammal_data], t_r,
                                                         template_file,
                                                         tempdir,
                                                         slice_timing=False)
        assert_true(os.path.isdir(registered_data.output_dir_))
        assert_true(os.path.isdir(registered_data.registered_func_))
        assert_true(os.path.isdir(registered_data.registered_anat_))

    if os.path.exists(write_dir):
        shutil.rmtree(write_dir)

    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)
