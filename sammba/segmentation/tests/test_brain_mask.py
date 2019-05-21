import os
from nose import with_setup
from nose.tools import assert_true, assert_equal, assert_greater
import nibabel
import numpy as np
from nilearn import masking
from nilearn.image import coord_transform
from nilearn.datasets.tests import test_utils as tst
from nilearn._utils.niimg_conversions import _check_same_fov
from sammba import testing_data
from sammba.segmentation import brain_mask


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_compute_histo_brain_mask():
    head_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    brain_mask_file = brain_mask.compute_histo_brain_mask(
        head_file, 400, write_dir=tst.tmpdir,
        verbose=False)
    assert_true(os.path.isfile(brain_mask_file))
    assert_true(_check_same_fov(nibabel.load(brain_mask_file),
                                nibabel.load(head_file)))


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_get_mask_measures():
    # Create ellipsoid image
    one_mask_data = np.ones((60, 60, 80))
    affine = np.eye(4)
    affine[0, 0] = .7
    affine[2, 2] = .5
    affine[:3, 3] = np.array([-2, 7, 3])
    one_mask_img = nibabel.Nifti1Image(one_mask_data, affine)
    i, j, k = np.where(one_mask_data <2)   
    voxels_coords = np.array(coord_transform(i, j, k, affine)).T
    center_coords = np.array([ 20.,  31.,  20.])
    positions = voxels_coords - center_coords
    angle = np.pi / 3
    c = np.cos(angle)
    s = np.sin(angle)
    rotation = np.array([[c, -s, 0], [s, c, 0.], [0, 0, 1]])
    lambd = np.eye(3)
    a = 6.
    b = 15.
    c = 16.
    lambd[0, 0] = 1. / a ** 2
    lambd[2, 2] = 1. / b ** 2
    lambd[1, 1] = 1. / c ** 2
    trans_sqrt = np.sqrt(lambd).dot(rotation)
    transformed_positions = trans_sqrt.dot(positions.T).T
    new_voxels = np.sum(transformed_positions ** 2, axis=1) <= 1
    mask_file = os.path.join(tst.tmpdir, 'ellipsoid.nii.gz')
    masking.unmask(new_voxels, one_mask_img).to_filename(mask_file)
    length_ap, length_rl, length_is, symmetry, volume = brain_mask._get_mask_measures(mask_file)
    np.testing.assert_array_almost_equal(length_ap, 2 * c)
    np.testing.assert_array_almost_equal(length_rl, 2 * b, decimal=0)
    np.testing.assert_array_almost_equal(length_is, 2 * a, decimal=0)
    assert_greater(symmetry, .99)
    #np.testing.assert_array_almost_equal(volume, 4. * np.pi /3. * a * b * c)


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_brain_extraction_report():
    head_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    expected_report = u"""\
                 AP length   RL length   IS length    symmetry      volume

fraction 0.20        13.50        9.90        4.50        0.91      381.32
  no fraction        17.40       12.30        5.40        1.00     1141.34

"""
    report = brain_mask.brain_extraction_report(head_file, 400,
                                                write_dir=tst.tmpdir,
                                                use_rats_tool=False)
    assert_equal(report, expected_report)
