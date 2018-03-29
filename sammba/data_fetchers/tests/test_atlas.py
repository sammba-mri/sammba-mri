import os
import json
import numpy as np
from numpy.testing import assert_array_almost_equal
from nose import with_setup
from nose.tools import assert_equal, assert_not_equal, assert_true
import nibabel
from nilearn.datasets import utils
from nilearn.datasets.tests import test_utils as tst
from nilearn._utils.niimg_conversions import check_niimg
from sammba.data_fetchers import atlas
from sammba import testing_data


def setup_mock():
    return tst.setup_mock(utils, atlas)


def teardown_mock():
    return tst.teardown_mock(utils, atlas)


@with_setup(setup_mock, teardown_mock)
@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_fetch_atlas_dorr_2008():
    datadir = os.path.join(tst.tmpdir, 'dorr_2008')
    os.mkdir(datadir)
    dummy = open(os.path.join(
        datadir, 'c57_brain_atlas_labels.csv'), 'w')
    dummy.write("\n1,amygdala,51,151\n27,fourth ventricle,118,118")
    dummy.close()

    # Default resolution
    bunch = atlas.fetch_atlas_dorr_2008(data_dir=tst.tmpdir, verbose=0)
    assert_equal(len(tst.mock_url_request.urls), 2)
    assert_equal(bunch['t2'],
                 os.path.join(datadir, 'Dorr_2008_average.nii.gz'))
    assert_equal(bunch['maps'],
                 os.path.join(datadir, 'Dorr_2008_labels.nii.gz'))

    # test resampling
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    anat_img = check_niimg(anat_file)
    anat_img.to_filename(bunch['t2'])
    anat_img = check_niimg(anat_file, dtype=int)
    anat_img.to_filename(bunch['maps'])

    bunch = atlas.fetch_atlas_dorr_2008(
        data_dir=tst.tmpdir, verbose=0, downsample='100')
    assert_equal(bunch['t2'],
                 os.path.join(datadir, 'Dorr_2008_average_100um.nii.gz'))
    assert_equal(bunch['maps'],
                 os.path.join(datadir, 'Dorr_2008_labels_100um.nii.gz'))
    assert_array_almost_equal(nibabel.load(bunch['t2']).header.get_zooms(),
                              (.1, .1, .1))
    assert_array_almost_equal(nibabel.load(bunch['maps']).header.get_zooms(),
                              (.1, .1, .1))
    assert_equal(nibabel.load(bunch['maps']).get_data().dtype, np.dtype(int))
    assert_equal(len(tst.mock_url_request.urls), 2)
    assert_equal(len(bunch['names']), 3)
    assert_equal(len(bunch['labels']), 3)

    # test with 'minc' format
    bunch = atlas.fetch_atlas_dorr_2008(data_dir=tst.tmpdir, verbose=0,
                                        image_format='minc')
    assert_equal(len(tst.mock_url_request.urls), 4)
    assert_equal(bunch['t2'],
                 os.path.join(datadir, 'male-female-mouse-atlas.mnc'))
    assert_equal(bunch['maps'],
                 os.path.join(datadir, 'c57_fixed_labels_resized.mnc'))
    assert_equal(len(bunch['names']), 3)
    assert_equal(len(bunch['labels']), 3)
    assert_not_equal(bunch.description, '')


@with_setup(setup_mock, teardown_mock)
@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_fetch_masks_dorr_2008():
    # create a dummy csv file for dorr labels
    datadir = os.path.join(tst.tmpdir, 'dorr_2008')
    os.mkdir(datadir)
    dummy_csv = open(os.path.join(
        datadir, 'c57_brain_atlas_labels.csv'), 'w')
    dummy_csv.write("Structure,right label,left label"
                    "\n1,amygdala,51,151\n19,corpus callosum,8,68"
                    "\n27,fourth ventricle,118,118"
                    "\n38,lateral ventricle,57,77")
    dummy_csv.close()
    dorr_atlas = atlas.fetch_atlas_dorr_2008(data_dir=tst.tmpdir, verbose=0)

    # create a dummy atlas image
    dummy_atlas_data = np.zeros((100, 100, 100))
    dummy_atlas_data[:10, :10, :10] = 51
    dummy_atlas_data[50:90, 50:90, 50:90] = 151
    dummy_atlas_data[10:20, :10, 10:30] = 8
    dummy_atlas_data[90:100, :10, 10:30] = 68
    dummy_atlas_data[10:20, 50:90, 10:20] = 118
    dummy_atlas_data[40:60, 30:40, 40:50] = 57
    dummy_atlas_data[60:70, 30:40, 40:50] = 77
    dummy_atlas_img = nibabel.Nifti1Image(dummy_atlas_data, np.eye(4))
    dummy_atlas_img.to_filename(dorr_atlas.maps)

    dorr_masks = atlas.fetch_masks_dorr_2008(data_dir=tst.tmpdir, verbose=0)
    assert_true(isinstance(dorr_masks.brain, nibabel.Nifti1Image))
    assert_true(isinstance(dorr_masks.gm, nibabel.Nifti1Image))
    assert_true(isinstance(dorr_masks.cc, nibabel.Nifti1Image))
    assert_true(isinstance(dorr_masks.ventricles, nibabel.Nifti1Image))


@with_setup(setup_mock, teardown_mock)
@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_fetch_atlas_waxholm_rat_2014():
    datadir = os.path.join(tst.tmpdir, 'waxholm_rat_2014')
    os.mkdir(datadir)
    # Create dummy json labels file
    json_filename = os.path.join(datadir, 'WHS_SD_rat_labels.json')
    with open(json_filename, 'w') as json_content:
        json.dump({"216": "Spinal cord"}, json_content)

    # default resolution
    bunch = atlas.fetch_atlas_waxholm_rat_2014(
        data_dir=tst.tmpdir, verbose=0)
    assert_equal(
        bunch['t2star'],
        os.path.join(datadir, 'WHS_SD_rat_T2star_v1_01_downsample3.nii.gz'))
    assert_equal(
        bunch['maps'],
        os.path.join(datadir, 'WHS_SD_rat_atlas_v1_01_downsample3.nii.gz'))
    assert_equal(len(tst.mock_url_request.urls), 2)

    assert_not_equal(bunch.description, '')

    # Downsampled 2 times
    bunch = atlas.fetch_atlas_waxholm_rat_2014(
        data_dir=tst.tmpdir, verbose=0, downsample='78')
    assert_equal(
        bunch['t2star'],
        os.path.join(datadir, 'WHS_SD_rat_T2star_v1_01_downsample2.nii.gz'))
    assert_equal(
        bunch['maps'],
        os.path.join(datadir, 'WHS_SD_rat_atlas_v1_01_downsample2.nii.gz'))
    assert_equal(len(tst.mock_url_request.urls), 4)

    # test resampling
    anat_file = os.path.join(os.path.dirname(testing_data.__file__),
                             'anat.nii.gz')
    anat_img = check_niimg(anat_file)
    anat_img.to_filename(bunch['t2star'])
    anat_img = check_niimg(anat_file, dtype=int)
    anat_img.to_filename(bunch['maps'])
    bunch = atlas.fetch_atlas_waxholm_rat_2014(
        data_dir=tst.tmpdir, verbose=0, downsample='200')
    assert_equal(len(tst.mock_url_request.urls), 4)
    assert_equal(
        bunch['t2star'],
        os.path.join(datadir, 'WHS_SD_rat_T2star_v1_01_200um.nii.gz'))
    assert_equal(
        bunch['maps'],
        os.path.join(datadir, 'WHS_SD_rat_atlas_v1_01_200um.nii.gz'))
    assert_array_almost_equal(nibabel.load(bunch['t2star']).header.get_zooms(),
                              (.2, .2, .2))
    assert_array_almost_equal(nibabel.load(bunch['maps']).header.get_zooms(),
                              (.2, .2, .2))


@with_setup(setup_mock, teardown_mock)
@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_fetch_atlas_lemur_mircen_2017():
    datadir = os.path.join(tst.tmpdir, 'mircen_2017')
    os.mkdir(datadir)
    dummy = open(os.path.join(
        datadir, 'MIRCen_mouselemur_atlas_labels.txt'), 'w')
    dummy.write('#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n'
                '    0     0    0    0        0  0  0    "Clear Label"\n'
                '    1   248  164  138        1  1  0    "hippocampus R"\n'
                '    2    96  238  253        1  1  0    "hippocampus L"\n')
    dummy.close()

    # Default resolution
    bunch = atlas.fetch_atlas_lemur_mircen_2017(data_dir=tst.tmpdir, verbose=0)
    assert_equal(len(tst.mock_url_request.urls), 2)
    assert_equal(bunch['t2'],
                 os.path.join(datadir, 'MIRCen_mouselemur_template.nii.gz'))
    assert_equal(bunch['maps'],
                 os.path.join(datadir, 'MIRCen_mouselemur_atlas.nii.gz'))
    assert_equal(nibabel.load(bunch['maps']).get_data().dtype, np.dtype(int))
    assert_equal(len(tst.mock_url_request.urls), 1)
    assert_equal(len(bunch['labels']), 3)
    assert_not_equal(bunch.description, '')
