import os
import json
from nose import with_setup
from nose.tools import assert_equal, assert_not_equal
from nilearn.datasets import utils
from nilearn.datasets.tests import test_utils as tst

from sammba.data_fetchers import struct


def setup_mock():
    return tst.setup_mock(utils, struct)


def teardown_mock():
    return tst.teardown_mock(utils, struct)


@with_setup(setup_mock, teardown_mock)
@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_fetch_atlas_dorr_2008():
    datadir = os.path.join(tst.tmpdir, 'dorr_2008')
    os.mkdir(datadir)
    dummy = open(os.path.join(
        datadir, 'c57_brain_atlas_labels.csv'), 'w')
    dummy.write("\n1,amygdala,51,151\n27,fourth ventricle,118,118")
    dummy.close()
    bunch = struct.fetch_atlas_dorr_2008(data_dir=tst.tmpdir, verbose=0)

    assert_equal(len(tst.mock_url_request.urls), 2)
    assert_equal(bunch['t2'],
                 os.path.join(datadir, 'male-female-mouse-atlas.mnc'))
    assert_equal(bunch['maps'],
                 os.path.join(datadir, 'c57_fixed_labels_resized.mnc'))
    assert_not_equal(bunch.description, '')


@with_setup(setup_mock, teardown_mock)
@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_fetch_atlas_waxholm_rat_2014():
    datadir = os.path.join(tst.tmpdir, 'waxholm_rat_2014')
    os.mkdir(datadir)
    # Create dummy json labels file
    json_filename = os.path.join(datadir, 'WHS_SD_rat_labels.json')
    with open(json_filename, 'w') as json_content:
        json.dump({"216": "Spinal cord"}, json_content)

    for downsample in ['2', '3']:
        bunch = struct.fetch_atlas_waxholm_rat_2014(
            data_dir=tst.tmpdir, verbose=0, downsample=downsample)
        downsampled_atlas = 'WHS_SD_rat_atlas_v1.01_' +\
                            'downsample{0}.nii.gz'.format(downsample)
        downsampled_t2star = 'WHS_SD_rat_T2star_v1.01_' +\
                             'downsample{0}.nii.gz'.format(downsample)
        assert_equal(bunch['t2star'],
                     os.path.join(datadir, downsampled_t2star))
        assert_equal(bunch['maps'],
                     os.path.join(datadir, downsampled_atlas))

    assert_equal(len(tst.mock_url_request.urls), 4)
    assert_not_equal(bunch.description, '')
