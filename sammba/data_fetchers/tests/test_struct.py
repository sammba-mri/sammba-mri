import os
import numpy as np
from nose import with_setup
from nose.tools import assert_equal, assert_true, assert_not_equal
from nilearn.datasets import utils
from nilearn.datasets.tests import test_utils as tst
from sammba.data_fetchers import struct


def setup_mock():
    return tst.setup_mock(utils, struct)


def teardown_mock():
    return tst.teardown_mock(utils, struct)


@with_setup(setup_mock, teardown_mock)
@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_fetch_lemur_mircen_2019_t2():
    local_url = "file://" + tst.datadir
    # Create a dummy urls file
    lemur_t2_dir = os.path.join(tst.tmpdir, 'mircen2019_t2')
    os.mkdir(lemur_t2_dir)
    ids = ['"sub-{:02d}"'.format(idx) for idx in range(1, 35)]
    pheno = np.vstack((ids, ['m'] * 34)).T
    np.savetxt(os.path.join(lemur_t2_dir, 'lemur_atlas_list_t2_bids.csv'),
               pheno, fmt='%s\t%s', header='animal_id\tgender')

    # By default one subject is loaded
    lemur_t2 = struct.fetch_lemur_mircen_2019_t2(data_dir=tst.tmpdir,
                                              url=local_url, verbose=0)
    assert_equal(len(tst.mock_url_request.urls), 2)
    assert_equal(len(lemur_t2.anat), 1)

    # Load 3 subjects
    tst.mock_url_request.reset()
    lemur_t2 = struct.fetch_lemur_mircen_2019_t2(data_dir=tst.tmpdir,
                                                 verbose=0,
                                                 url=local_url,
                                                 subjects=[0, 10, 33])

    # subject 0 has already been downloaded
    assert_equal(len(tst.mock_url_request.urls), 4)
    assert_equal(len(lemur_t2.anat), 3)

    # returned phenotypic data will be an array
    assert_true(isinstance(lemur_t2.pheno, np.recarray))
    np.testing.assert_array_equal(lemur_t2.pheno.animal_id,
                                  ['"sub-01"', '"sub-11"', '"sub-34"'])
    assert_not_equal(lemur_t2.description, '')    
