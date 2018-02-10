import os
import numpy as np
import json
from nose import with_setup
from nose.tools import assert_equal, assert_true, assert_not_equal
from nilearn.datasets import utils
from nilearn.datasets.tests import test_utils as tst

from sammba.data_fetchers import func


def setup_mock():
    return tst.setup_mock(utils, func)


def teardown_mock():
    return tst.teardown_mock(utils, func)


@with_setup(setup_mock, teardown_mock)
@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_fetch_zurich_test_retest():
    local_url = "file://" + tst.datadir
    subjects_ids = ['1366', '1367', '1368', '1369', '1371', '1378', '1380',
                    '1402', '1403', '1404', '1405', '1406', '1407', '1411',
                    '1412']
    baseline_ids = [s + '_baseline' for s in subjects_ids]
    post_ids = [s + '_post' for s in subjects_ids]
    result = [{u'insert_date': u'2007-10-04 12:23:04.0',
               u'xsiType': u'xnat:mrSessionData',
               u'URI': u'/data/experiments/CENTRAL_E{0:05f}'.format(n),
               u'label': s,
               u'project': u'CSD_MRI_MOUSE', u'date': u'',
               u'ID': u'CENTRAL_E{0:05f}'.format(n)}
              for n, s in enumerate(baseline_ids + post_ids)]

    urls_table = {'ResultSet': {'totalRecords': '3982',
                                'Result': result,
                                'title': 'Matching experiments'}}
    # Create a dummy urls file
    zurich_dir = os.path.join(tst.tmpdir, 'zurich_retest')
    os.mkdir(zurich_dir)
    json_file = os.path.join(zurich_dir, 'experiments.html')
    with open(json_file, 'wb', encoding='utf8') as f:
        json.dump(urls_table, f)

    # First session, all subjects
    zurich = func.fetch_zurich_test_retest(data_dir=tst.tmpdir, url=local_url,
                                           verbose=0)
    assert_equal(len(tst.mock_url_request.urls), 30)
    assert_equal(len(zurich.func), 15)
    assert_equal(len(zurich.anat), 15)
    assert_true(np.all(np.asarray(zurich.session) == 1))

    # Both sessions, 12 subjects
    tst.mock_url_request.reset()
    zurich = func.fetch_zurich_test_retest(data_dir=tst.tmpdir, verbose=0,
                                           sessions=[1, 2], url=local_url,
                                           subjects=range(12))

    # Session 1 has already been downloaded
    assert_equal(len(tst.mock_url_request.urls), 24)
    assert_equal(len(zurich.func), 24)
    assert_equal(len(zurich.anat), 24)
    s = np.asarray(zurich.session)
    assert_true(np.all(s[:12] == 1))
    assert_true(np.all(s[12:24] == 2))
    assert_not_equal(zurich.description, '')


@with_setup(setup_mock, teardown_mock)
@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_fetch_zurich_anesthesiant():
    local_url = "file://" + tst.datadir
    subjects_ids = ['iso2273', 'iso2274', 'iso2238', 'iso2239', 'iso2250',
                    'iso2270', 'iso3294', 'iso3296', 'med2259', 'med2241',
                    'med2247', 'med2251', 'med2256', 'med2257', 'mi272871',
                    'mi273299', 'mi273457', 'mi273458', 'mi273459',
                    'mi273460', 'mi273461', 'mi273300', 'medHalfDose',
                    'medHalfDose1', 'medHalfDose2', 'medHalfDose3',
                    'iso1c3perc', 'iso1c3perc', 'iso1c5perc',
                    'iso2870_1c5perc']
    result = [{u'insert_date': u'2007-10-04 12:23:04.0',
               u'xsiType': u'xnat:mrSessionData',
               u'URI': u'/data/experiments/CENTRAL_E{0:05f}'.format(n),
               u'label': s,
               u'project': u'fMRI_ane_mouse', u'date': u'',
               u'ID': u'CENTRAL_E{0:05f}'.format(n)}
              for n, s in enumerate(subjects_ids)]

    urls_table = {'ResultSet': {'totalRecords': '3982',
                                'Result': result,
                                'title': 'Matching experiments'}}
    # Create a dummy urls file
    anest_dir = os.path.join(tst.tmpdir, 'zurich_anest')
    os.mkdir(anest_dir)
    json_file = os.path.join(anest_dir, 'experiments.html')
    with open(json_file, 'wb', encoding='utf8') as f:
        json.dump(urls_table, f)

    # 1 subject
    anest = func.fetch_zurich_anesthesiant(subjects=[0], url=local_url,
                                           data_dir=tst.tmpdir, verbose=0)
    assert_equal(len(tst.mock_url_request.urls), 2)
    assert_equal(len(anest.func), 1)
    assert_equal(len(anest.anesthesiant), 1)
    assert_true(anest.anesthesiant[0] == 'Iso1')

    # First subject has already been downloaded
    anest = func.fetch_zurich_anesthesiant(data_dir=tst.tmpdir, verbose=0)
    assert_equal(len(tst.mock_url_request.urls), 58)
    assert_equal(len(anest.func), 30)
    s = np.asarray(anest.anesthesiant)
    assert_true(np.all(s[:8] == 'Iso1'))
    assert_true(np.all(s[8:14] == 'Med'))
    assert_true(np.all(s[14:22] == 'Med-Iso'))
    assert_true(np.all(s[22:26] == 'Med-half'))
    assert_true(np.all(s[26:28] == 'Iso1pt3'))
    assert_true(np.all(s[28:] == 'Iso1pt5'))
    assert_not_equal(anest.description, '')
