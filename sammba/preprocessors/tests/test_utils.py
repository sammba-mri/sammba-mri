import os
from numpy.testing import assert_array_equal
from sammba.preprocessors import utils
from sammba import testing_data
import nibabel
import tempfile


def test_correct_affines():
    in_file = os.path.join(os.path.dirname(testing_data.__file__),
                           'func.nii.gz')
    in_header = nibabel.load(in_file).header
    in_sform, in_scode = in_header.get_sform(coded=True)
    tempdir = tempfile.mkdtemp()
    out_file = os.path.join(tempdir, 'func.nii.gz')
    utils.correct_affines(in_file, out_file)
    out_header = nibabel.load(in_file).header
    out_sform, out_scode = out_header.get_sform(coded=True)
    out_qform, out_qcode = out_header.get_qform(coded=True)

    assert_array_equal(out_sform, in_sform)
    assert_array_equal(out_qform, out_sform)
    assert_array_equal(out_qcode, out_scode)

    if os.path.exists(out_file):
        os.remove(out_file)
    if os.path.exists(tempdir):
        os.removedirs(tempdir)
