import os
import tempfile
from numpy.testing import assert_array_equal
import nibabel
from nilearn._utils.testing import assert_raises_regex
from sammba.preprocessors import utils
from sammba import testing_data


def test_reset_affines():
    in_file = os.path.join(os.path.dirname(testing_data.__file__),
                           'func.nii.gz')
    in_header = nibabel.load(in_file).header
    in_sform, in_scode = in_header.get_sform(coded=True)
    tempdir = tempfile.mkdtemp()
    out_file = os.path.join(tempdir, 'func.nii.gz')
    utils._reset_affines(in_file, out_file)
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


def test_create_pipeline_graph():
    # Check error is raised if unkown pipeline name
    assert_raises_regex(NotImplementedError, 'Pipeline name must be one of',
                        utils.create_pipeline_graph, 'rigid-body', '')

    # Check error is raised if wrong graph kind
    assert_raises_regex(ValueError, "Graph kind must be one of ",
                        utils.create_pipeline_graph,
                        'rigid-body_registration', '', graph_kind='original')

    # Check graph file is created
    tempdir = tempfile.mkdtemp()
    graph_file = os.path.join(tempdir, 'tmp_graph.dot')
    utils.create_pipeline_graph('rigid-body_registration', graph_file)
    assert(os.path.exists(graph_file))

    os.remove(graph_file)
    os.removedirs(tempdir)

    # Check error is raised if writing directory dosn't exist
    assert_raises_regex(IOError, 'directory not existant',
                        utils.create_pipeline_graph,
                        'rigid-body_registration', graph_file)