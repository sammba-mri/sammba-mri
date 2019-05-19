import os
from nose import with_setup
from nilearn.datasets.tests import test_utils as tst
from nilearn._utils.testing import assert_raises_regex
from sammba import graphs


@with_setup(tst.setup_tmpdata, tst.teardown_tmpdata)
def test_create_pipeline_graph():
    # Check error is raised if unkown pipeline name
    assert_raises_regex(NotImplementedError, 'Pipeline name must be one of',
                        graphs.create_pipeline_graph, 'rigid-body', '')

    # Check error is raised if wrong graph kind
    assert_raises_regex(ValueError, "Graph kind must be one of ",
                        graphs.create_pipeline_graph,
                        'anats_to_common_rigid', '', graph_kind='original')

    # Check graph file is created
    graph_file = os.path.join(tst.tmpdir, 'tmp_graph.png')
    graphs.create_pipeline_graph('anats_to_common_rigid', graph_file)
    assert(os.path.exists(graph_file))
