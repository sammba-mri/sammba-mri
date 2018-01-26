import sys
sys.path.remove('/usr/lib/python2.7/dist-packages')
sys.path.append('/usr/lib/python2.7/dist-packages')
import os
import shutil
import tempfile
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from nose.tools import assert_true, assert_false
import nibabel
from nilearn._utils.testing import assert_raises_regex
from sammba.externals.nipype.interfaces import afni, fsl
from sammba.registration import utils
from sammba import testing_data


def test_reset_affines():
    # test qform is replaced by sfrom with default parameters
    in_file = os.path.join(os.path.dirname(testing_data.__file__),
                           'anat.nii.gz')
    source_file = os.path.join(os.path.dirname(testing_data.__file__),
                               'anat.nii.gz')
    tempdir = tempfile.mkdtemp()
    in_file = os.path.join(tempdir, 'in_anat.nii.gz')
    shutil.copy(source_file, in_file)
    in_header = nibabel.load(in_file).header
    in_sform, in_scode = in_header.get_sform(coded=True)
    out_file = os.path.join(tempdir, 'out_anat.nii.gz')
    utils._reset_affines(in_file, out_file)
    out_header = nibabel.load(out_file).header
    out_sform, out_scode = out_header.get_sform(coded=True)
    out_qform, out_qcode = out_header.get_qform(coded=True)

    assert_array_equal(out_sform, in_sform)
    assert_array_almost_equal(out_qform, out_sform)
    assert_array_equal(out_qcode, out_scode)

    # same test, overwriting
    utils._reset_affines(in_file, out_file, overwrite=True)
    out_header = nibabel.load(out_file).header
    assert_array_equal(out_sform, in_sform)
    assert_array_almost_equal(out_qform, out_sform)
    assert_array_equal(out_qcode, out_scode)

    # previous test, with verbose=0
    utils._reset_affines(in_file, out_file, overwrite=True, verbose=0)
    out_header = nibabel.load(out_file).header
    assert_array_equal(out_sform, in_sform)
    assert_array_almost_equal(out_qform, out_sform)
    assert_array_equal(out_qcode, out_scode)

    if os.path.exists(in_file):
        os.remove(in_file)
    if os.path.exists(out_file):
        os.remove(out_file)

    if os.path.exists(tempdir):
        os.removedirs(tempdir)


def test_check_same_obliquity():
    img_filename1 = os.path.join(os.path.dirname(testing_data.__file__),
                                 'anat.nii.gz')
    img_filename2 = os.path.join(os.path.dirname(testing_data.__file__),
                                 'func.nii.gz')
    assert_true(utils._check_same_obliquity(img_filename1, img_filename1))
    assert_false(utils._check_same_obliquity(img_filename1, img_filename2))


def test_fix_obliquity():
    target_filename = os.path.join(os.path.dirname(testing_data.__file__),
                                   'anat.nii.gz')
    if afni.Info().version():
        # ZCutUp removes the obliquity
        tempdir = tempfile.mkdtemp()
        tmp_filename = os.path.join(tempdir, 'img_test_obliquity.nii.gz')
        slicer = afni.ZCutUp().run
        _ = slicer(in_file=target_filename,
                   keep='0 27',
                   out_file=tmp_filename)
        assert_false(
            utils._check_same_obliquity(tmp_filename, target_filename))
        utils.fix_obliquity(tmp_filename, target_filename)
        assert_true(
            utils._check_same_obliquity(tmp_filename, target_filename))

        if os.path.exists(tmp_filename):
            os.remove(tmp_filename)
        if os.path.exists(tempdir):
            os.removedirs(tempdir)


def test_copy_geometry():
    func_filename = os.path.join(os.path.dirname(testing_data.__file__),
                                 'func.nii.gz')
    if afni.Info().version():
        # Allineate changes the geometry
        tempdir = tempfile.mkdtemp()
        tmp_filename = os.path.join(tempdir, 'test_geometry.1D')
        volreg = afni.Volreg().run
        out_volreg = volreg(
            in_file=func_filename,
            outputtype='NIFTI_GZ',
            oned_file=os.path.join(tempdir, 'test_geometry.1Dfile.1D'),
            oned_matrix_save=os.path.join(tempdir, 'test_geometry.aff12.1D'))
        allineate = afni.Allineate().run
        out_allineate = allineate(
            in_file=func_filename,
            master=func_filename,
            in_matrix=out_volreg.outputs.oned_matrix_save,
            out_file=os.path.join(tempdir, 'test_geometry.nii.gz'))

        header = nibabel.load(out_volreg.outputs.out_file).header
        allineate_header = nibabel.load(out_allineate.outputs.out_file).header
        allineate_equal_fields = [np.allclose(v1, v2)
                                  for v1, v2 in zip(header.values(),
                                                    allineate_header.values())]
        assert_false(np.alltrue(allineate_equal_fields))
        changed_filename = utils.copy_geometry(
            out_volreg.outputs.out_file, out_allineate.outputs.out_file,
            in_place=False)
        copy_geom = fsl.CopyGeom().run
        out_copy_geom = copy_geom(dest_file=out_allineate.outputs.out_file,
                                  in_file=out_volreg.outputs.out_file)
        changed_header = nibabel.load(changed_filename).header
        fsl_cpgeom_header = nibabel.load(out_copy_geom.outputs.out_file).header
        copied_equal_fields = [np.allclose(v1, v2)
                               for v1, v2 in zip(changed_header.values(),
                                                 fsl_cpgeom_header.values())]
        assert(np.alltrue(copied_equal_fields))

        if os.path.exists(tmp_filename):
            os.remove(tmp_filename)
        if os.path.exists(tempdir):
            os.removedirs(tempdir)


def test_create_pipeline_graph():
    # Check error is raised if unkown pipeline name
    assert_raises_regex(NotImplementedError, 'Pipeline name must be one of',
                        utils.create_pipeline_graph, 'rigid-body', '')

    # Check error is raised if wrong graph kind
    assert_raises_regex(ValueError, "Graph kind must be one of ",
                        utils.create_pipeline_graph,
                        'anats_to_common_rigid', '', graph_kind='original')

    # Check graph file is created
    tempdir = tempfile.mkdtemp()
    graph_file = os.path.join(tempdir, 'tmp_graph.png')
    utils.create_pipeline_graph('anats_to_common_rigid', graph_file)
    assert(os.path.exists(graph_file))

    graph_file_root, _ = os.path.splitext(graph_file)
    if os.path.exists(graph_file):
        os.remove(graph_file)
    if os.path.exists(graph_file_root):
        os.remove(graph_file_root)
    if os.path.exists(tempdir):
        os.removedirs(tempdir)
