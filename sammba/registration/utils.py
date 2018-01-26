import os
import shutil
import numpy as np
import nibabel
from sammba.externals.nipype.caching import Memory
from sammba.externals.nipype.utils.filemanip import fname_presuffix
import sammba.externals.nipype.pipeline.engine as pe
from sammba.externals.nipype.interfaces import afni, fsl
from sammba.interfaces import segmentation


def _reset_affines(in_file, out_file, overwrite=False, axes_to_permute=None,
                   axes_to_flip=None, xyzscale=None, center_mass=None,
                   verbose=1):
    """Sets the qform equal to the sform in the header, with optionally
       rescaling, setting the image center to 0, permuting or/and flipping axes
    """
    if not os.path.isfile(out_file) or overwrite:
        shutil.copy(in_file, out_file)
    else:
        return

    if verbose:
        terminal_output = 'stream'
    else:
        terminal_output = 'none'

    in_file = out_file
    if xyzscale is not None:
        refit = afni.Refit()
        refit.inputs.in_file = in_file
        refit.inputs.xyzscale = xyzscale
        refit.set_default_terminal_output(terminal_output)
        result = refit.run()
        in_file = result.outputs.out_file

    if center_mass is not None:
        set_center_mass = afni.CenterMass()
        set_center_mass.inputs.in_file = in_file
        set_center_mass.inputs.cm_file = fname_presuffix(out_file,
                                                         suffix='.txt',
                                                         use_ext=False)
        set_center_mass.inputs.set_cm = center_mass
#        set_center_mass.set_default_terminal_output(terminal_output) # XXX BUG
        result = set_center_mass.run()
        in_file = result.outputs.out_file

    img = nibabel.load(in_file)
    header = img.header.copy()
    sform, code = header.get_sform(coded=True)

    if axes_to_flip:
        for axis in axes_to_flip:
            sform[axis] *= -1

    if axes_to_permute:
        for (axis1, axis2) in axes_to_permute:
            sform[[axis1, axis2]] = sform[[axis2, axis1]]

    header.set_sform(sform)
    header.set_qform(sform, int(code))
    nibabel.Nifti1Image(img.get_data(), sform, header).to_filename(out_file)


def fix_obliquity(to_fix_filename, reference_filename, caching=False,
                  caching_dir=None, overwrite=False,
                  verbose=True):
    if caching_dir is None:
        caching_dir = os.getcwd()
    if caching:
        memory = Memory(caching_dir)

    if verbose:
        terminal_output = 'allatonce'
    else:
        terminal_output = 'none'

    if overwrite:
        environ = {'AFNI_DECONFLICT': 'OVERWRITE'}
    else:
        environ = {}

    if caching:
        copy = memory.cache(afni.Copy)
        refit = memory.cache(afni.Refit)
        copy.interface().set_default_terminal_output(terminal_output)
        refit.interface().set_default_terminal_output(terminal_output)
    else:
        copy = afni.Copy(terminal_output=terminal_output).run
        refit = afni.Refit(terminal_output=terminal_output).run

    tmp_folder = os.path.join(caching_dir, 'tmp')
    if not os.path.isdir(tmp_folder):
        os.makedirs(tmp_folder)

    reference_basename = os.path.basename(reference_filename)
    orig_reference_filename = fname_presuffix(os.path.join(
        tmp_folder, reference_basename), suffix='+orig.BRIK',
        use_ext=False)
    if not os.path.isfile(orig_reference_filename) or overwrite:
        out_copy_oblique = copy(in_file=reference_filename,
                                out_file=orig_reference_filename,
                                environ=environ)
        orig_reference_filename = out_copy_oblique.outputs.out_file

    to_fix_basename = os.path.basename(to_fix_filename)
    orig_to_fix_filename = fname_presuffix(os.path.join(
        tmp_folder, to_fix_basename), suffix='+orig.BRIK',
        use_ext=False)
    if not os.path.isfile(orig_to_fix_filename) or overwrite:
        out_copy = copy(in_file=to_fix_filename,
                        out_file=orig_to_fix_filename,
                        environ=environ)
        orig_to_fix_filename = out_copy.outputs.out_file

    out_refit = refit(in_file=orig_to_fix_filename,
                      atrcopy=(orig_reference_filename,
                               'IJK_TO_DICOM_REAL'))

    out_copy = copy(in_file=out_refit.outputs.out_file,
                    environ={'AFNI_DECONFLICT': 'OVERWRITE'},
                    out_file=to_fix_filename)

    if overwrite:
        shutil.rmtree(tmp_folder)
        if caching:
            memory.clear_previous_run()

    return out_copy.outputs.out_file


def _check_same_obliquity(img_filename1, img_filename2):
    headers_values = []
    for img_filename in [img_filename1, img_filename2]:
        img = nibabel.load(img_filename)
        header = img.header
        header_values = [header['pixdim'][:4]]
        # TODO: check that 'qform_code', 'sform_code' can differ
        for key in ['quatern_b', 'quatern_c', 'quatern_d', 'qoffset_x',
                    'qoffset_y', 'qoffset_z', 'srow_x', 'srow_y', 'srow_z']:
            header_values.append(header[key])
        headers_values.append(header_values)
    equal_fields = [np.allclose(v1, v2)
                    for v1, v2 in zip(headers_values[0], headers_values[1])]

    return np.alltrue(equal_fields)


def copy_geometry(filename_to_copy, filename_to_change, out_filename=None,
                  copy_shape=True, in_place=True):
    """ Mimics FSL command fslcpgeom to copy geometry information from header.

    filename_to_copy : str
        Path to the image with the header information to copy.

    to_change_filename : str
        Path to the image with the new geometry.

    out_filename : str or None, optional
        Path to save the image with the changed header to.

    copy_shape : bool, optional
        If False, image data shape is not copied.

    in_place : bool, optional
        If False, a new image is created with the copied geometry information.
    """
    img_to_copy = nibabel.load(filename_to_copy)
    img_to_change = nibabel.load(filename_to_change)
    header_to_copy = img_to_copy.header
    header_to_change = img_to_change.header
    new_header = header_to_change.copy()
    geometry_keys = ['sform_code', 'qform_code', 'quatern_b',
                     'quatern_c', 'quatern_d', 'qoffset_x',
                     'qoffset_y', 'qoffset_z', 'srow_x', 'srow_y', 'srow_z']
    if copy_shape:
        geometry_keys += ['dim']
    for key in geometry_keys:
        new_header[key] = header_to_copy[key]

    new_header['pixdim'][:4] = header_to_copy['pixdim'][:4]
    new_img = nibabel.Nifti1Image(img_to_change.get_data(), img_to_copy.affine,
                                  header=new_header)

    if not in_place:
        if out_filename is None:
            out_filename = fname_presuffix(filename_to_change,
                                           suffix='copied_geom')
    else:
        out_filename = filename_to_change
    new_img.to_filename(out_filename)
    return out_filename


def check_same_geometry(header1, header2):
    unchecked_fields = ['descrip', 'pixdim', 'scl_slope', 'scl_inter',
                        'xyzt_units', 'qform_code', 'regular']
    equal_values = []
    for item1, item2 in zip(header1.items(), header2.items()):
        if item1[0] in unchecked_fields:
            continue
        if item1[1].dtype in [np.dtype('float'), np.dtype('float32'),
                              np.dtype('float64'), np.dtype('float128')]:
            equal_values.append(np.allclose(item1[1], item2[1]))
        else:
            equal_values.append(item1[1] == item2[1])

    equal_values.append(np.allclose(header1['pixdim'][:4],
                                    header2['pixdim'][:4]))
    return [np.alltrue(a) for a in equal_values]


def create_pipeline_graph(pipeline_name, graph_file,
                          graph_kind='hierarchical'):
    """Creates pipeline graph for a given piepline.

    Parameters
    ----------
    pipeline_name : one of {'anat_to_common_rigid', 'anat_to_common_affine',
        'anat_to_common_nonlinear'}
        Pipeline name.

    graph_file : str.
        Path to save the graph image to.

    graph_kind : one of {'orig', 'hierarchical', 'flat', 'exec', 'colored'}, optional.
        The kind of the graph, passed to
        sammba.externals.nipype.pipeline.workflows.Workflow().write_graph
    """
    pipeline_names = ['anats_to_common_rigid', 'anats_to_common_affine',
                      'anats_to_common_nonlinear']
    if pipeline_name not in pipeline_names:
        raise NotImplementedError(
            'Pipeline name must be one of {0}, you entered {1}'.format(
                pipeline_names, pipeline_name))
    graph_kinds = ['orig', 'hierarchical', 'flat', 'exec', 'colored']
    if graph_kind not in graph_kinds:
        raise ValueError(
            'Graph kind must be one of {0}, you entered {1}'.format(
                graph_kinds, graph_kind))

    workflow = pe.Workflow(name=pipeline_name)

    #######################################################################
    # Specify rigid body registration pipeline steps
    unifize = pe.Node(interface=afni.Unifize(), name='bias_correct')
    clip_level = pe.Node(interface=afni.ClipLevel(),
                         name='compute_mask_threshold')
    compute_mask = pe.Node(interface=segmentation.MathMorphoMask(),
                           name='compute_brain_mask')
    apply_mask = pe.Node(interface=fsl.ApplyMask(), name='apply_brain_mask')
    center_mass = pe.Node(interface=afni.CenterMass(),
                          name='compute_and_set_cm_in_header')
    refit_copy = pe.Node(afni.Refit(), name='copy_cm_in_header')
    tcat1 = pe.Node(afni.TCat(), name='concatenate_across_individuals1')
    tstat1 = pe.Node(afni.TStat(), name='compute_average1')
    undump = pe.Node(afni.Undump(), name='create_empty_template')
    refit_set = pe.Node(afni.Refit(), name='set_cm_in_header')
    resample1 = pe.Node(afni.Resample(), name='resample1')
    resample2 = pe.Node(afni.Resample(), name='resample2')
    shift_rotate = pe.Node(afni.Allineate(), name='shift_rotate')
    apply_allineate1 = pe.Node(afni.Allineate(), name='apply_transform1')
    tcat2 = pe.Node(afni.TCat(), name='concatenate_across_individuals2')
    tstat2 = pe.Node(afni.TStat(), name='compute_average2')
    tcat3 = pe.Node(afni.TCat(), name='concatenate_across_individuals3')
    tstat3 = pe.Node(afni.TStat(), name='compute_average3')

    workflow.add_nodes([unifize, clip_level, compute_mask, apply_mask,
                        center_mass,
                        refit_copy, tcat1, tstat1, undump, refit_set,
                        resample1, resample2, shift_rotate, apply_allineate1,
                        tcat2, tstat2, tcat3, tstat3])

    #######################################################################
    # and connections
    workflow.connect(unifize, 'out_file', clip_level, 'in_file')
    workflow.connect(clip_level, 'clip_val',
                     compute_mask, 'intensity_threshold')
    workflow.connect(unifize, 'out_file', compute_mask, 'in_file')
    workflow.connect(compute_mask, 'out_file', apply_mask, 'mask_file')
    workflow.connect(apply_mask, 'out_file',
                     center_mass, 'in_file')
    workflow.connect(unifize, 'out_file', refit_copy, 'in_file')
    workflow.connect(center_mass, 'out_file',
                     refit_copy, 'duporigin_file')
    workflow.connect(center_mass, 'out_file', tcat1, 'in_files')
    workflow.connect(tcat1, 'out_file', tstat1, 'in_file')
    workflow.connect(tstat1, 'out_file', undump, 'in_file')
    workflow.connect(undump, 'out_file', refit_set, 'in_file')
    workflow.connect(refit_set, 'out_file', resample1, 'master')
    workflow.connect(refit_copy, 'out_file', resample1, 'in_file')
    workflow.connect(refit_set, 'out_file', resample2, 'master')
    workflow.connect(center_mass, 'out_file', resample2, 'in_file')
    workflow.connect(resample2, 'out_file', tcat2, 'in_files')
    workflow.connect(tcat2, 'out_file', tstat2, 'in_file')
    workflow.connect(tstat2, 'out_file', shift_rotate, 'reference')
    workflow.connect(resample2, 'out_file', shift_rotate, 'in_file')
    workflow.connect(tstat2, 'out_file', apply_allineate1, 'master')
    workflow.connect(resample1, 'out_file',
                     apply_allineate1, 'in_file')
    workflow.connect(shift_rotate, 'out_matrix',
                     apply_allineate1, 'in_matrix')
    workflow.connect(apply_allineate1, 'out_file', tcat3, 'in_files')
    workflow.connect(tcat3, 'out_file', tstat3, 'in_file')
    if pipeline_name in ['anats_to_common_affine',
                         'anat_to_common_nonlinear']:
        mask = pe.Node(afni.MaskTool(), name='generate_count_mask')
        allineate = pe.Node(afni.Allineate(), name='allineate')
        catmatvec = pe.Node(afni.CatMatvec(), name='concatenate_transforms')
        apply_allineate2 = pe.Node(afni.Allineate(), name='apply_transform2')
        tcat3 = pe.Node(
            afni.TCat(), name='concatenate_across_individuals4')
        tstat3 = pe.Node(afni.TStat(), name='compute_average4')

        workflow.add_nodes([mask, allineate, catmatvec, apply_allineate2,
                            tcat3, tstat3])

        workflow.connect(tcat2, 'out_file', mask, 'in_file')
        workflow.connect(mask, 'out_file', allineate, 'weight')
        workflow.connect(apply_allineate1, 'out_file',
                         allineate, 'in_file')
        workflow.connect(allineate, 'out_matrix',
                         catmatvec, 'in_file')
        #XXX how can we enter multiple files ? 
        workflow.connect(catmatvec, 'out_file',
                         apply_allineate2, 'in_matrix')
        workflow.connect(resample1, 'out_file',
                         apply_allineate2, 'in_file')
        workflow.connect(apply_allineate2, 'out_file', tcat3, 'in_files')
        workflow.connect(tcat3, 'out_file', tstat3, 'in_file')

    if pipeline_name == 'anats_to_common_nonlinear':
        pass

    graph_file_root, graph_file_ext = os.path.splitext(graph_file)
    if graph_file_ext:
        _ = workflow.write_graph(graph2use=graph_kind,
                                 format=graph_file_ext[1:],
                                 dotfilename=graph_file_root)
    else:
        _ = workflow.write_graph(graph2use=graph_kind,
                                 dotfilename=graph_file_root)
