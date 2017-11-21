import os
import shutil
import nibabel
from sammba.externals.nipype.caching import Memory
from sammba.externals.nipype.utils.filemanip import fname_presuffix
import sammba.externals.nipype.pipeline.engine as pe
from sammba.externals.nipype.interfaces import afni, fsl
from sammba.interfaces import RatsMM


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
    header = img.header
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
                  caching_dir=None, clear_memory=False, overwrite=False):
    if caching:
        memory = Memory(caching_dir)

    if caching_dir is None:
        caching_dir = os.getcwd()

    if overwrite:
        environ = {'AFNI_DECONFLICT': 'OVERWRITE'}
    else:
        environ = {}

    if caching:
        copy = memory.cache(afni.Copy)
    else:
        copy = afni.Copy().run

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

    if caching:
        refit = memory.cache(afni.Refit)
    else:
        refit = afni.Refit().run

    out_refit = refit(in_file=orig_to_fix_filename,
                      atrcopy=(orig_reference_filename,
                               'IJK_TO_DICOM_REAL'))

    out_copy = copy(in_file=out_refit.outputs.out_file,
                    environ={'AFNI_DECONFLICT': 'OVERWRITE'},
                    out_file=to_fix_filename)

    if clear_memory:
        shutil.rmtree(tmp_folder)
        memory.clear_previous_run()

    return out_copy.outputs.out_file


def create_pipeline_graph(pipeline_name, graph_file,
                          graph_kind='hierarchical'):
    """
    Parameters
    ----------
    pipeline_name : one of {'rigid-body_registration',
                            'affine registration',
                            'nonlinear registration'}

    graph_file : Str, path to save the graph image

    graph_kind : one of {'orig', 'hierarchical', 'flat', 'exec', 'colored'}
        'orig': creates a top level graph without expanding internal
        workflow nodes;
        'flat': expands workflow nodes recursively;
        'hierarchical': expands workflow nodes recursively with a notion
        on hierarchy;
        'colored': expands workflow nodes recursively with a
        notion on hierarchy in color;
        'exec': expands workflows to depict iterables
    """
    pipeline_names = ['rigid_anat_to_common', 'affine_anat_to_common',
                      'nonlinear_anat_to_common']
    if pipeline_name not in pipeline_names:
        raise NotImplementedError(
            'Pipeline name must be one of {0}, you entered {1}'.format(
                pipeline_names, pipeline_name))
    graph_kinds = ['orig', 'hierarchical', 'flat', 'exec', 'colored']
    if graph_kind not in graph_kinds:
        raise ValueError(
            'Graph kind must be one of {0}, you entered {1}'.format(
                graph_kinds, graph_kind))

    graph_file = os.path.abspath(graph_file)
    write_dir = os.path.dirname(graph_file)
    if not os.path.isdir(write_dir):
        raise IOError('{0} directory not existant'.format(write_dir))
    workflow = pe.Workflow(name=pipeline_name)

    #######################################################################
    # Specify rigid body registration pipeline steps
    unifize_node = pe.Node(interface=afni.Unifize(), name='bias_correct')
    clip_level_node = pe.Node(interface=afni.ClipLevel(),
                              name='compute_mask_threshold')
    rats_node = pe.Node(interface=RatsMM(), name='compute_brain_mask')
    apply_mask_node = pe.Node(interface=fsl.ApplyMask(),
                              name='apply_brain_mask')
    center_mass_node = pe.Node(interface=afni.CenterMass(),
                               name='compute_and_set_cm_in_header')
    refit_copy_node = pe.Node(afni.Refit(), name='copy_cm_in_header')
    tcat_node1 = pe.Node(afni.TCat(), name='concatenate_across_individuals1')
    tstat_node1 = pe.Node(afni.TStat(), name='compute_average1')
    undump_node = pe.Node(afni.Undump(), name='create_empty_template')
    refit_set_node = pe.Node(afni.Refit(), name='set_cm_in_header')
    resample_node1 = pe.Node(afni.Resample(), name='resample1')
    resample_node2 = pe.Node(afni.Resample(), name='resample2')
    shift_rotate_node = pe.Node(afni.Allineate(), name='shift_rotate')
    apply_allineate_node1 = pe.Node(afni.Allineate(),
                                    name='apply_transform1')
    tcat_node2 = pe.Node(afni.TCat(), name='concatenate_across_individuals2')
    tstat_node2 = pe.Node(afni.TStat(), name='compute_average2')
    tcat_node3 = pe.Node(afni.TCat(), name='concatenate_across_individuals3')
    tstat_node3 = pe.Node(afni.TStat(), name='compute_average3')

    workflow.add_nodes([unifize_node, clip_level_node, rats_node,
                        apply_mask_node, center_mass_node, refit_copy_node,
                        tcat_node1, tstat_node1,
                        undump_node, refit_set_node, resample_node1,
                        resample_node2,
                        shift_rotate_node, apply_allineate_node1,
                        tcat_node2, tstat_node2, tcat_node3, tstat_node3])

    #######################################################################
    # and connections
    workflow.connect(unifize_node, 'out_file', clip_level_node, 'in_file')
    workflow.connect(clip_level_node, 'clip_val',
                     rats_node, 'intensity_threshold')
    workflow.connect(unifize_node, 'out_file', rats_node, 'in_file')
    workflow.connect(rats_node, 'out_file', apply_mask_node, 'mask_file')
    workflow.connect(apply_mask_node, 'out_file',
                     center_mass_node, 'in_file')
    workflow.connect(unifize_node, 'out_file', refit_copy_node, 'in_file')
    workflow.connect(center_mass_node, 'out_file',
                     refit_copy_node, 'duporigin_file')
    workflow.connect(center_mass_node, 'out_file', tcat_node1, 'in_files')
    workflow.connect(tcat_node1, 'out_file', tstat_node1, 'in_file')
    workflow.connect(tstat_node1, 'out_file', undump_node, 'in_file')
    workflow.connect(undump_node, 'out_file', refit_set_node, 'in_file')
    workflow.connect(refit_set_node, 'out_file', resample_node1, 'master')
    workflow.connect(refit_set_node, 'out_file', resample_node1, 'in_file')
    workflow.connect(refit_copy_node, 'out_file', resample_node2, 'master')
    workflow.connect(center_mass_node, 'out_file', resample_node2, 'in_file')
    workflow.connect(resample_node2, 'out_file', tcat_node2, 'in_files')
    workflow.connect(tcat_node2, 'out_file', tstat_node2, 'in_file')
    workflow.connect(tstat_node2, 'out_file', shift_rotate_node, 'reference')
    workflow.connect(resample_node2, 'out_file', shift_rotate_node, 'in_file')
    workflow.connect(tstat_node2, 'out_file', apply_allineate_node1, 'master')
    workflow.connect(resample_node1, 'out_file',
                     apply_allineate_node1, 'in_file')
    workflow.connect(shift_rotate_node, 'out_matrix',
                     apply_allineate_node1, 'in_matrix')
    workflow.connect(apply_allineate_node1, 'out_file', tcat_node3, 'in_files')
    workflow.connect(tcat_node3, 'out_file', tstat_node3, 'in_file')
    if pipeline_name in ['affine_registration', 'nonlinear_registration']:
        mask_node = pe.Node(afni.MaskTool(), name='generate_count_mask')
        allineate_node3 = pe.Node(afni.Allineate(), name='affine_register')
        catmatvec_node = pe.Node(afni.CatMatvec(), name='concatenate_transforms')
        allineate_node4 = pe.Node(afni.Allineate(), name='apply_transform_on_head')
        allineate_node5 = pe.Node(afni.Allineate(),
                                  name='apply_transform_on_brain')
        tcat_node3 = pe.Node(
            afni.TCat(), name='concatenate_affine_registred_across_individuals')
        tstat_node3 = pe.Node(afni.TStat(),
                              name='compute_affine_registred_average')
    
        workflow.add_nodes([mask_node, allineate_node3, catmatvec_node,
                            allineate_node4, allineate_node5,
                            tcat_node3, tstat_node3])

    if pipeline_name == 'nonlinear_anat_to_common':
        pass

    current_dir = os.getcwd()
    os.chdir(write_dir)
    _ = workflow.write_graph(graph2use=graph_kind, dotfilename=graph_file)
    os.chdir(current_dir)
