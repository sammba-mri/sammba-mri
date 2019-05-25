import os
from nipype.interfaces import afni
import nipype.pipeline.engine as pe
from sammba.segmentation import interfaces


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
        nipype.pipeline.workflows.Workflow().write_graph
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
    compute_mask = pe.Node(interface=interfaces.MathMorphoMask(),
                           name='compute_brain_mask')
    apply_mask = pe.Node(interface=afni.Calc(), name='apply_brain_mask')
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
    workflow.connect(compute_mask, 'out_file', apply_mask, 'in_file_a')
    workflow.connect(unifize, 'out_file', apply_mask, 'in_file_b')
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
