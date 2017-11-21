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
    unifize = pe.Node(interface=afni.Unifize(), name='bias_correct')
    clip_level = pe.Node(interface=afni.ClipLevel(),
                         name='compute_mask_threshold')
    rats = pe.Node(interface=RatsMM(), name='compute_brain_mask')
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

    workflow.adds([unifize, clip_level, rats, apply_mask, center_mass,
                   refit_copy, tcat1, tstat1, undump, refit_set, resample1,
                   resample2, shift_rotate, apply_allineate1, tcat2, tstat2,
                   tcat3, tstat3])

    #######################################################################
    # and connections
    workflow.connect(unifize, 'out_file', clip_level, 'in_file')
    workflow.connect(clip_level, 'clip_val',
                     rats, 'intensity_threshold')
    workflow.connect(unifize, 'out_file', rats, 'in_file')
    workflow.connect(rats, 'out_file', apply_mask, 'mask_file')
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
    if pipeline_name in ['affine_registration', 'nonlinear_registration']:
        mask = pe.Node(afni.MaskTool(), name='generate_count_mask')
        allineate = pe.Node(afni.Allineate(), name='allineate')
        catmatvec = pe.Node(afni.CatMatvec(), name='concatenate_transforms')
        apply_allineate3 = pe.Node(afni.Allineate(), name='apply_transform3')
        apply_allineate4 = pe.Node(afni.Allineate(), name='apply_transform4')
        tcat3 = pe.Node(
            afni.TCat(), name='concatenate_across_individuals4')
        tstat3 = pe.Node(afni.TStat(), name='compute_average4')

        workflow.adds([mask, allineate, catmatvec, apply_allineate3,
                       apply_allineate4, tcat3, tstat3])

        workflow.connect(tcat2, 'out_file', mask, 'in_file')
        workflow.connect(mask, 'out_file', allineate, 'weight')
        workflow.connect(apply_allineate1, 'out_file',
                         allineate, 'in_file')
        workflow.connect(allineate, 'out_matrix',
                         catmatvec, 'in_file')
        #XXX how can we enter multiple files ? 
        workflow.connect(catmatvec, 'out_file',
                         apply_allineate3, 'in_matrix')
        workflow.connect(resample1, 'out_file',
                         apply_allineate3, 'in_file')
        workflow.connect(apply_allineate3, 'out_file', tcat3, 'in_files')
        workflow.connect(tcat3, 'out_file', tstat3, 'in_file')

    if pipeline_name == 'nonlinear_anat_to_common':
        pass

    current_dir = os.getcwd()
    os.chdir(write_dir)
    _ = workflow.write_graph(graph2use=graph_kind, dotfilename=graph_file)
    os.chdir(current_dir)
