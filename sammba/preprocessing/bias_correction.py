import os
from ..externals.nipype.caching import Memory
from ..externals.nipype.interfaces import afni, ants, fsl
from ..externals.nipype.utils.filemanip import fname_presuffix
from .utils import compute_n4_max_shrink


def _ants_n4(in_file, write_dir=None, caching=False,
             terminal_output='allatonce', verbose=True, environ=None):
    if write_dir is None:
        write_dir = os.path.dirname(in_file)

    if environ is None:
        environ = {'AFNI_DECONFLICT': 'OVERWRITE'}

    if caching:
        memory = Memory(write_dir)
        bias_correct = memory.cache(ants.N4BiasFieldCorrection)
        copy = memory.cache(afni.Copy)
        copy_geom = memory.cache(fsl.CopyGeom)
        bias_correct.interface().set_default_terminal_output(terminal_output)
        copy.interface().set_default_terminal_output(terminal_output)
    else:
        bias_correct = ants.N4BiasFieldCorrection(
            terminal_output=terminal_output).run
        copy = afni.Copy(terminal_output=terminal_output).run
        copy_geom = fsl.CopyGeom(terminal_output=terminal_output).run

    out_bias_correct = bias_correct(
        input_image=in_file,
        shrink_factor=compute_n4_max_shrink(in_file),
        verbose=verbose,
        output_image=fname_presuffix(in_file, suffix='_n4_deoblique',
                                     newpath=write_dir))
    if False:
        out_copy = copy(
            in_file=out_bias_correct.outputs.output_image,
            out_file=fname_presuffix(in_file,
                                     suffix='_n4',
                                     newpath=write_dir),
            environ=environ,
            verb=verbose)
        out_copy_geom = copy_geom(dest_file=out_copy.outputs.out_file,
                                  in_file=in_file)

    out_copy_geom = copy_geom(dest_file=out_bias_correct.outputs.output_image,
                              in_file=in_file)

    return out_copy_geom.outputs.out_file


def _afni_unifize(in_file, write_dir=None, out_file=None, caching=False,
                  terminal_output='allatonce', verbose=True, environ=None,
                  **unifize_kwargs):
    if write_dir is None:
        write_dir = os.path.dirname(in_file)

    if environ is None:
        environ = {'AFNI_DECONFLICT': 'OVERWRITE'}

    if caching:
        memory = Memory(write_dir)
        copy_geom = memory.cache(fsl.CopyGeom)
        unifize = memory.cache(afni.Unifize)
        copy = memory.cache(afni.Copy)
        unifize.interface().set_default_terminal_output(terminal_output)
        copy.interface().set_default_terminal_output(terminal_output)
    else:
        copy_geom = fsl.CopyGeom(terminal_output=terminal_output).run
        unifize = afni.Unifize(terminal_output=terminal_output).run
        copy = afni.Copy(terminal_output=terminal_output).run

    if out_file is None:
        out_file = fname_presuffix(in_file,
                                   suffix='_unifized',
                                   newpath=write_dir)
    out_unifize = unifize(in_file=in_file,
                          out_file=out_file,
                          environ=environ,
                          quiet=not(verbose),
                          **unifize_kwargs)
    if False:
        out_copy = copy(
            in_file=out_unifize.outputs.output_image,
            out_file=fname_presuffix(in_file,
                                     suffix='_n4',
                                     newpath=write_dir),
            environ=environ,
            verb=verbose)
        out_copy_geom = copy_geom(dest_file=out_copy.outputs.out_file,
                                  in_file=in_file)

    out_copy_geom = copy_geom(dest_file=out_unifize.outputs.out_file,
                              in_file=in_file)
    return out_copy_geom.outputs.out_file
