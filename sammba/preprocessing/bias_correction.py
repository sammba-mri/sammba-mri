import os
import nibabel
from ..externals.nipype.caching import Memory
from ..externals.nipype.interfaces import afni, ants, fsl
from ..externals.nipype.utils.filemanip import fname_presuffix

def _compute_n4_max_shrink(in_file):
    """ Computes the maximal allowed shrink factor for ANTS
    N4BiasFieldCorrection.

    Note
    -----
    To lessen computation time, N4BiasFieldCorrection can resample the input
    image. The shrink factor, specified as a single integer, describes this
    resampling. The default shrink factor is 4 which is only applied to the
    first two or three dimensions assumed spatial. The spacing for each
    dimension is computed as
    dimension - shrink - dimension % shrink
    N4BiasFieldCorrection raises and error when one obtained spacing is not
    positive.
    """
    img = nibabel.load(in_file)
    spatial_shapes = img.shape[:3]
    return min(4, min(spatial_shapes) / 2)


def ants_n4(in_file, write_dir=None, caching=False,
            terminal_output='allatonce', verbose=True, environ=None,
            copy_geometry=True):
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

    unbiased_file = fname_presuffix(in_file, suffix='_n4',
                                    newpath=write_dir)
    if copy_geometry:
        output_image = fname_presuffix(in_file, suffix='_n4_rough_geom',
                                       newpath=write_dir)
    else:                                     
        output_image = unbiased_file

    out_bias_correct = bias_correct(
        input_image=in_file,
        shrink_factor=_compute_n4_max_shrink(in_file),
        verbose=verbose,
        output_image=output_image)

    if copy_geometry:
        out_copy = copy(
            in_file=out_bias_correct.outputs.output_image,
            out_file=unbiased_file,
            environ=environ)
        out_copy_geom = copy_geom(dest_file=out_copy.outputs.out_file,
                                  in_file=in_file)
    return unbiased_file


def afni_unifize(in_file, write_dir=None, out_file=None, caching=False,
                 terminal_output='allatonce', verbose=True, environ=None,
                 copy_geometry=False,
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
    if copy_geometry:
        unifized_file = fname_presuffix(in_file,
                                        suffix='_unifized_rough_geom',
                                        newpath=write_dir)
    else:
        unifized_file = out_file
        
    out_unifize = unifize(in_file=in_file,
                          out_file=unifized_file,
                          environ=environ,
                          quiet=not(verbose),
                          **unifize_kwargs)

    if copy_geometry:
        out_copy = copy(
            in_file=out_unifize.outputs.out_file,
            out_file=out_file,
            environ=environ)
        out_copy_geom = copy_geom(dest_file=out_copy.outputs.out_file,
                                  in_file=in_file)
    return out_file
