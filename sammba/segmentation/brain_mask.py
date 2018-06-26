import os
import numpy as np
import nibabel
from nilearn.image.resampling import coord_transform
from ..externals.nipype.caching import Memory
from ..externals.nipype.interfaces import afni, fsl
from ..externals.nipype.utils.filemanip import fname_presuffix
from . import interfaces


def _compute_volume(mask_img):
    affine_det = np.abs(np.linalg.det(mask_img.affine[:3, :3]))
    mask_data = mask_img.get_data()
    n_voxels_mask = np.sum(mask_data > 0).astype(float)
    return n_voxels_mask * affine_det


def _compute_mask_measures(mask_file):
    """ Outputs the mask

    Parameters
    ----------
    mask_file : str
        Path to the mask image
    """
    # TODO: symmetry, length and width
    mask_img = nibabel.load(mask_file)
    volume = _compute_volume(mask_img)

    mask_data = mask_img.get_data()
    i, j, k = np.where(mask_data != 0)
    voxels_coords = np.array(coord_transform(i, j, k, mask_img.affine)).T
    x_range, y_range, z_range = voxels_coords.max(axis=0) - \
        voxels_coords.min(axis=0)
    return x_range, y_range, z_range, volume


def brain_extraction_report(head_file, brain_volume, write_dir=None,
                            clipping_fractions=[.2, None], use_rats_tool=True,
                            caching=False,
                            terminal_output='allatonce', digits=2):
    if use_rats_tool:
        compute_brain_mask = compute_morpho_brain_mask
    else:
        compute_brain_mask = compute_histo_brain_mask

    masks_measures = []
    for cl_frac in clipping_fractions:
        if cl_frac is None:
            unifize = False
            unifize_kwargs = {}
        else:
            unifize = True
            unifize_kwargs = {'cl_frac': cl_frac}

        brain_mask_file = compute_brain_mask(head_file, brain_volume,
                                             unifize=unifize,
                                             write_dir=write_dir,
                                             caching=caching,
                                             terminal_output=terminal_output,
                                             **unifize_kwargs)

    masks_measures.append(_compute_mask_measures(brain_mask_file))

    target_names = ['{0:0.2f}'.format(cl_frac) if cl_frac is not None
                    else 'None' for cl_frac in clipping_fractions]
    name_width = max(len(cn) for cn in target_names)
    width = max(name_width, digits)

    headers = ["x extent", "y extend", "z extend", "volume"]
    head_fmt = u'{:>{width}s} ' + u' {:>9}' * len(headers)
    report = head_fmt.format(u'', *headers, width=width)
    report += u'\n\n'

    row_fmt = u'{:>{width}s} ' + u' {:>9.{digits}f}' * 3 + u' {:>9}\n'
    x_extent, y_extent, z_extent, volume = zip*(masks_measures)
    rows = zip(target_names, x_extent, y_extent, z_extent, volume)
    for row in rows:
        report += row_fmt.format(*row, width=width, digits=digits)

    report += u'\n'

    return report

    return report


def compute_morpho_brain_mask(head_file, brain_volume, write_dir=None,
                              unifize=True,
                              caching=False,
                              terminal_output='allatonce', **unifize_kwargs):
    """
    Parameters
    ----------
    brain_volume : int
        Volume of the brain used for brain extraction.
        Typically 400 for mouse and 1800 for rat.

    unifize : bool, optional
        If True, brain mask is computed using RATS Mathematical Morphology.
        Otherwise, a histogram-based brain segmentation is used.

    caching : bool, optional
        Wether or not to use caching.

    unifize_kwargs : dict, optional
        Is passed to sammba.externals.nipype.interfaces.afni.Unifize.

    Returns
    -------
    path to brain extracted image.

    Notes
    -----
    RATS tool is used for brain extraction and has to be cited. For more
    information, see
    `RATS <http://www.iibi.uiowa.edu/content/rats-overview/>`_
    """
    if write_dir is None:
        write_dir = os.path.dirname(head_file)

    if interfaces.Info().version() is None:
        raise ValueError('Can not locate Rats')
    environ = {'AFNI_DECONFLICT': 'OVERWRITE'}

    if caching:
        memory = Memory(write_dir)
        clip_level = memory.cache(afni.ClipLevel)
        compute_mask = memory.cache(interfaces.MathMorphoMask)
        unifize = memory.cache(afni.Unifize)
        for step in [compute_mask, unifize]:
            step.interface().set_default_terminal_output(terminal_output)
    else:
        clip_level = afni.ClipLevel().run
        compute_mask = interfaces.MathMorphoMask(
            terminal_output=terminal_output).run
        unifize = afni.Unifize(terminal_output=terminal_output).run

    if unifize:
        if unifize_kwargs is None:
            unifize_kwargs = {}

        file_to_mask = _afni_bias_correct(
            head_file, write_dir,
            out_file=fname_presuffix(head_file,
                                     suffix='_unifized_for_extraction',
                                     newpath=write_dir),
            caching=caching,
            terminal_output=terminal_output,
            environ=environ, **unifize_kwargs)
    else:
        file_to_mask = head_file

    out_clip_level = clip_level(in_file=file_to_mask)
    out_compute_mask = compute_mask(
        in_file=file_to_mask,
        out_file=fname_presuffix(file_to_mask,
                                 suffix='_rats_brain_mask',
                                 newpath=write_dir),
        volume_threshold=brain_volume,
        intensity_threshold=int(out_clip_level.outputs.clip_val))

    # Remove intermediate output
    if not caching and unifize:
        os.remove(file_to_mask)

    return out_compute_mask.outputs.out_file


def compute_histo_brain_mask(head_file, brain_volume, write_dir=None, unifize=True,
                       caching=False,
                       terminal_output='allatonce',
                       lower_cutoff=.2, upper_cutoff=.85, closing=0,
                       connected=True, dilation_size=(1, 1, 2), opening=5,
                       **unifize_kwargs):
    """
    Parameters
    ----------
    brain_volume : int
        Volume of the brain used for brain extraction.
        Typically 400 for mouse and 1800 for rat.

    caching : bool, optional
        Wether or not to use caching.

    unifize_kwargs : dict, optional
        Is passed to sammba.externals.nipype.interfaces.afni.Unifize.

    Returns
    -------
    path to brain extracted image.

    """
    if write_dir is None:
        write_dir = os.path.dirname(head_file)

    environ = {'AFNI_DECONFLICT': 'OVERWRITE'}
    if caching:
        memory = Memory(write_dir)
        clip_level = memory.cache(afni.ClipLevel)
        compute_mask = memory.cache(interfaces.HistogramMask)
        unifize = memory.cache(afni.Unifize)
        for step in [compute_mask, unifize]:
            step.interface().set_default_terminal_output(terminal_output)
    else:
        clip_level = afni.ClipLevel().run
        compute_mask = interfaces.HistogramMask(
            terminal_output=terminal_output).run
        unifize = afni.Unifize(terminal_output=terminal_output).run

    if unifize:
        if unifize_kwargs is None:
            unifize_kwargs = {}

        file_to_mask = _afni_bias_correct(
            head_file, write_dir,
            out_file=fname_presuffix(head_file,
                                     suffix='_unifized_for_extraction',
                                     newpath=write_dir),
            caching=caching,
            terminal_output=terminal_output,
            environ=environ, **unifize_kwargs)
    else:
        file_to_mask = head_file

    out_clip_level = clip_level(in_file=file_to_mask)
    out_compute_mask = compute_mask(
        in_file=file_to_mask,
        out_file=fname_presuffix(file_to_mask,
                                 suffix='_histo_brain_mask',
                                 newpath=write_dir),
        volume_threshold=brain_volume,
        intensity_threshold=int(out_clip_level.outputs.clip_val),
        lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff, closing=closing,
        connected=connected, dilation_size=dilation_size, opening=opening)

    # Remove intermediate output
    if not caching and unifize:
        os.remove(file_to_mask)

    return out_compute_mask.outputs.out_file


def _apply_mask(head_file, mask_file, write_dir=None,
                caching=False, terminal_output='allatonce'):
    """
    Parameters
    ----------
    head_file : str
        Path to the image to mask.

    mask_file : str
        Path to the image mask to apply.

    write_dir : str or None, optional
        Path to the folder to save the output file to. If None, the folder
        of the head file is used.

    caching : bool, optional
        Wether or not to use caching.

    terminal_output : one of {'stream', 'allatonce', 'file', 'none'}
        Control terminal output :
            'stream' : displays to terminal immediately,
            'allatonce' : waits till command is finished to display output,
            'file' : writes output to file
            'none' : output is ignored

    Returns
    -------
    path to brain extracted image.
    """
    if write_dir is None:
        write_dir = os.path.dirname(head_file)

    if caching:
        memory = Memory(write_dir)
        apply_mask = memory.cache(fsl.ApplyMask)
        apply_mask.interface().set_default_terminal_output(terminal_output)
    else:
        apply_mask = fsl.ApplyMask(terminal_output=terminal_output).run

    # Check mask is binary
    mask_img = nibabel.load(mask_file)
    mask = mask_img.get_data()
    values = np.unique(mask)
    if len(values) == 2:
        # If there are 2 different values, one of them must be 0 (background)
        if not 0 in values:
            raise ValueError('Background of the mask must be represented with'
                             '0. Given mask contains: %s.' % values)
    elif len(values) != 2:
        # If there are more than 2 values, the mask is invalid
        raise ValueError('Given mask is not made of 2 values: %s'
                         '. Cannot interpret as true or false' % values)
    try:
        np.testing.assert_array_equal(nibabel.load(mask_file).affine,
                                      nibabel.load(head_file).affine)
    except AssertionError:
        raise ValueError('Given mask {0} and file {1} do not have the same '
                         'affine'.format(mask_file, head_file))

    out_apply_mask = apply_mask(in_file=head_file,
                                mask_file=mask_file,
                                out_file=fname_presuffix(head_file,
                                                         suffix='_masked',
                                                         newpath=write_dir))
    return out_apply_mask.outputs.out_file
