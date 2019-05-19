import os
import shutil
import tempfile
import numpy as np
import nibabel
from nilearn import image
from sammba.externals.nipype.caching import Memory
from sammba.externals.nipype.utils.filemanip import fname_presuffix
from sammba.externals.nipype.interfaces import afni


def _get_afni_output_type(in_file):
    """ Extracts AFNI outputtype name from file extension
    """
    fname, ext = os.path.splitext(in_file)
    if '.gz' in ext:
        _, ext2 = os.path.splitext(fname)
        ext = ext2 + ext
    if ext == '.nii':
        output_type = 'NIFTI'
    elif ext == '.nii.gz':
        output_type = 'NIFTI_GZ'
    else:
        raise ValueError('Unknown extension for {}'.format(in_file))

    return output_type


def _reset_affines(in_file, out_file, overwrite=False, axes_to_permute=None,
                   axes_to_flip=None, xyzscale=None, center_mass=None,
                   verbose=1):
    """Sets the qform equal to the sform in the header, with optionally
       rescaling, setting the image center to 0, permuting or/and flipping axes
    """
    if not os.path.isfile(out_file) or overwrite:
        shutil.copy(in_file, out_file)
    else:
        raise ValueError('{} already existed'.format(out_file))

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
                  verbose=True, environ=None):
    if caching_dir is None:
        caching_dir = os.getcwd()
    if caching:
        memory = Memory(caching_dir)

    if verbose:
        terminal_output = 'allatonce'
    else:
        terminal_output = 'none'

    if environ is None:
        if caching:
            environ = {}
        else:
            environ = {'AFNI_DECONFLICT': 'OVERWRITE'}

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
    out_copy_oblique = copy(in_file=reference_filename,
                            out_file=orig_reference_filename,
                            environ=environ)
    orig_reference_filename = out_copy_oblique.outputs.out_file

    to_fix_basename = os.path.basename(to_fix_filename)
    orig_to_fix_filename = fname_presuffix(os.path.join(
        tmp_folder, to_fix_basename), suffix='+orig.BRIK',
        use_ext=False)
    out_copy = copy(in_file=to_fix_filename,
                    out_file=orig_to_fix_filename,
                    environ=environ)
    orig_to_fix_filename = out_copy.outputs.out_file

    out_refit = refit(in_file=orig_to_fix_filename,
                      atrcopy=(orig_reference_filename,
                               'IJK_TO_DICOM_REAL'))

    out_copy = copy(in_file=out_refit.outputs.out_file,
                    out_file=fname_presuffix(to_fix_filename,
                                             suffix='_oblique'),
                    environ=environ)

    if not caching:
        shutil.rmtree(tmp_folder)

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
                  copy_shape=True, in_place=True, allow_resampling=False):
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

    allow_resampling : bool, optional
        If True, resampling is done when the two images have different shapes.
    """
    img_to_copy = nibabel.load(filename_to_copy)
    img_to_change = nibabel.load(filename_to_change)
    header_to_copy = img_to_copy.header
    header_to_change = img_to_change.header
    new_header = header_to_change.copy()
    geometry_keys = ['sform_code', 'qform_code', 'quatern_b',
                     'quatern_c', 'quatern_d', 'qoffset_x',
                     'qoffset_y', 'qoffset_z', 'srow_x', 'srow_y', 'srow_z']
    data_to_change = img_to_change.get_data()
    if copy_shape:
        geometry_keys += ['dim']
        geometry_keys += ['dim_info']
        geometry_keys += ['slice_end']
        target_shape = img_to_copy.get_data().shape
        if data_to_change.shape != target_shape:
            if allow_resampling:
                img_to_change = image.resample_to_img(img_to_change,
                                                      img_to_copy)
                data_to_copy = img_to_change.get_data()
                data_to_change = np.squeeze(
                    data_to_copy[..., :header_to_copy['slice_end'] + 1])
            else:
                raise ValueError('images have different shapes: {0} and {1}. '
                                 'You need to set `allow_resampling` to True. '
                                 ''.format(data_to_change.shape, target_shape))
    for key in geometry_keys:
        new_header[key] = header_to_copy[key]

    new_header['pixdim'][:4] = header_to_copy['pixdim'][:4]
    new_img = nibabel.Nifti1Image(data_to_change, img_to_copy.affine,
                                  header=new_header)

    if not in_place:
        if out_filename is None:
            out_filename = fname_presuffix(filename_to_change,
                                           suffix='copied_geom')
    else:
        out_filename = filename_to_change
    new_img.to_filename(out_filename)
    return out_filename


def _check_same_geometry(img_filename1, img_filename2):
    unchecked_fields = ['descrip', 'pixdim', 'scl_slope', 'scl_inter',
                        'xyzt_units', 'qform_code', 'regular']
    header1 = nibabel.load(img_filename1).header
    header2 = nibabel.load(img_filename2).header
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
    return np.alltrue([np.alltrue(a) for a in equal_values])

