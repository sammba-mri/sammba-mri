import os
from nipype.interfaces import afni
import nibabel
import shutil
from nipype.caching import Memory
from nipype.utils.filemanip import fname_presuffix


def correct_affines(in_file, out_file, in_place=False, axes_to_permute=None,
                    axes_to_flip=None, xyzscale=None, cm_file=None):
    """Sets the qform equal to the sform in the header, with optionally
       rescaling, setting the image center to 0, permuting or/and flipping axes
    """
    shutil.copy(in_file, out_file)
    in_file = out_file
    if xyzscale is not None:
        refit = afni.Refit()
        refit.inputs.in_file = in_file
        refit.inputs.xyzscale = .1
        result = refit.run()
        in_file = result.outputs.out_file

    if cm_file is not None:
        center_mass = afni.CenterMass()
        center_mass.inputs.in_file = in_file
        center_mass.inputs.cm_file = cm_file
        center_mass.inputs.set_cm = (0, 0, 0)
        result = center_mass.run()
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
    header.set_qform(sform, int(code), strip_shears=False)
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
