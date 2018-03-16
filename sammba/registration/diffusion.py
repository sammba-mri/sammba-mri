import os
import numpy as np
from nipype.workflows.data import get_flirt_schedule
from ..externals.nipype.caching import Memory
from ..externals.nipype.interfaces import afni
from ..externals.nipype.utils.filemanip import fname_presuffix
from .base import (BaseSession, _delete_orientation, _rigid_body_register,
                   _warp, _per_slice_qwarp)
from nipype.workflows.dmri.fsl.utils import (b0_average,
                                                        hmc_split,
                                                        _checkinitxfm,
                                                        extract_bval,
                                                        rotate_bvecs,
                                                        insert_mat,
                                                        recompose_dwi,
                                                        recompose_xfm)
from nipype.workflows.dmri.fsl.artifacts import _xfm_jacobian
from ..externals.nipype.interfaces import fsl, ants


def _to_float(in_file, out_file=None):
    import os
    import nibabel
    import numpy as np

    if out_file is None:
        fname, ext = os.path.splitext(os.path.basename(in_file))
        if ext == ".gz":
            fname, ext2 = os.path.splitext(fname)
            ext = ext2 + ext
        out_file = os.path.abspath("%s_float%s" % (fname, ext))

    img = nibabel.load(in_file)
    data = img.get_data().astype(np.float)
    float_img = nibabel.Nifti1Image(data, img.affine, img.header)
    float_img.to_filename(out_file)
    return out_file


def _dwi_flirt(reference, in_file, ref_mask, in_bval, in_xfms=None,
               excl_nodiff=False, flirt_param={}):
    """
    Adapted from nipype dMRI workflows to use it as a function.
    """

    bias_correct = ants.N4BiasFieldCorrection().run
    out_bias_correct = bias_correct(input_image=reference,
                                    mask_image=ref_mask,
                                    dimension=3)

    apply_mask = fsl.maths.ApplyMask().run
    out_apply_mask_b0 = apply_mask(
        in_file=out_bias_correct.outputs.output_image, mask_file=ref_mask)

    dilate = fsl.maths.MathsCommand().run
    out_dilate = dilate(in_file=ref_mask,
                        nan2zeros=True,
                        args='-kernel sphere 5 -dilM')

    init_xfms = _checkinitxfm(in_bval, excl_nodiff, in_xfms=in_xfms)

    split = fsl.Split().run
    out_split = split(in_file=in_file, dimension='t')
    preproc_files = []
    matrix_files = []
    for splitted_file, init_xfm in zip(out_split.outputs.out_files, init_xfms):
        apply_mask = fsl.ApplyMask().run
        out_apply_mask_dwi = apply_mask(in_file=in_file,
                                        mask_file=out_dilate.outputs.out_file)

        flirt = fsl.FLIRT(**flirt_param).run
        out_flirt = flirt(in_weight=out_dilate.outputs.out_file,
                          ref_weight=out_dilate.outputs.out_file,
                          reference=out_apply_mask_b0.outputs.out_file,
                          in_file=out_apply_mask_dwi.outputs.out_file,
                          in_matrix_file=init_xfm)
        threshold = fsl.Threshold().run
        matrix_files.append(out_flirt.outputs.out_matrix_file)
        out_threshold = threshold(in_file=out_flirt.outputs.out_file,
                                  thresh=0.0)
        preproc_files.append(out_threshold.outputs.out_file)

    merge = fsl.Merge(dimension='t').run
    out_merge = merge(in_files=preproc_files)
    return out_merge.outputs.merged_file, matrix_files


def _correct_head_motion(dwi_file, bvecs_file, bvals_file, brain_mask_file,
                         reference_frame=0):
    """
    Adapted from nipype dMRI workflows to use it as a function.
    Correct for head motion artifacts in dMRI sequences.
    It takes a series of diffusion weighted images and rigidly co-registers
    them to one reference image. Finally, the `b`-matrix is rotated accordingly
    [Leemans09]_ making use of the rotation matrix obtained by FLIRT.

    Search angles have been limited to 4 degrees, based on results in
    [Yendiki13]_.

    A list of rigid transformation matrices is provided, so that transforms
    can be chained.
    This is useful to correct for artifacts with only one interpolation process
    (as previously discussed `here
    <https://github.com/nipy/nipype/pull/530#issuecomment-14505042>`_),
    and also to compute nuisance regressors as proposed by [Yendiki13]_.

    .. warning:: This workflow rotates the `b`-vectors, so please be advised
      that not all the dicom converters ensure the consistency between the
      resulting nifti orientation and the gradients table (e.g. dcm2nii
      checks it).

    .. admonition:: References

      .. [Leemans09] Leemans A, and Jones DK, `The B-matrix must be rotated
        when correcting for subject motion in DTI data
        <http://dx.doi.org/10.1002/mrm.21890>`_,
        Magn Reson Med. 61(6):1336-49. 2009. doi: 10.1002/mrm.21890.

      .. [Yendiki13] Yendiki A et al., `Spurious group differences due to head
        motion in a diffusion MRI study
        <http://dx.doi.org/10.1016/j.neuroimage.2013.11.027>`_.
        Neuroimage. 21(88C):79-90. 2013. doi: 10.1016/j.neuroimage.2013.11.027

    Parameters
    ----------
    dwi_file : str
         Path to the dwi file

    bvecs_file : str
         Path to the gradients file (b-vectors)

    bvals_file : str
         Path to the b-values file

    brain_mask_file : str
        Path to the weights mask of reference image (a file with data range
        in [0.0, 1.0], indicating the weight of each voxel when computing
        the metric.

    reference_frame : int, optional
        Index of the b0 volume that should be taken as reference.

    Returns
    -------
    output : tuple of str, str and list of str
        Path to the corrected DW image, Path to the rotated gradient vectors
        table and list of paths to the transformation matrices.

    """
    params = dict(dof=6, save_log=True, no_search=True,
                  schedule=get_flirt_schedule('hmc'))

    [out_ref, out_mov, out_bval, volid] = hmc_split(
        in_file=dwi_file, in_bval=bvals_file, ref_num=reference_frame)

    motion_corrected_file, out_matrices = _dwi_flirt(out_ref, out_mov,
                                                     brain_mask_file,
                                                     out_bval, in_xfms=None,
                                                     excl_nodiff=False,
                                                     flirt_param=params)

    all_matrices = insert_mat(inlist=out_matrices, volid=volid)
    rotated_bvecs = rotate_bvecs(in_bvec=bvecs_file, in_matrix=all_matrices)

    return motion_corrected_file, rotated_bvecs, all_matrices


def _correct_eddy_currents(dwi_file, bvals_file, brain_mask_file,
                           initialization_matrices):
    """
    Adapted from nipype dMRI workflows to use it as a function.
    Corrects for artifacts induced by Eddy currents in dMRI sequences.
    It takes a series of diffusion weighted images and linearly co-registers
    them to one reference image (the average of all b0s in the dataset).

    DWIs are also modulated by the determinant of the Jacobian as indicated by
    [Jones10]_ and [Rohde04]_.

    A list of rigid transformation matrices can be provided, sourcing from a
    :func:`.hmc_pipeline` workflow, to initialize registrations in a *motion
    free* framework.

    A list of affine transformation matrices is available as output, so that
    transforms can be chained (discussion
    `here <https://github.com/nipy/nipype/pull/530#issuecomment-14505042>`_).

    Parameters
    ----------
    dwi_file : str
         Path to the dwi file

    bvecs_file : str
         Path to the gradients file (b-vectors)

    bvals_file : str
         Path to the b-values file

    brain_mask_file : str
        Path to the weights mask of reference image (a file with data range
        in [0.0, 1.0], indicating the weight of each voxel when computing
        the metric.

    initialization_matrices : list of str
        List of Paths to matrices to initialize registration (typically
        from head-motion correction).


    Returns
    -------
    output : tuple of str and list of str
        Path to the corrected DW image and list of paths to the transformation
        matrices.

    .. admonition:: References

      .. [Jones10] Jones DK, `The signal intensity must be modulated by the
        determinant of the Jacobian when correcting for eddy currents in
        diffusion MRI
        <http://cds.ismrm.org/protected/10MProceedings/files/1644_129.pdf>`_,
        Proc. ISMRM 18th Annual Meeting, (2010).

      .. [Rohde04] Rohde et al., `Comprehensive Approach for Correction of
        Motion and Distortion in Diffusion-Weighted MRI
        <http://stbb.nichd.nih.gov/pdf/com_app_cor_mri04.pdf>`_, MRM
        51:103-114 (2004).
    """
    params = dict(dof=12, no_search=True, interp='spline', bgvalue=0,
                  schedule=get_flirt_schedule('ecc'))
    reference_file = b0_average(dwi_file, bvals_file)
    extracted_dwis_file = extract_bval(dwi_file, bvals_file, b='diff')

    ecc_corrected_file, out_matrices = _dwi_flirt(reference_file,
                                                  extracted_dwis_file,
                                                  brain_mask_file,
                                                  bvals_file,
                                                  in_xfms=initialization_matrices,
                                                  excl_nodiff=False,
                                                  flirt_param=params)

    recomposed_matrices = recompose_xfm(in_bval=bvals_file,
                                        in_xfms=out_matrices)
 
    split = fsl.Split(dimension='t').run
    out_split = split(ecc_corrected_file)
    multiply = fsl.BinaryMaths(operation='mul').run
    threshold = fsl.Threshold(thresh=0.0).run
    corrected_files = []
    for splitted_file, matrix in zip(out_split.outputs.out_files,
                                     out_matrices):
        for operand_value in _xfm_jacobian(matrix):
            out_multiply = multiply(operand_value=operand_value,
                                    in_file=splitted_file)
        out_threshold = threshold(out_multiply.outputs.out_file)
        corrected_files.append(out_threshold.outputs.out_file)

    ecc_file = recompose_dwi(in_dwi=dwi_file, in_bval=bvals_file,
                             in_corrected=corrected_files)

    return ecc_file, recomposed_matrices


class DWISession(BaseSession):
    """
    Encapsulation for diffusion data, relative to preprocessing.

    Parameters
    ----------
    dwi : str
        Path to the DW image

    anat : str
        Path to anatomical image

    brain_volume : int, optional
        Volume of the brain used for brain extraction.
        Typically 400 for mouse and 1650 for rat.

    output_dir : str, optional
        Path to the output directory. If not specified, current directory is
        used. Final and intermediate images are stored in the subdirectory
        `animal_id` of the given `output_dir`.
    """
    def __init__(self, dwi=None, bvals=None, anat=None,
                 brain_volume=None, output_dir=None):
        self.dwi = dwi
        self.bvals = bvals
        self.anat = anat
        self.brain_volume = brain_volume
        self.output_dir = output_dir

    def _check_inputs(self):
        if not os.path.isfile(self.dwi):
            raise IOError('dwi must be an existing image file,'
                          'you gave {0}'.format(self.dwi))

        if not os.path.isfile(self.bvals):
            raise IOError('bvals must be an existing b-values file,'
                          'you gave {0}'.format(self.bvals))

        if not os.path.isfile(self.anat):
            raise IOError('anat must be an existing image file,'
                          'you gave {0}'.format(self.anat))

    def coregister(self, delete_orientation=False, use_rats_tool=True,
                   max_b=10.,
                   prior_rigid_body_registration=False,
                   caching=False, voxel_size_x=.1, voxel_size_y=.1,
                   verbose=True, **environ_kwargs):
        """
        Coregistration of the animal's anatomical to the average of
        B0 scans from the DW image.
        Be aware that motion and eddy-current correction are not achieved here.
        The DW image is aligned to the anatomical, first with a
        rigid body registration and then on a per-slice basis (only a fine
        correction, this is mostly for correction of EPI distortion).

        Parameters
        ----------
        delete_orientation : bool, optional
            If True, set maximal resolution to .1 for anat image and .2 for DW
            image (while keeping the same proportions) and set origin to zero.
            Useful if orientation meta-data in headers is corrupted.

        use_rats_tool : bool, optional
            If True, brain mask is computed using RATS Mathematical Morphology.
            Otherwise, a histogram-based brain segmentation is used.

        max_b : float, optional
            Value under which b-values are assumed zero.

        prior_rigid_body_registration : bool, optional
            If True, a rigid-body registration of the anat to the diffusion is
            performed prior to the warp. Useful if the images headers have
            missing/wrong information.

        voxel_size_x : float, optional
            Resampling resolution for the x-axis, in mm.

        voxel_size_y : float, optional
            Resampling resolution for the y-axis, in mm.

        caching : bool, optional
            Wether or not to use caching.

        verbose : bool, optional
            If True, all steps are verbose. Note that caching implies some
            verbosity in any case.

        environ_kwargs : extra arguments keywords
            Extra arguments keywords, passed to all interfaces environ
            variable.

        Returns
        -------
        The following attributes are added
            - `coreg_dwi_` : str
                             Path to paths to the coregistered DW image.
            - `coreg_anat_` : str
                              Path to paths to the coregistered anatomical
                              image.
            - `coreg_transform_` : str
                                   Path to the transform from anat to DW image.
        Notes
        -----
        If `use_rats_tool` is turned on, RATS tool is used for brain extraction
        and has to be cited. For more information, see
        `RATS <http://www.iibi.uiowa.edu/content/rats-overview/>`_
        """
        self._check_inputs()
        self._set_output_dir()

        environ = {'AFNI_DECONFLICT': 'OVERWRITE'}
        for (key, value) in environ_kwargs.items():
            environ[key] = value

        if verbose:
            terminal_output = 'allatonce'
        else:
            terminal_output = 'none'

        if caching:
            memory = Memory(self.output_dir)
            tcat_subbrick = memory.cache(afni.TCatSubBrick)
            tstat = memory.cache(afni.TStat)
            unifize = memory.cache(afni.Unifize)
            catmatvec = memory.cache(afni.CatMatvec)
            for step in [tstat, tcat_subbrick, unifize]:
                step.interface().set_default_terminal_output(terminal_output)
            overwrite = False
        else:
            tcat_subbrick = afni.TCatSubBrick(terminal_output=terminal_output).run
            tstat = afni.TStat(terminal_output=terminal_output).run
            unifize = afni.Unifize(terminal_output=terminal_output).run
            catmatvec = afni.CatMatvec().run
            overwrite = True

        output_files = []

        if delete_orientation:
            dwi_filename = _delete_orientation(self.dwi,
                                               write_dir=self.output_dir,
                                               min_zoom=.2,
                                               caching=caching,
                                               verbose=verbose)
            anat_filename = _delete_orientation(self.anat,
                                                write_dir=self.output_dir,
                                                min_zoom=.1,
                                                caching=caching,
                                                verbose=verbose)
        else:
            dwi_filename = self.dwi
            anat_filename = self.anat

        ####################
        # Average B0 volumes
        ####################
        b_values = np.loadtxt(self.bvals)
        b0_frames = np.argwhere(b_values <= max_b).flatten().tolist()
        out_tcat_subbrick = tcat_subbrick(
            in_files=[(dwi_filename, "'{}'".format(b0_frames))],
            out_file=fname_presuffix(dwi_filename, suffix='_b0s',
                                     newpath=self.output_dir),
            environ=environ)
        out_tstat = tstat(
            in_file=out_tcat_subbrick.outputs.out_file,
            out_file=fname_presuffix(dwi_filename, suffix='_b0s_mean',
                                     newpath=self.output_dir))
        b0_filename = out_tstat.outputs.out_file
        output_files.extend([out_tcat_subbrick.outputs.out_file, b0_filename])

        ###########################################
        # Correct anat and DWI for intensity bias #
        ###########################################
        # Correct the B0 average for intensities bias
        out_bias_correct = unifize(
            in_file=b0_filename,
            out_file=fname_presuffix(b0_filename, suffix='_unifized'),
            environ=environ)
        unbiased_b0_filename = out_bias_correct.outputs.out_file

        # Bias correct the antomical image
        out_unifize = unifize(
            in_file=anat_filename,
            out_file=fname_presuffix(anat_filename, suffix='_unifized',
                                     newpath=self.output_dir),
            environ=environ)
        unbiased_anat_filename = out_unifize.outputs.out_file

        # Update outputs
        output_files.extend([unbiased_b0_filename,
                             unbiased_anat_filename])

        ###########################################
        # Rigid-body registration anat -> mean B0 #
        ###########################################
        if prior_rigid_body_registration:
            if self.brain_volume is None:
                raise ValueError("'brain_volume' value is needed to perform "
                                 "rigid-body registration")
            allineated_anat_filename, rigid_transform_file = \
                _rigid_body_register(unbiased_anat_filename,
                                     unbiased_b0_filename,
                                     self.brain_volume,
                                     use_rats_tool=use_rats_tool,
                                     caching=caching, verbose=verbose,
                                     terminal_output=terminal_output,
                                     environ=environ)
            output_files.extend([rigid_transform_file,
                                 allineated_anat_filename])
        else:
            allineated_anat_filename = unbiased_anat_filename

        ##########################################
        # Nonlinear registration anat -> mean B0 #
        ##########################################
        registered_anat_oblique_filename, mat_filename, warp_output_files =\
            _warp(allineated_anat_filename, unbiased_b0_filename,
                  caching=caching,
                  terminal_output=terminal_output, overwrite=overwrite,
                  verbose=verbose,
                  environ=environ)

        # Save the anat to mean B0 tranform
        transform_filename = fname_presuffix(registered_anat_oblique_filename,
                                             suffix='_anat_to_b0.aff12.1D',
                                             use_ext=False)
        _ = catmatvec(in_file=[(mat_filename, 'ONELINE')],
                      oneline=True,
                      out_file=transform_filename)
        output_files.extend(warp_output_files)

        #####################################################
        # Per-slice non-linear registration mean B0 -> anat #
        #####################################################
        warped_b0_filename, warp_filenames, warped_dwi_filename =\
            _per_slice_qwarp(unbiased_b0_filename,
                             registered_anat_oblique_filename,
                             voxel_size_x, voxel_size_y,
                             apply_to_file=dwi_filename,
                             write_dir=self.output_dir,
                             overwrite=overwrite, verbose=verbose,
                             caching=caching, terminal_output=terminal_output,
                             environ=environ)
        output_files.append(warped_b0_filename)

        if not caching:
            for out_file in output_files:
                os.remove(out_file)

        # Update the diffusion session
        self.coreg_dwi_ = warped_dwi_filename
        self.coreg_anat_ = registered_anat_oblique_filename
        self.coreg_transform_ = transform_filename

        return self
