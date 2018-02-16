import os
import numpy as np
from nilearn._utils.compat import _basestring
from ..externals.nipype.caching import Memory
from ..externals.nipype.interfaces import afni
from ..externals.nipype.utils.filemanip import fname_presuffix
from .base import (BaseSession, extract_brain, _rigid_body_register, _warp,
                   _per_slice_qwarp)


class DWISession(BaseSession):
    """
    Encapsulation for diffusion data, relative to preprocessing.

    Parameters
    ----------
    dwi : str
        Path to the DW image

    anat : str
        Path to anatomical image

    animal_id : str
        Animal id

    brain_volume : int
        Volume of the brain used for brain extraction.
        Typically 400 for mouse and 1800 for rat.

    output_dir : str, optional
        Path to the output directory. If not specified, current directory is
        used. Final and intermediate images are stored in the subdirectory
        `animal_id` of the given `output_dir`.
    """
    def __init__(self, dwi=None, bvals=None, anat=None,
                 brain_volume=None,
                 animal_id=None, output_dir=None):
        self.dwi = dwi
        self.bvals = bvals
        self.anat = anat
        self.brain_volume = brain_volume
        self.animal_id = animal_id
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

        if self.brain_volume is None:
            raise ValueError('you must provide the expected brain volume.')

        if not isinstance(self.animal_id, _basestring):
            raise ValueError('animal_id must be a string, you provided '
                             '{0}'.format(self.animal_id))

    def coregister(self, use_rats_tool=True, max_b=10.,
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
            Extra arguments keywords, passed to interfaces environ variable.

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
        dwi_filename = self.dwi
        anat_filename = self.anat

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
            unifize.interface().set_default_terminal_output(terminal_output)
            overwrite = False
        else:
            tcat_subbrick = afni.TCatSubBrick(terminal_output=terminal_output).run
            tstat = afni.TStat(terminal_output=terminal_output).run
            unifize = afni.Unifize(terminal_output=terminal_output).run
            catmatvec = afni.CatMatvec().run
            overwrite = True

        self._check_inputs()
        self._set_output_dir()
        images_dir = os.path.join(os.path.abspath(self.output_dir),
                                  self.animal_id)
        current_dir = os.getcwd()
        os.chdir(images_dir)
        output_files = []

        ####################
        # Average B0 volumes
        ####################
        b_values = np.loadtxt(self.bvals)
        b0_frames = np.argwhere(b_values <= max_b).flatten().tolist()
        out_tcat_subbrick = tcat_subbrick(in_files=[(dwi_filename,
                                                    "'{}'".format(b0_frames))],
                                          outputtype='NIFTI_GZ',
                                          environ=environ)
        out_tstat = tstat(in_file=out_tcat_subbrick.outputs.out_file,
                          outputtype='NIFTI_GZ')
        b0_filename = out_tstat.outputs.out_file

        ##########################################
        # Corret anat and DWI for intensity bias #
        ##########################################
        # Correct the B0 average for intensities bias
        out_bias_correct = unifize(in_file=b0_filename,
                                   outputtype='NIFTI_GZ', environ=environ)
        unbiased_b0_filename = out_bias_correct.outputs.out_file

        # Bias correct the antomical image
        out_unifize = unifize(in_file=anat_filename, outputtype='NIFTI_GZ',
                              environ=environ)
        unbiased_anat_filename = out_unifize.outputs.out_file

        # Update outputs
        output_files.extend([unbiased_b0_filename,
                             unbiased_anat_filename])

        ###########################################
        # Rigid-body registration anat -> mean B0 #
        ###########################################
        if prior_rigid_body_registration:
            anat_brain_filename = extract_brain(
                anat_filename, self.output_dir, self.brain_volume,
                caching=caching, use_rats_tool=use_rats_tool,
                terminal_output=terminal_output, environ=environ)
            b0_brain_filename = extract_brain(
                b0_filename, self.output_dir,
                self.brain_volume, caching=caching, use_rats_tool=use_rats_tool,
                terminal_output=terminal_output, environ=environ)
            allineated_anat_filename, rigid_transform_file = \
                _rigid_body_register(unbiased_anat_filename,
                                     anat_brain_filename,
                                     unbiased_b0_filename,
                                     b0_brain_filename, self.output_dir,
                                     caching=caching,
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
                  self.output_dir, caching=caching,
                  terminal_output=terminal_output, overwrite=overwrite,
                  environ=environ)

        # Concatenate all the anat to mean B0 tranforms
        output_files.extend(warp_output_files)
        transform_filename = fname_presuffix(registered_anat_oblique_filename,
                                             suffix='_anat_to_b0.aff12.1D',
                                             use_ext=False)
        _ = catmatvec(in_file=[(mat_filename, 'ONELINE')],
                      oneline=True,
                      out_file=transform_filename)

        #####################################################
        # Per-slice non-linear registration mean B0 -> anat #
        #####################################################
        warped_b0_filename, warp_filenames, warped_dwi_filename =\
            _per_slice_qwarp(unbiased_b0_filename,
                             registered_anat_oblique_filename,
                             self.output_dir, voxel_size_x, voxel_size_y,
                             apply_to_file=dwi_filename,
                             overwrite=overwrite,
                             caching=caching, terminal_output=terminal_output,
                             environ=environ)
        # Update the outputs
        output_files.append(warped_b0_filename)
        if not caching:
            for out_file in output_files:
                os.remove(out_file)

        # Update the diffusion data
        setattr(self, "coreg_dwi_", warped_dwi_filename)
        setattr(self, "coreg_anat_", registered_anat_oblique_filename)
        setattr(self, "coreg_transform_", transform_filename)
        os.chdir(current_dir)
