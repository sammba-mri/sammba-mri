import os
from sammba.externals.nipype.caching import Memory
from sammba.externals.nipype.interfaces import afni, ants
from sammba.externals.nipype.utils.filemanip import fname_presuffix
from .struct import anats_to_template
from .utils import compute_n4_max_shrink
from .base import (BaseSession, extract_brain, _rigid_body_register, _warp,
                   _per_slice_qwarp, _transform_to_template)


class PerfSession(BaseSession):
    """
    Encapsulation for perfusion data, relative to preprocessing.

    Parameters
    ----------
    perf : str
        Path to the functional 4D image

    anat : str
        Path to anatomical image

    brain_volume : int, optional
        Volume of the brain used for brain extraction.
        Typically 400 for mouse and 1650 for rat.

    output_dir : str, optional
        Path to the output directory. If not specified, current directory is
        used.
    """

    def __init__(self, perf=None, m0=None, anat=None, brain_volume=None,
                 t_r=None, output_dir=None):
        self.perf = perf
        self.m0 = m0
        self.anat = anat
        self.brain_volume = brain_volume
        self.output_dir = output_dir

    def _check_inputs(self):
        if not os.path.isfile(self.perf):
            raise IOError('perf must be an existing image file,'
                          'you gave {0}'.format(self.perf))

        if not os.path.isfile(self.m0):
            raise IOError('m0 must be an existing image file,'
                          'you gave {0}'.format(self.m0))

        if not os.path.isfile(self.anat):
            raise IOError('anat must be an existing image file,'
                          'you gave {0}'.format(self.anat))

    def coregister(self, use_rats_tool=True,
                   prior_rigid_body_registration=False,
                   voxel_size_x=.1, voxel_size_y=.1, caching=False,
                   verbose=True, **environ_kwargs):
        """
        Coregistration of the subject's functional and anatomical images.
        The functional volume is aligned to the anatomical, first with a
        rigid body registration and then on a per-slice basis (only a fine
        correction, this is mostly for correction of EPI distortion).

        Parameters
        ----------
        use_rats_tool : bool, optional
            If True, brain mask is computed using RATS Mathematical Morphology.
            Otherwise, a histogram-based brain segmentation is used.

        prior_rigid_body_registration : bool, optional
            If True, a rigid-body registration of the anat to the func is
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
            - `coreg_perf_` : str
                              Path to paths to the coregistered perfusion
                              image.
            - `coreg_anat_` : str
                              Path to paths to the coregistered functional
                              image.
            - `coreg_transform_` : str
                                   Path to the transform from anat to func.
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
            unifize = memory.cache(afni.Unifize)
            bias_correct = memory.cache(ants.N4BiasFieldCorrection)
            catmatvec = memory.cache(afni.CatMatvec)
            unifize.interface().set_default_terminal_output(terminal_output)
            overwrite = False
        else:
            unifize = afni.Unifize(terminal_output=terminal_output).run
            bias_correct = ants.N4BiasFieldCorrection(
                terminal_output=terminal_output).run
            catmatvec = afni.CatMatvec().run
            overwrite = True

        output_files = []

        ###########################################
        # Corret anat and M0 for intensity bias #
        ###########################################
        # Bias correct the antomical image
        out_unifize = unifize(in_file=self.anat,
                              out_file=fname_presuffix(self.anat,
                                                       suffix='_unifized',
                                                       newpath=self.output_dir),
                              environ=environ)
        unbiased_anat_filename = out_unifize.outputs.out_file

        # Correct the functional average for intensities bias
        out_bias_correct = bias_correct(
            input_image=self.m0,
            shrink_factor=compute_n4_max_shrink(self.m0),
            output_image=fname_presuffix(self.m0, suffix='_unbiased',
                                         newpath=self.output_dir))

        unbiased_m0_filename = out_bias_correct.outputs.output_image

        # Update outputs
        output_files.extend([unbiased_m0_filename,
                             unbiased_anat_filename])

        #############################################
        # Rigid-body registration anat -> mean func #
        #############################################
        if prior_rigid_body_registration:
            if self.brain_volume is None:
                raise ValueError("'brain_volume' value is needed to perform "
                                 "rigid-body registration")
            allineated_anat_filename, rigid_transform_file = \
                _rigid_body_register(unbiased_anat_filename,
                                     unbiased_m0_filename,
                                     self.output_dir, self.brain_volume,
                                     use_rats_tool=use_rats_tool,
                                     caching=caching,
                                     terminal_output=terminal_output,
                                     environ=environ)
            output_files.extend([rigid_transform_file,
                                 allineated_anat_filename])
        else:
            allineated_anat_filename = unbiased_anat_filename

        ############################################
        # Nonlinear registration anat -> mean func #
        ############################################
        registered_anat_oblique_filename, mat_filename =\
            _warp(allineated_anat_filename, unbiased_m0_filename,
                  caching=caching, verbose=verbose,
                  terminal_output=terminal_output, overwrite=overwrite,
                  environ=environ)

        # Concatenate all the anat to func tranforms
        output_files.append(mat_filename)
        transform_filename = fname_presuffix(registered_anat_oblique_filename,
                                             suffix='_anat_to_func.aff12.1D',
                                             use_ext=False)
        _ = catmatvec(in_file=[(mat_filename, 'ONELINE')],
                      oneline=True,
                      out_file=transform_filename)

        ##################################################
        # Per-slice non-linear registration func -> anat #
        ##################################################
        warped_m0_filename, warp_filenames, warped_perf_filename =\
            _per_slice_qwarp(unbiased_m0_filename,
                             registered_anat_oblique_filename,
                             voxel_size_x, voxel_size_y,
                             apply_to_file=self.perf,
                             verbose=verbose,
                             write_dir=self.output_dir, overwrite=overwrite,
                             caching=caching, terminal_output=terminal_output,
                             environ=environ)
        print(output_files)

        # Remove the intermediate outputs
        if not caching:
            for out_file in output_files:
                os.remove(out_file)

        # Update the fmri data
        setattr(self, "coreg_perf_", warped_perf_filename)
        setattr(self, "coreg_m0_", warped_m0_filename)
        setattr(self, "coreg_anat_", registered_anat_oblique_filename)
        setattr(self, "coreg_transform_", transform_filename)
        setattr(self, "coreg_warps_", warp_filenames)

    def register_to_template(self, head_template_filename,
                             brain_template_filename=None,
                             dilated_head_mask_filename=None,
                             prior_rigid_body_registration=False,
                             slice_timing=True,
                             func_voxel_size=None,
                             maxlev=None,
                             caching=False, verbose=True):
        """ Registration of subject's functional and anatomical images to
        a given template.

        Parameters
        ----------
        head_template_filename : str
            Template to register the functional to.

        brain_template_filename : str, optional
            Path to a brain template, passed to
            sammba.registration.anats_to_template

        dilated_head_mask_filename : str, optional
            Path to a dilated head mask, passed to
            sammba.registration.anats_to_template

        func_voxel_size : 3-tuple of floats, optional
            Voxel size of the registered functional, in mm.

        maxlev : int or None, optional
            Maximal level for the warp when registering anat to template.
            Passed to
            sammba.registration.anats_to_template

        caching : bool, optional
            Wether or not to use caching.

        verbose : bool, optional
            If True, all steps are verbose. Note that caching implies some
            verbosity in any case.

        Returns
        -------
        The following attributes are added/updated
            - `template_` : str
                           Path to the given registration template.
            - `registered_func_` : str
                                   Path to the funct registered to template.
            - `registered_anat_` : str
                                   Path to the anat registered to template.

        See also
        --------
        sammba.registration.anats_to_template
        """
        self._check_inputs()
        if not hasattr(self, 'coreg_m0_') or not hasattr(self,
                                                         'coreg_transform_'):
            raise ValueError('Anatomical and M0 have not been '
                             'coregistered. Please use `coreg` function first')

        # XXX do a function for creating new attributes ?
        setattr(self, "template_", head_template_filename)
        anats_registration = anats_to_template(
            [self.anat],
            head_template_filename,
            self.output_dir,
            self.brain_volume,
            brain_template_filename=brain_template_filename,
            dilated_head_mask_filename=dilated_head_mask_filename,
            maxlev=maxlev,
            caching=caching, verbose=verbose)
        setattr(self, "registered_anat_", anats_registration.registered[0])

        normalized_m0_filename = _transform_to_template(
            self.coreg_m0_,
            head_template_filename,
            self.output_dir,
            [self.coreg_transform_, anats_registration.pre_transforms[0],
             anats_registration.transforms[0]],
            voxel_size=func_voxel_size, caching=caching, verbose=verbose)

        setattr(self, "registered_m0_", normalized_m0_filename)
        normalized_perf_filename = _transform_to_template(
            self.coreg_perf_,
            head_template_filename,
            self.output_dir,
            [self.coreg_transform_, anats_registration.pre_transforms[0],
             anats_registration.transforms[0]],
            voxel_size=func_voxel_size, caching=caching, verbose=verbose)

        setattr(self, "registered_perf_", normalized_perf_filename)
