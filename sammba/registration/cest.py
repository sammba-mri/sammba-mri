import os
from nilearn._utils.compat import _basestring
from ..externals.nipype.caching import Memory
from ..externals.nipype.interfaces import afni
from ..externals.nipype.utils.filemanip import fname_presuffix
from .struct import anats_to_template
from .multimodal import (extract_brain, _rigid_body_register, _warp,
                         _transform_to_template)


class CESTSession(object):
    """
    Encapsulation for CEST data, relative to preprocessing.

    Parameters
    ----------
    cest : str
        Path to the CEST image

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

    def __init__(self, cest=None, anat=None, brain_volume=None,
                 animal_id=None, output_dir=None):
        self.cest = cest
        self.anat = anat
        self.brain_volume = brain_volume
        self.animal_id = animal_id
        self.output_dir = output_dir

    def _set_items(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def _check_inputs(self):
        if not os.path.isfile(self.cest):
            raise IOError('func must be an existing image file,'
                          'you gave {0}'.format(self.cest))

        if not os.path.isfile(self.anat):
            raise IOError('anat must be an existing image file,'
                          'you gave {0}'.format(self.anat))

        if self.brain_volume is None:
            raise ValueError('you must provide the expected brain volume.')

        if not isinstance(self.animal_id, _basestring):
            raise ValueError('animal_id must be a string, you provided '
                             '{0}'.format(self.animal_id))

    def _set_output_dir(self):
        if self.output_dir is None:
            self.output_dir = os.getcwd()

        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)

    def coregister(self, use_rats_tool=True,
                   prior_rigid_body_registration=False,
                   caching=False, voxel_size_x=.1, voxel_size_y=.1,
                   verbose=True, **environ_kwargs):
        """
        Coregistration of the animal's CEST and anatomical images.
        The CEST volume is aligned to the anatomical, first with a
        rigid body registration and then a non linear warp.

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
            - `coreg_cest_` : str
                              Path to paths to the coregistered CEST image.
            - `coreg_anat_` : str
                              Path to paths to the coregistered functional
                              image.
            - `coreg_transform_` : str
                                   Path to the transform from anat to CEST.
        Notes
        -----
        If `use_rats_tool` is turned on, RATS tool is used for brain extraction
        and has to be cited. For more information, see
        `RATS <http://www.iibi.uiowa.edu/content/rats-overview/>`_
        """
        cest_filename = self.cest
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
            unifize = memory.cache(afni.Unifize)
            catmatvec = memory.cache(afni.CatMatvec)
            unifize.interface().set_default_terminal_output(terminal_output)
            overwrite = False
        else:
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

        ###########################################
        # Corret anat and func for intensity bias #
        ###########################################
        # Correct the functional average for intensities bias
        out_bias_correct = unifize(in_file=cest_filename,
                                   outputtype='NIFTI_GZ', environ=environ)
        unbiased_cest_filename = out_bias_correct.outputs.out_file

        # Bias correct the antomical image
        out_unifize = unifize(in_file=anat_filename, outputtype='NIFTI_GZ',
                              environ=environ)
        unbiased_anat_filename = out_unifize.outputs.out_file

        # Update outputs
        output_files.extend([unbiased_cest_filename,
                             unbiased_anat_filename])

        ########################################
        # Rigid-body registration anat -> cest #
        ########################################
        if prior_rigid_body_registration:
            anat_brain_filename = extract_brain(
                unbiased_anat_filename, self.output_dir, self.brain_volume,
                caching=False, use_rats_tool=use_rats_tool,
                terminal_output=terminal_output, environ=environ)
            func_brain_filename = extract_brain(
                unbiased_cest_filename, self.output_dir,
                self.brain_volume, caching=False, use_rats_tool=use_rats_tool,
                terminal_output=terminal_output, environ=environ)
            allineated_anat_filename, rigid_transform_file = \
                _rigid_body_register(unbiased_anat_filename,
                                     anat_brain_filename,
                                     unbiased_cest_filename,
                                     func_brain_filename, self.output_dir,
                                     self.brain_volume, caching=False,
                                     terminal_output=terminal_output,
                                     environ=environ)
            output_files.extend([anat_brain_filename,
                                 func_brain_filename,
                                 rigid_transform_file,
                                 allineated_anat_filename])
        else:
            allineated_anat_filename = unbiased_anat_filename

        ############################################
        # Nonlinear registration anat -> mean func #
        ############################################
        registered_anat_oblique_filename, mat_filename, warp_output_files =\
            _warp(allineated_anat_filename, unbiased_cest_filename,
                  self.output_dir, caching=caching,
                  terminal_output=terminal_output, overwrite=overwrite,
                  environ=environ)

        # Concatenate all the anat to func tranforms
        output_files.extend(warp_output_files)
        transform_filename = fname_presuffix(registered_anat_oblique_filename,
                                             suffix='_anat_to_cest.aff12.1D',
                                             use_ext=False)
        _ = catmatvec(in_file=[(mat_filename, 'ONELINE')],
                      oneline=True,
                      out_file=transform_filename)

        if not caching:
            for out_file in output_files:
                os.remove(out_file)

        # Update the fmri data
        setattr(self, "coreg_anat_", registered_anat_oblique_filename)
        setattr(self, "coreg_transform_", transform_filename)
        os.chdir(current_dir)

    def register_to_template(self, head_template_filename,
                             brain_template_filename=None,
                             dilated_head_mask_filename=None,
                             prior_rigid_body_registration=False,
                             slice_timing=True,
                             cest_voxel_size=None,
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

        cest_voxel_size : 3-tuple of floats, optional
            Voxel size of the registered CEST, in mm.

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
            - `registered_cest_` : str
                                   Path to the CEST registered to template.
            - `registered_anat_` : str
                                   Path to the anat registered to template.

        See also
        --------
        sammba.registration.anats_to_template
        """
        self._check_inputs()
        images_dir = os.path.join(os.path.abspath(self.output_dir),
                                  self.animal_id)

        if not hasattr(self, 'coreg_transform_'):
            raise ValueError('Anatomical image has not been registered '
                             'to CEAT. Please use `coreg` function first')

        # XXX do a function for creating new attributes ?
        setattr(self, "template_", head_template_filename)
        anats_registration = anats_to_template(
            [self.anat],
            head_template_filename,
            images_dir,
            self.brain_volume,
            brain_template_filename=brain_template_filename,
            dilated_head_mask_filename=dilated_head_mask_filename,
            maxlev=maxlev,
            caching=caching, verbose=verbose)
        setattr(self, "registered_anat_", anats_registration.registered[0])

        normalized_cest_filename = _transform_to_template(
            self.cest,
            head_template_filename,
            images_dir,
            [self.coreg_transform_, anats_registration.pre_transforms[0],
             anats_registration.transforms[0]],
            voxel_size=cest_voxel_size, caching=caching, verbose=verbose)

        setattr(self, "registered_cest_", normalized_cest_filename)
