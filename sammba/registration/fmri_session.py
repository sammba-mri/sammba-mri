import os
from nilearn._utils.compat import _basestring
from ..externals.nipype.caching import Memory
from ..externals.nipype.interfaces import afni, fsl
from ..externals.nipype.utils.filemanip import fname_presuffix
from .struct import anats_to_template
from .multi_modality import (_mask_outside_brain, _rigid_body_register, _warp,
                             _per_slice_qwarp)


def _realign(func_filename, write_dir, caching=False,
             terminal_output='allatonce', environ={}):
    if caching:
        memory = Memory(write_dir)
        clip_level = memory.cache(afni.ClipLevel)
        threshold = memory.cache(fsl.Threshold)
        volreg = memory.cache(afni.Volreg)
        allineate = memory.cache(afni.Allineate)
        copy_geom = memory.cache(fsl.CopyGeom)
        tstat = memory.cache(afni.TStat)
        for step in [volreg, allineate, copy_geom,
                     tstat, threshold]:
            step.interface().set_default_terminal_output(terminal_output)
    else:
        clip_level = afni.ClipLevel().run
        threshold = fsl.Threshold(terminal_output=terminal_output).run
        volreg = afni.Volreg(terminal_output=terminal_output).run
        allineate = afni.Allineate(terminal_output=terminal_output).run
        copy_geom = fsl.CopyGeom(terminal_output=terminal_output).run
        tstat = afni.TStat(terminal_output=terminal_output).run

    out_clip_level = clip_level(in_file=func_filename)
    out_threshold = threshold(in_file=func_filename,
                              thresh=out_clip_level.outputs.clip_val)
    thresholded_filename = out_threshold.outputs.out_file

    out_volreg = volreg(  # XXX dfile not saved
        in_file=thresholded_filename,
        outputtype='NIFTI_GZ',
        environ=environ,
        oned_file=fname_presuffix(thresholded_filename,
                                  suffix='_realign.1Dfile.1D', use_ext=False),
        oned_matrix_save=fname_presuffix(thresholded_filename,
                                         suffix='_realign.aff12.1D',
                                         use_ext=False))

    # Apply the registration to the whole head
    out_allineate = allineate(in_file=func_filename,
                              master=func_filename,
                              in_matrix=out_volreg.outputs.oned_matrix_save,
                              out_file=fname_presuffix(func_filename,
                                                       suffix='_realign'),
                              environ=environ)

    # 3dAllineate removes the obliquity. This is not a good way to readd it as
    # removes motion correction info in the header if it were an AFNI file...as
    # it happens it's NIfTI which does not store that so irrelevant!
    out_copy_geom = copy_geom(dest_file=out_allineate.outputs.out_file,
                              in_file=out_volreg.outputs.out_file)

    allineated_filename = out_copy_geom.outputs.out_file

    # Create a (hopefully) nice mean image for use in the registration
    out_tstat = tstat(in_file=allineated_filename, args='-mean',
                      outputtype='NIFTI_GZ', environ=environ)

    # Remove intermediate outputs
    if not caching:
        for output_file in [thresholded_filename,
                            out_volreg.outputs.oned_matrix_save,
                            out_volreg.outputs.out_file,
                            out_volreg.outputs.md1d_file]:
            os.remove(output_file)
    return (allineated_filename, out_tstat.outputs.out_file,
            out_volreg.outputs.oned_file)


def _func_to_template(func_coreg_filename, template_filename, write_dir,
                      func_to_anat_oned_filename,
                      anat_to_template_oned_filename,
                      anat_to_template_warp_filename,
                      caching=False, verbose=True):
    """ Applies successive transforms to coregistered functional to put it in
    template space.

    Parameters
    ----------
    coreg_func_filename : str
        Path to functional volume, coregistered to a common space with the
        anatomical volume.

    template_filename : str
        Template to register the functional to.

    func_to_anat_oned_filename : str
        Path to the affine 1D transform from functional to coregistration
        space.

    anat_to_template_oned_filename : str
        Path to the affine 1D transform from anatomical to template space.

    anat_to_template_warp_filename : str
        Path to the warp transform from anatomical to template space.

    caching : bool, optional
        Wether or not to use caching.

    verbose : bool, optional
        If True, all steps are verbose. Note that caching implies some
        verbosity in any case.
    """
    if verbose:
        terminal_output = 'allatonce'
    else:
        terminal_output = 'none'

    if caching:
        memory = Memory(write_dir)
        warp_apply = memory.cache(afni.NwarpApply)
        warp_apply.interface().set_default_terminal_output(terminal_output)
    else:
        warp_apply = afni.NwarpApply(terminal_output=terminal_output).run

    current_dir = os.getcwd()
    os.chdir(write_dir)
    normalized_filename = fname_presuffix(func_coreg_filename,
                                          suffix='_normalized')
    warp = "'{0} {1} {2}'".format(anat_to_template_warp_filename,
                                  anat_to_template_oned_filename,
                                  func_to_anat_oned_filename)
    _ = warp_apply(in_file=func_coreg_filename,
                   master=template_filename,
                   warp=warp,
                   out_file=normalized_filename)
    os.chdir(current_dir)
    return normalized_filename


class FMRISession(object):
    """
    Encapsulation for fMRI data, relative to preprocessing.

    Parameters
    ----------
    func : str
        Path to the functional 4D image

    anat : str
        Path to anatomical image

    animal_id : str
        Animal id
    """

    def __init__(self, func=None, anat=None, animal_id=None):
        self.func = func
        self.anat = anat
        self.animal_id = animal_id

    def _set_items(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def _check_inputs(self):
        if self.func is not None:
            if not os.path.isfile(self.func):
                raise IOError('func must be an existing image file,'
                              'you gave {0}'.format(self.func))

        if self.anat is not None:
            if not os.path.isfile(self.anat):
                raise IOError('anat must be an existing image file,'
                              'you gave {0}'.format(self.anat))

        if not isinstance(self.animal_id, _basestring):
            raise ValueError('animal_id must be a string, you provided '
                             '{0}'.format(self.animal_id))

    def _set_output_dir_(self, output_dir):
        setattr(self, 'output_dir_', output_dir)
        if hasattr(self, 'output_dir_'):
            if not os.path.isdir(self.output_dir_):
                os.makedirs(self.output_dir_)

    def coregister(self, t_r, brain_volume, write_dir, slice_timing=True,
                   prior_rigid_body_registration=False,
                   caching=False, voxel_size_x=.1, voxel_size_y=.1,
                   verbose=True, **environ_kwargs):
        """
        Coregistration of the subject's functional and anatomical images.
        The functional volume is aligned to the anatomical, first with a
        rigid body registration and then on a per-slice basis (only a fine
        correction, this is mostly for correction of EPI distortion).

        Parameters
        ----------
        t_r : float
            Repetition time for the EPI, in seconds.

        brain_volume : int
            Volumes of the brain as passed to Rats_MM brain extraction tool.
            Typically 400 for mouse and 1800 for rat.

        write_dir : str
            Directory to save the output and temporary images.

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
            - `output_dir_` : str
                              Path to the output directory.
            - `coreg_func_` : str
                              Path to paths to the coregistered functional
                              image.
            - `coreg_anat_` : str
                              Path to paths to the coregistered functional
                              image.
            - `coreg_transform_` : str
                                   Path to the transform from anat to func.
        """
        func_filename = self.func
        anat_filename = self.anat

        environ = {'AFNI_DECONFLICT': 'OVERWRITE'}
        for (key, value) in environ_kwargs.items():
            environ[key] = value

        if verbose:
            terminal_output = 'allatonce'
        else:
            terminal_output = 'none'

        if caching:
            memory = Memory(write_dir)
            tshift = memory.cache(afni.TShift)
            unifize = memory.cache(afni.Unifize)
            catmatvec = memory.cache(afni.CatMatvec)
            for step in [tshift, unifize]:
                step.interface().set_default_terminal_output(terminal_output)
        else:
            tshift = afni.TShift(terminal_output=terminal_output).run
            unifize = afni.Unifize(terminal_output=terminal_output).run
            catmatvec = afni.CatMatvec().run

        self._check_inputs()
        output_dir = os.path.join(os.path.abspath(write_dir),
                                  self.animal_id)
        self._set_output_dir_(output_dir)
        current_dir = os.getcwd()
        os.chdir(output_dir)
        output_files = []

        #######################################
        # Correct functional for slice timing #
        #######################################
        if slice_timing:
            out_tshift = tshift(in_file=func_filename,
                                outputtype='NIFTI_GZ',
                                tpattern='altplus',
                                tr=str(t_r),
                                environ=environ)
            func_filename = out_tshift.outputs.out_file
            output_files.append(func_filename)

        ################################################
        # Register functional volumes to the first one #
        ################################################
        allineated_filename, mean_aligned_filename, _ = \
            _realign(func_filename, self.output_dir_, caching=caching,
                     terminal_output=terminal_output, environ=environ)

        ###########################################
        # Corret anat and func for intensity bias #
        ###########################################
        # Correct the functional average for intensities bias
        out_bias_correct = unifize(in_file=mean_aligned_filename,
                                   outputtype='NIFTI_GZ', environ=environ)
        unbiased_mean_func_filename = out_bias_correct.outputs.out_file

        # Bias correct the antomical image
        out_unifize = unifize(in_file=anat_filename, outputtype='NIFTI_GZ',
                              environ=environ)
        unbiased_anat_filename = out_unifize.outputs.out_file

        # Update outputs
        output_files.extend([unbiased_mean_func_filename,
                             unbiased_anat_filename])

        #############################################
        # Rigid-body registration anat -> mean func #
        #############################################
        if prior_rigid_body_registration:
            anat_brain_filename = _mask_outside_brain(
                unbiased_anat_filename, write_dir, brain_volume, caching=False,
                terminal_output=terminal_output, environ=environ)
            func_brain_filename = _mask_outside_brain(
                unbiased_mean_func_filename, write_dir, brain_volume,
                caching=False, terminal_output=terminal_output,
                environ=environ)
            allineated_anat_filename, rigid_transform_file = \
                _rigid_body_register(unbiased_anat_filename,
                                     anat_brain_filename,
                                     unbiased_mean_func_filename,
                                     func_brain_filename, write_dir,
                                     brain_volume, caching=False,
                                     terminal_output=terminal_output,
                                     environ=environ)
            output_files.append(allineated_anat_filename)
        else:
            allineated_anat_filename = unbiased_anat_filename

        ############################################
        # Nonlinear registration anat -> mean func #
        ############################################
        registered_anat_oblique_filename, mat_filename, warp_output_files =\
            _warp(allineated_anat_filename, unbiased_mean_func_filename,
                  self.output_dir_, caching=caching,
                  terminal_output=terminal_output,
                  environ=environ)
        # Concatenate all the anat to func tranforms
        output_files.extend(warp_output_files)
        transform_filename = fname_presuffix(registered_anat_oblique_filename,
                                             suffix='_anat_to_func.aff12.1D',
                                             use_ext=False)
        if prior_rigid_body_registration:
            _ = catmatvec(in_file=[(mat_filename, 'ONELINE'),
                                   (rigid_transform_file, 'ONELINE')],
                          oneline=True,
                          out_file=transform_filename)
        else:
            _ = catmatvec(in_file=[(mat_filename, 'ONELINE')],
                          oneline=True,
                          out_file=transform_filename)

        ##################################################
        # Per-slice non-linear registration func -> anat #
        ##################################################
        warped_mean_func_filename, warp_filenames, warped_func_filename =\
            _per_slice_qwarp(unbiased_mean_func_filename,
                             registered_anat_oblique_filename,
                             self.output_dir_, voxel_size_x, voxel_size_y,
                             apply_to_filename=allineated_filename,
                             caching=caching, terminal_output=terminal_output,
                             environ=environ)
        # Update the outputs
        output_files.append(warped_mean_func_filename)
        if not caching:
            for out_file in output_files:
                os.remove(out_file)

        # Update the fmri data
        setattr(self, "coreg_func_", warped_func_filename)
        setattr(self, "coreg_anat_", registered_anat_oblique_filename)
        setattr(self, "coreg_transform_", transform_filename)
        os.chdir(current_dir)

    def register_to_template(self, t_r, brain_volume, head_template_filename,
                             write_dir, brain_template_filename=None,
                             dilated_head_mask_filename=None,
                             prior_rigid_body_registration=False,
                             slice_timing=True,
                             maxlev=2,
                             caching=False, verbose=True):
        """ Registration of subject's functional and anatomical images to
        a given template.

        Parameters
        ----------
        t_r : float
            Repetition time for the EPI, in seconds.

        head_template_filename : str
            Template to register the functional to.

        brain_volume : int
            Volumes of the brain as passed to Rats_MM brain extraction tool.
            Typically 400 for mouse and 1800 for rat.

        write_dir : str
            Path to the affine 1D transform from anatomical to template space.

        brain_template_filename : str, optional
            Path to a brain template, passed to
            sammba.registration.anats_to_template

        dilated_head_mask_filename : str, optional
            Path to a dilated head mask, passed to
            sammba.registration.anats_to_template

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
        The following attributes are added
            - `template_` : str
                           Path to the given registration template.
            - `output_dir_` : str
                              Path to the output directory for each animal.
            - `coreg_func_` : str
                              Path to the coregistered functional image.
            - `coreg_anat_` : str
                              Path to the coregistered anatomical image.
            - `coreg_transform_` : str
                                   Path to the transform from anat to func.
            - `registered_func_` : str
                                   Path to the funct registered to template.
            - `registered_anat_` : str
                                   Path to the anat registered to template.

        See also
        --------
        sammba.registration.coregister_fmri_session,
        sammba.registration.anats_to_template
        """
        self._check_inputs()
        animal_output_dir = os.path.join(os.path.abspath(write_dir),
                                         self.animal_id)
        self._set_output_dir_(animal_output_dir)
        # XXX do a function for creating new attributes ?
        setattr(self, "template_", head_template_filename)

        if not hasattr(self, 'coreg_func_') or not hasattr(self,
                                                           'coreg_transform_'):
            self.coregister(
                t_r, brain_volume, write_dir,
                prior_rigid_body_registration=prior_rigid_body_registration,
                slice_timing=slice_timing,
                caching=caching, verbose=verbose)

        anats_registration = anats_to_template(
            [self.anat],
            head_template_filename,
            self.output_dir_,
            brain_volume,
            brain_template_filename=brain_template_filename,
            dilated_head_mask_filename=dilated_head_mask_filename,
            maxlev=maxlev,
            caching=caching, verbose=verbose)
        setattr(self, "registered_anat_", anats_registration.registered[0])

        normalized_func_filename = _func_to_template(
            self.coreg_func_,
            head_template_filename,
            self.output_dir_,
            self.coreg_transform_,
            anats_registration.pre_transforms[0],
            anats_registration.transforms[0],
            caching=caching, verbose=verbose)

        setattr(self, "registered_func_", normalized_func_filename)
