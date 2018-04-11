import os
from .base import (compute_brain_mask, _bias_correct,
                   _afni_bias_correct, _apply_mask, _transform_to_template,
                   _apply_perslice_warp, _apply_transforms,
                   mask_report)
from .perfusion import coregister as coregister_perf
from .func import _realign, _slice_time
from .func import coregister as coregister_func
from .struct import anats_to_template
from sklearn.base import BaseEstimator, TransformerMixin


class AnatCoregistrator(BaseEstimator, TransformerMixin):
    """
    Encapsulation for anatomical data, relative to registration to a template.

    Parameters
    ----------
    anat : str
        Path to the anatomical image.

    brain_volume : int
        Volume of the brain used for brain extraction.
        Typically 400 for mouse and 1650 for rat.

    output_dir : str, optional
        Path to the output directory. If not specified, current directory is
        used.

    caching : bool, optional
        If True, caching is used for all the registration steps.

    verbose : int, optional
        Verbosity level. Note that caching implies some
        verbosity in any case.

    mask_clipping_fraction : float, optional
        Clip level fraction is passed to
        sammba.externals.nipype.interfaces.afni.Unifize, to tune
        the bias correction step done prior to brain mask segementation.
        Only values between 0.1 and 0.9 are accepted. Smaller fractions tend to
        make the mask larger.
    """
    def __init__(self, anat=None, brain_volume=None,
                 output_dir=None, caching=False,
                 verbose=True, use_rats_tool=True,
                 mask_clipping_fraction=.1):
        self.anat = anat
        self.brain_volume = brain_volume
        self.output_dir = output_dir
        self.use_rats_tool = use_rats_tool
        self.caching = caching
        self.verbose = verbose
        self.mask_clipping_fraction = mask_clipping_fraction

    def _check_inputs(self):
        if not os.path.isfile(self.anat):
            raise IOError('`anat` must be an existing '
                          'image file, you gave {0}'.format(self.anat))

        if self.brain_volume is None:
            raise IOError('`brain_volume` must be provided')

        if self.mask_clipping_fraction is not None:
            if (self.mask_clipping_fraction < .1 or
                    self.mask_clipping_fraction > .9):
                raise ValueError("'mask_clipping_fraction' must be between 0.1"
                                 "and 0.9, you provided {}"
                                 "".format(self.mask_clipping_fraction))

    def _set_output_dir(self):
        if self.output_dir is None:
            self.output_dir = os.getcwd()

        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)

    def check_segmentation(self):
        self.fit()
        print(mask_report(self.brain_, self.brain_volume))

    def segment(self):
        unifized_anat_file = _afni_bias_correct(
            self.anat, write_dir=self.output_dir,
            terminal_output=self.terminal_output, caching=self.caching)
        setattr(self, '_unifized_anat', unifized_anat_file)

        if self.mask_clipping_fraction:
            brain_mask_file = compute_brain_mask(
                self.anat, self.brain_volume, write_dir=self.output_dir,
                caching=self.caching,
                terminal_output=self.terminal_output,
                use_rats_tool=self.use_rats_tool,
                cl_frac=self.mask_clipping_fraction)
        else:
            brain_mask_file = compute_brain_mask(
                self.anat, self.brain_volume, write_dir=self.output_dir,
                caching=self.caching,
                terminal_output=self.terminal_output,
                use_rats_tool=self.use_rats_tool,
                bias_correct=False)

        brain_file = _apply_mask(unifized_anat_file, brain_mask_file,
                                 write_dir=self.output_dir,
                                 caching=self.caching,
                                 terminal_output=self.terminal_output)
        setattr(self, 'brain_', brain_file)

    def fit(self, y=None):
        self._check_inputs()
        if self.verbose:
            self.terminal_output = 'stream'
        else:
            self.terminal_output = 'none'
        self._set_output_dir()
        self.segment()
        return self

    def _check_fitted(self):
        if not hasattr(self, 'brain_'):
            raise ValueError('It seems that %s has not been fitted. '
                             'You must call fit() before calling transform().'
                             % self.__class__.__name__)

    def transform(self, in_file, modality, slice_timing=True,
                  t_r=None, prior_rigid_body_registration=False):
        """ Prepare and perform coregistration.

        Parameters
        ----------
        in_file : str
            Path to the raw modality image.

        modality : one of {'perf', 'func'}
            Name of the MRI modality.

        slice_timing : bool, optional
            If True, slice timing correction is performed

        t_r : float, optional
            Repetition time, only needed for slice timing correction.

        prior_rigid_body_registration : bool, optional
            If True, a rigid-body registration of the anat to the modality is
            performed prior to the warp. Useful if the images headers have
            missing/wrong information.

        Returns
        -------
        out_file : str
            Path to the modality image after coregistration.
        """
        self._check_fitted()
        for attribute in ['modality_brain_', '_unbiased_modality', 'modality']:
            if hasattr(self, attribute):
                delattr(self, attribute)

        if modality == 'perf':
            to_coregister_file = in_file
        elif modality == 'func':
            if slice_timing:
                # Correct functional for slice timing
                if t_r is None:
                    raise ValueError("'t_r' is needed for slice timing "
                                     "correction")
                func_file = _slice_time(in_file, t_r, self.output_dir,
                                        caching=self.caching)
            else:
                func_file = in_file

            # Register functional volumes to the first one #
            allineated_file, mean_aligned_file, _ = \
                _realign(func_file, write_dir=self.output_dir,
                         caching=self.caching,
                         terminal_output=self.terminal_output)
            to_coregister_file = mean_aligned_file
        else:
            raise ValueError("Only 'func' and 'perf' modalities are "
                             "implemented")

        setattr(self, 'modality', modality)
        unbiased_file = _bias_correct(to_coregister_file,
                                      write_dir=self.output_dir,
                                      terminal_output=self.terminal_output,
                                      caching=self.caching)

        if prior_rigid_body_registration:
            if self.mask_clipping_fraction:
                brain_mask_file = compute_brain_mask(
                    to_coregister_file, self.brain_volume,
                    write_dir=self.output_dir,
                    caching=self.caching,
                    terminal_output=self.terminal_output,
                    use_rats_tool=self.use_rats_tool,
                    cl_frac=self.mask_clipping_fraction)
            else:
                brain_mask_file = compute_brain_mask(
                    to_coregister_file, self.brain_volume,
                    write_dir=self.output_dir,
                    caching=self.caching,
                    terminal_output=self.terminal_output,
                    use_rats_tool=self.use_rats_tool,
                    bias_correct=False)

            brain_file = _apply_mask(unbiased_file, brain_mask_file,
                                     write_dir=self.output_dir,
                                     caching=self.caching,
                                     terminal_output=self.terminal_output)
            setattr(self, 'modality_brain_', brain_file)
        else:
            setattr(self, 'modality_brain_', None)

        if modality == 'func':
            coregistration = coregister_func(
                self._unifized_anat, unbiased_file,
                self.output_dir,
                anat_brain_file=self.brain_,
                func_brain_file=self.modality_brain_,
                prior_rigid_body_registration=prior_rigid_body_registration,
                caching=self.caching)
            setattr(self, '_coreg_warps', coregistration.coreg_warps_)
            coreg_modality_file = _apply_perslice_warp(
                allineated_file, self._coreg_warps, .1, .1,
                write_dir=self.output_dir, caching=self.caching)
        elif modality == 'perf':
            coregistration = coregister_perf(
                self._unifized_anat, unbiased_file,
                self.output_dir,
                anat_brain_file=self.brain_,
                m0_brain_file=self.modality_brain_,
                prior_rigid_body_registration=prior_rigid_body_registration,
                caching=self.caching)
            setattr(self, '_coreg_warps', coregistration.coreg_warps_)
            coreg_modality_file = coregistration.coreg_m0_

        setattr(self, 'coreg_anat_', coregistration.coreg_anat_)
        setattr(self, '_coreg_transform', coregistration.coreg_transform_)

        return coreg_modality_file

    def fit_transform(self, in_file, modality, slice_timing=None,
                      t_r=None, prior_rigid_body_registration=None,
                      **fit_params):
        """Fit to data, then transform it
#        Prepare and perform coregistration.

        Parameters
        ----------
        in_file : str
            Path to the raw modality image.

        modality : one of {'perf', 'func'}
            Name of the MRI modality.

        slice_timing : bool, optional
            If True, slice timing correction is performed

        t_r : float, optional
            Repetition time, only needed for slice timing correction.

        prior_rigid_body_registration : bool, optional
            If True, a rigid-body registration of the anat to the modality is
            performed prior to the warp. Useful if the images headers have
            missing/wrong information.

        Returns
        -------
        out_file : str
            Path to the modality image after coregistration.
        """
        return self.fit(**fit_params).transform(
            in_file,
            modality=modality,
            slice_timing=slice_timing,
            t_r=t_r,
            prior_rigid_body_registration=prior_rigid_body_registration)


class TemplateRegistrator(BaseEstimator, TransformerMixin):
    """
    Encapsulation for anatomical data, relative to registration to a template.

    Parameters
    ----------
    anat : str
        Path to the anatomical image.

    head_template : str
        Path to the template image.

    brain_volume : int
        Volume of the brain used for brain extraction.
        Typically 400 for mouse and 1650 for rat.

    brain_template_mask : str, optional
        Path to the template brain mask image.

    dilated_template_mask : str, optional
        Path to a dilated head mask. Note that this must be compliant with the
        the given head template. If None, the mask is set to the non-background
        voxels of the head template after one dilation.

    output_dir : str, optional
        Path to the output directory. If not specified, current directory is
        used.

    caching : bool, optional
        If True, caching is used for all the registration steps.

    verbose : int, optional
        Verbosity level. Note that caching implies some
        verbosity in any case.

    clip_level_fraction : float, optional
        Clip level fraction is passed to
        sammba.externals.nipype.interfaces.afni.Unifize, to tune
        the bias correction step done prior to brain mask segementation.
        Only values between 0.1 and 0.9 are accepted. Smaller fractions tend to
        make the mask larger.
    """

    def __init__(self, template=None,
                 brain_extracted_template=None, brain_volume=None,
                 dilated_template_mask=None, output_dir=None, caching=False,
                 verbose=True, use_rats_tool=True,
                 mask_clipping_fraction=.1):
        self.template = template
        self.brain_extracted_template = brain_extracted_template
        self.dilated_template_mask = dilated_template_mask
        self.brain_volume = brain_volume
        self.output_dir = output_dir
        self.use_rats_tool = use_rats_tool
        self.caching = caching
        self.verbose = verbose
        self.mask_clipping_fraction = mask_clipping_fraction

    def _check_inputs(self):
        if not os.path.isfile(self.template):
            raise IOError('`template` must be an existing '
                          'image file, you gave {0}'.format(self.template))

        if self.brain_volume is None:
            raise IOError('`brain_volume` must be provided')

        if self.mask_clipping_fraction is not None:
            if (self.mask_clipping_fraction < .1 or
                    self.mask_clipping_fraction > .9):
                raise ValueError("'mask_clipping_fraction' must be between 0.1"
                                 "and 0.9, you provided {}"
                                 "".format(self.mask_clipping_fraction))

    def _check_fitted(self):
        if not hasattr(self, 'terminal_output'):
            raise ValueError('It seems that %s has not been fitted. '
                             'You must call fit() before calling transform().'
                             % self.__class__.__name__)

    def _set_output_dir(self):
        if self.output_dir is None:
            self.output_dir = os.getcwd()

        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)

    def check_segmentation(self, in_file):
        _, brain_file = self.segment(in_file)
        print(mask_report(brain_file, self.brain_volume))

    def fit(self, y=None):
        self._check_inputs()
        self._set_output_dir()
        if not self.brain_extracted_template:
            brain_extracted_template = self.segment(self.template,
                                                    unifize=False)
            setattr(self, 'brain_extracted_template', brain_extracted_template)

        if self.verbose:
            self.terminal_output = 'stream'
        else:
            self.terminal_output = 'none'

        return self

    def segment(self, in_file, unifize=True):
        if unifize:
            file_to_mask = _afni_bias_correct(
                in_file, write_dir=self.output_dir,
                terminal_output=self.terminal_output, caching=self.caching)
        else:
            file_to_mask = in_file

        if self.mask_clipping_fraction:
            brain_mask_file = compute_brain_mask(
                in_file, self.brain_volume, write_dir=self.output_dir,
                caching=self.caching,
                terminal_output=self.terminal_output,
                use_rats_tool=self.use_rats_tool,
                cl_frac=self.mask_clipping_fraction)
        else:
            brain_mask_file = compute_brain_mask(
                in_file, self.brain_volume, write_dir=self.output_dir,
                caching=self.caching,
                terminal_output=self.terminal_output,
                use_rats_tool=self.use_rats_tool,
                bias_correct=False)

        brain_file = _apply_mask(file_to_mask, brain_mask_file,
                                 write_dir=self.output_dir,
                                 caching=self.caching,
                                 terminal_output=self.terminal_output)
        return file_to_mask, brain_file

    def transform(self, anat_file, y=None):
        """ Estimates noramlization from anatomical to template space.
        """
        self._check_fitted()

        for attribute in ['anat_', '_unifized_anat_', 'anat_brain_',
                          '_normalization_transform',
                          '_normalization_pretransform']:
            if hasattr(self, attribute):
                delattr(self, attribute)

        setattr(self, 'anat_', anat_file)
        file_to_mask, brain_file = self.segment(anat_file)
        setattr(self, '_unifized_anat_', file_to_mask)
        setattr(self, 'anat_brain_', brain_file)

        normalization = anats_to_template(
            [self._unifized_anat_], [self.anat_brain_], self.template,
            self.brain_extracted_template, write_dir=self.output_dir,
            dilated_head_mask_filename=self.dilated_template_mask,
            caching=self.caching, verbose=self.verbose, maxlev=0)  # XXX maxlev

        setattr(self, '_normalization_pretransform',
                normalization.pretransforms[0])
        setattr(self, '_normalization_transform',
                normalization.transforms[0])
        return normalization.registered[0]

    def inverse_transform(self, in_file, interpolation='wsinc5'):
        """Use provided loadings to compute corresponding linear component
        combination in whole-brain voxel space
        Parameters
        ----------
        loadings: list of numpy array (n_samples x n_components)
            Component signals to tranform back into voxel signals
        Returns
        -------
        reconstructed_imgs: list of nibabel.Nifti1Image
            For each loading, reconstructed Nifti1Image
        """
        self._check_fitted()
        self._check_transformed()

        inverted_in_file = _apply_transforms(in_file, self.anat_,
                                             self.output_dir,
                                             [self._normalization_pretransform,
                                              self._normalization_transform],
                                             inverse=True,
                                             interpolation=interpolation,
                                             caching=self.caching,
                                             verbose=self.verbose)
        return inverted_in_file

    def transform_modality(self, in_file, modality,
                           slice_timing=True, t_r=None,
                           prior_rigid_body_registration=False,
                           voxel_size=None):
        self._check_fitted()
        self._check_transformed()

        for attribute in ['modality_brain_', '_unbiased_modality', 'modality']:
            if hasattr(self, attribute):
                delattr(self, attribute)

        if modality == 'perf':
            to_coregister_file = in_file
        elif modality == 'func':
            if slice_timing:
                # Correct functional for slice timing
                if t_r is None:
                    raise ValueError("'t_r' is needed for slice timing "
                                     "correction")
                func_file = _slice_time(in_file, t_r, self.output_dir,
                                        caching=self.caching)
            else:
                func_file = in_file

            # Register functional volumes to the first one #
            allineated_file, mean_aligned_file, _ = \
                _realign(func_file, write_dir=self.output_dir,
                         caching=self.caching,
                         terminal_output=self.terminal_output)
            to_coregister_file = mean_aligned_file
        else:
            raise ValueError("Only 'func' and 'perf' modalities are "
                             "implemented")

        setattr(self, 'modality', modality)
        unbiased_file = _bias_correct(to_coregister_file,
                                      write_dir=self.output_dir,
                                      terminal_output=self.terminal_output,
                                      caching=self.caching)

        if prior_rigid_body_registration:
            if self.mask_clipping_fraction:
                brain_mask_file = compute_brain_mask(
                    to_coregister_file, self.brain_volume,
                    write_dir=self.output_dir,
                    caching=self.caching,
                    terminal_output=self.terminal_output,
                    use_rats_tool=self.use_rats_tool,
                    cl_frac=self.mask_clipping_fraction)
            else:
                brain_mask_file = compute_brain_mask(
                    to_coregister_file, self.brain_volume,
                    write_dir=self.output_dir,
                    caching=self.caching,
                    terminal_output=self.terminal_output,
                    use_rats_tool=self.use_rats_tool,
                    bias_correct=False)

            brain_file = _apply_mask(unbiased_file, brain_mask_file,
                                     write_dir=self.output_dir,
                                     caching=self.caching,
                                     terminal_output=self.terminal_output)
            setattr(self, 'modality_brain_', brain_file)
        else:
            setattr(self, 'modality_brain_', None)

        if modality == 'func':
            coregistration = coregister_func(
                self._unifized_anat_, unbiased_file,
                self.output_dir,
                anat_brain_file=self.anat_brain_,
                func_brain_file=self.modality_brain_,
                prior_rigid_body_registration=prior_rigid_body_registration,
                caching=self.caching)
            setattr(self, '_coreg_warps', coregistration.coreg_warps_)
            coreg_modality_file = _apply_perslice_warp(
                allineated_file, self._coreg_warps, .1, .1,
                write_dir=self.output_dir, caching=self.caching)
        elif modality == 'perf':
            coregistration = coregister_perf(
                self._unifized_anat_, unbiased_file,
                self.output_dir,
                anat_brain_file=self.anat_brain_,
                m0_brain_file=self.modality_brain_,
                prior_rigid_body_registration=prior_rigid_body_registration,
                caching=self.caching)
            setattr(self, '_coreg_warps', coregistration.coreg_warps_)
            coreg_modality_file = coregistration.coreg_m0_

        setattr(self, 'coreg_modality_', coreg_modality_file)
        setattr(self, 'coreg_anat_', coregistration.coreg_anat_)
        setattr(self, '_coreg_transform', coregistration.coreg_transform_)

        normalized_file = _transform_to_template(
            coreg_modality_file, self.template, self.output_dir,
            self._coreg_transform, self._normalization_pretransform,
            self._normalization_transform,
            voxel_size=voxel_size, caching=self.caching)

        return normalized_file

    def _check_transformed(self):
        if not hasattr(self, '_normalization_transform'):
            raise ValueError(
                'It seems that %s has not been transformed. '
                'You must call transform() before calling apply_to().'
                % self.__class__.__name__)

    def inverse_transform_modality(self, in_file, interpolation='wsinc5'):
        """Use provided loadings to compute corresponding linear component
        combination in whole-brain voxel space
        Parameters
        ----------
        loadings: list of numpy array (n_samples x n_components)
            Component signals to tranform back into voxel signals
        Returns
        -------
        reconstructed_imgs: list of nibabel.Nifti1Image
            For each loading, reconstructed Nifti1Image
        """
        self._check_fitted()
        self._check_transformed()
        if not hasattr(self, '_coreg_transform'):
            raise ValueError(
                'It seems that %s has not been modality transformed. '
                'You must call transform_modality() before calling '
                'inverse_transform_modality().'
                % self.__class__.__name__)
        inverted_in_file = _apply_transforms(in_file, self.coreg_anat_,
                                             self.output_dir,
                                             [self._coreg_transform,
                                              self._normalization_pretransform,
                                              self._normalization_transform],
                                             inverse=True,
                                             interpolation=interpolation,
                                             caching=self.caching,
                                             verbose=self.verbose)
        return inverted_in_file
