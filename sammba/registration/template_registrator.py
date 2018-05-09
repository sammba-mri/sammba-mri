import os
from .base import (compute_brain_mask, _bias_correct,
                   _afni_bias_correct, _apply_mask,
                   _apply_perslice_warp, _apply_transforms)
from .perfusion import coregister as coregister_perf
from .func import _realign, _slice_time
from .func import coregister as coregister_func
from .struct import anats_to_template
from .base_registrator import BaseRegistrator


class TemplateRegistrator(BaseRegistrator):
    """
    Encapsulation for anatomical data, relative to registration to a template.

    Parameters
    ----------
    template : str
        Path to the head template image.

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

    mask_clipping_fraction : float, optional
        Clip level fraction is passed to
        sammba.externals.nipype.interfaces.afni.Unifize, to tune
        the bias correction step done prior to brain mask segementation.
        Only values between 0.1 and 0.9 are accepted. Smaller fractions tend to
        make the mask larger.
    """

    def __init__(self, template=None,
                 template_brain_mask=None, brain_volume=None,
                 dilated_template_mask=None, output_dir=None, caching=False,
                 verbose=True, use_rats_tool=True,
                 mask_clipping_fraction=.2, convergence=0.005):
        self.template = template
        self.template_brain_mask = template_brain_mask
        self.dilated_template_mask = dilated_template_mask
        self.brain_volume = brain_volume
        self.output_dir = output_dir
        self.use_rats_tool = use_rats_tool
        self.caching = caching
        self.verbose = verbose
        self.mask_clipping_fraction = mask_clipping_fraction
        self.convergence = convergence

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

    def _check_anat_fitted(self):
        if not hasattr(self, '_normalization_transform'):
            raise ValueError('It seems that %s has not been fitted. '
                             'You must call fit_anat() before calling '
                             'transform_anat() or fit_modality().'
                             % self.__class__.__name__)

    def fit_anat(self, anat_file, brain_mask_file=None):
        """ Estimates registration from anatomical to template space.
        """
        self._fit()
        if not self.template_brain_mask:
            _, self.template_brain_ = self.segment(self.template,
                                                   unifize=False)
        else:
            self.template_brain_ = _apply_mask(
                self.template, self.template_brain_mask,
                write_dir=self.output_dir,
                caching=self.caching,
                terminal_output=self.terminal_output)

        self.anat_ = anat_file
        if brain_mask_file is None:
            self._unifized_anat_, self.anat_brain_ = self.segment(self.anat_)
        else:
            self._unifized_anat_ = _afni_bias_correct(
                self.anat_, write_dir=self.output_dir,
                terminal_output=self.terminal_output, caching=self.caching)
            self.anat_brain_ = _apply_mask(
                self._unifized_anat_, brain_mask_file,
                write_dir=self.output_dir,
                caching=self.caching,
                terminal_output=self.terminal_output)
        
        normalization = anats_to_template(
            [self._unifized_anat_], [self.anat_brain_], self.template,
            self.template_brain_, write_dir=self.output_dir,
            dilated_head_mask_filename=self.dilated_template_mask,
            caching=self.caching, verbose=self.verbose, maxlev=None,
            convergence=self.convergence)  # XXX maxlev

        self.registered_anat = normalization.registered[0]
        self._normalization_pretransform = normalization.pretransforms[0]
        self._normalization_transform = normalization.transforms[0]
        return self

    def transform_anat(self, in_file, interpolation='wsinc5'):
        """ Transforms the given in_file from anatomical space to template
            space.
        """
        self._check_anat_fitted()
        transformed_file = _apply_transforms(in_file, self.anat_,
                                             self.output_dir,
                                             [self._normalization_transform,
                                              self._normalization_pretransform],
                                             interpolation=interpolation,
                                             caching=self.caching,
                                             verbose=self.verbose)
        return transformed_file

    def fit_transform_anat(self, anat_file, apply_to_file, **fit_params):
        """ Estimates registration from anatomical space to template space
            then transforms the given file.
        """
        return self.fit(anat_file, **fit_params).transform(apply_to_file)

    def _check_modality_fitted(self, modality):
        coreg_transform = '_{}_to_anat_transform'.format(modality)
        if not hasattr(self, coreg_transform):
            raise ValueError('It seems that %s has not been fitted. '
                             'You must call fit_modality() before calling '
                             'transform_modality().' % self.__class__.__name__)

    def fit_modality(self, in_file, modality, slice_timing=True, t_r=None,
                     prior_rigid_body_registration=False, voxel_size=None):
        """ Estimates registration from the space of a given modality to
            the template space.
        """
        self._check_anat_fitted()

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
        else:
            brain_file = None

        if modality == 'func':
            self.func_brain_ = brain_file
            coregistration = coregister_func(
                self._unifized_anat_, unbiased_file,
                self.output_dir,
                anat_brain_file=self.anat_brain_,
                func_brain_file=brain_file,
                prior_rigid_body_registration=prior_rigid_body_registration,
                caching=self.caching)
            self._func_undistort_warps = coregistration.coreg_warps_
            self.anat_in_func_space = coregistration.coreg_anat_
            self._func_to_anat_transform = coregistration.coreg_transform_
            self.undistorted_func = _apply_perslice_warp(
                allineated_file, self._func_undistort_warps, .1, .1,
                write_dir=self.output_dir, caching=self.caching)

            self.registered_func_ = _apply_transforms(
                self.undistorted_func, self.template, self.output_dir,
                [self._normalization_transform,
                 self._normalization_pretransform,
                 self._func_to_anat_transform],
                voxel_size=voxel_size, caching=self.caching)
        elif modality == 'perf':
            self.perf_brain_ = brain_file
            coregistration = coregister_perf(
                self._unifized_anat_, unbiased_file,
                self.output_dir,
                anat_brain_file=self.anat_brain_,
                m0_brain_file=brain_file,
                prior_rigid_body_registration=prior_rigid_body_registration,
                caching=self.caching)
            self.undistorted_perf = coregistration.coreg_m0_
            self._perf_undistort_warps = coregistration.coreg_warps_
            self.anat_in_perf_space = coregistration.coreg_anat_
            self._perf_to_anat_transform = coregistration.coreg_transform_

            self.registered_perf_ = _apply_transforms(
                self.undistorted_perf, self.template, self.output_dir,
                [self._normalization_transform,
                 self._normalization_pretransform,
                 self._perf_to_anat_transform],
                voxel_size=voxel_size, caching=self.caching)

        return self

    def transform_modality(self, in_file, modality, interpolation='wsinc5',
                           voxel_size=None):
        """Transforms the given file from the space of the given modality to
            the template space then .
        """
        modality_undistort_warps = self.__getattribute__(
            '_{}_undistort_warps'.format(modality))
        coreg_transform_file = self.__getattribute__(
            '_{}_to_anat_transform'.format(modality))
        undistorted_file = _apply_perslice_warp(in_file,
                                                modality_undistort_warps,
                                                .1,
                                                .1,
                                                write_dir=self.output_dir,
                                                caching=self.caching)
        normalized_file = _apply_transforms(
            undistorted_file, self.template, self.output_dir,
            [self._normalization_transform, self._normalization_pretransform,
             coreg_transform_file],
            voxel_size=voxel_size, caching=self.caching)
        return normalized_file

    def fit_transform_modality(self, modality_file, modality, apply_to_file,
                               interpolation='wsinc5',
                               voxel_size=None, **fit_params):
        """ Estimates registration from the space of the given modality to the
            template space then transforms the given file.
        """
        self._check_anat_fitted()
        return self.fit_modality(
            modality_file, modality, **fit_params).transform_modaliy(
            apply_to_file, interpolation=interpolation, voxel_size=voxel_size)

    def inverse_transform_modality(self, in_file, modality,
                                   interpolation='wsinc5'):
        """ Trasnforms the given file from template space to modality space.
        """
        self._check_anat_fitted()
        self._check_modality_fitted()
        coreg_transform_file = self.__getattribute__(
            '_{}_to_anat_transform'.format(modality))
        modality_file = self.__getattribute__(modality)

        inverted_file = _apply_transforms(in_file, modality_file,
                                          self.output_dir,
                                          [self._normalization_transform,
                                           self._normalization_pretransform,
                                           coreg_transform_file],
                                          inverse=True,
                                          interpolation=interpolation,
                                          caching=self.caching,
                                          verbose=self.verbose)
        return inverted_file
