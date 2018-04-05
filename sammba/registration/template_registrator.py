import os
from .base import (BaseRegistrator, compute_brain_mask, _bias_correct,
                   _afni_bias_correct, _apply_mask, _transform_to_template,
                   _apply_perslice_warp)
from .perfusion import coregister as coregister_perf
from .func import _realign, _slice_time
from .func import coregister as coregister_func
from .struct import anats_to_template


class TemplateRegistrator(BaseRegistrator):
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

    def __init__(self, anat=None, template=None,
                 brain_extracted_template=None, brain_volume=None,
                 dilated_template_mask=None, output_dir=None, caching=False,
                 verbose=True, use_rats_tool=True,
                 clip_level_fraction=.1):
        self.anat = anat
        self.template = template
        self.brain_extracted_template = brain_extracted_template
        self.dilated_template_mask = dilated_template_mask
        self.brain_volume = brain_volume
        self.output_dir = output_dir
        self.use_rats_tool = use_rats_tool
        self.caching = caching
        self.verbose = verbose
        self.clip_level_fraction = clip_level_fraction

    def _check_inputs(self):
        if not os.path.isfile(self.anat):
            raise IOError('`anat` must be an existing '
                          'image file, you gave {0}'.format(self.anat))

        if not os.path.isfile(self.template):
            raise IOError('`template` must be an existing '
                          'image file, you gave {0}'.format(self.template))

        if self.brain_volume is None:
            raise IOError('`brain_volume` must be provided')

        if self.clip_level_fraction is not None:
            if self.clip_level_fraction < .1 or self.clip_level_fraction > .9:
                raise ValueError("'clip_level_fraction' must be between 0.1"
                                 "and 0.9, you provided {}".format(
                                 self.clip_level_fraction))

    def fit(self):
        self._check_inputs()
#        self._output_type = _get_afni_output_type(self.anat)
        if self.verbose:
            self.terminal_output = 'stream'
        else:
            self.terminal_output = 'none'
        self._set_output_dir()

    def segment(self):
        self.fit()
        unifized_anat_file = _afni_bias_correct(
            self.anat, write_dir=self.output_dir,
            terminal_output=self.terminal_output, caching=self.caching)
        setattr(self, '_unifized_anat', unifized_anat_file)

        if self.clip_level_fraction:
            brain_mask_file = compute_brain_mask(
                self.anat, self.brain_volume, write_dir=self.output_dir,
                caching=self.caching,
                terminal_output=self.terminal_output,
                use_rats_tool=self.use_rats_tool,
                cl_frac=self.clip_level_fraction)
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

    def normalize(self):
        """ Estimates noramlization from anatomical to template space.
        """
        if not hasattr(self, 'brain_'):
            raise ValueError('anatomical image has not been segmented')

        for attribute in ['normalized_anat', '_normalization_transform',
                          '_normalization_pretransform']:
            if hasattr(self, attribute):
                delattr(self, attribute)

        normalization = anats_to_template(
            [self._unifized_anat], [self.brain_], self.template,
            self.brain_extracted_template, write_dir=self.output_dir,
            dilated_head_mask_filename=self.dilated_template_mask,
            caching=self.caching, verbose=self.verbose, maxlev=0)

        setattr(self, 'normalized_anat', normalization.registered[0])
        setattr(self, '_normalization_pretransform',
                normalization.pretransforms[0])
        setattr(self, '_normalization_transform',
                normalization.transforms[0])
        return self

    def coregister_modality(self, in_file, modality, clip_level_fraction=None,
                            slice_timing=True, t_r=None,
                            prior_rigid_body_registration=False):
        for attribute in ['modality_brain_', '_unbiased_modality', 'modality']:
            if hasattr(self, attribute):
                delattr(self, attribute)

        setattr(self, 'modality', modality)
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
            if clip_level_fraction:
                brain_mask_file = compute_brain_mask(
                    to_coregister_file, self.brain_volume,
                    write_dir=self.output_dir,
                    caching=self.caching,
                    terminal_output=self.terminal_output,
                    use_rats_tool=self.use_rats_tool,
                    cl_frac=clip_level_fraction)
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

    def normalize_modality(self, modality_file, modality, voxel_size=None,
                           clip_level_fraction=None,
                           slice_timing=True, t_r=None,
                           prior_rigid_body_registration=False):
        """ Coregisters the anatomical and the given modality to the same space
        """
        if not hasattr(self, '_normalization_transform'):
            raise ValueError('anatomical image has not been normalized')

        coreg_modality_file = self.coregister_modality(
            modality_file, modality,
            clip_level_fraction=clip_level_fraction,
            slice_timing=slice_timing, t_r=t_r,
            prior_rigid_body_registration=prior_rigid_body_registration)
        normalized_file = _transform_to_template(
            coreg_modality_file, self.template, self.output_dir,
            self._coreg_transform, self._normalization_pretransform,
            self._normalization_transform,
            voxel_size=voxel_size, caching=self.caching)
        return normalized_file

    def apply_modality_coregistration(self, apply_to_file):
        """ Applies modality coregristration to a file in the modality space.
        """
        if not hasattr(self, '_coreg_warps'):
            raise ValueError('modality has not been coregistrated')

        coreg_apply_to_file = _apply_perslice_warp(apply_to_file,
                                                   self._coreg_warps,
                                                   .1,
                                                   .1,
                                                   write_dir=self.output_dir,
                                                   caching=self.caching)
        return coreg_apply_to_file

    def apply_modality_normalization(self, apply_to_file, voxel_size=None):
        """ Applies modality normalization to a file in the modality space.
        """
        if not hasattr(self, '_coreg_warps'):
            raise ValueError('modality has not been coregistrated')
        if not hasattr(self, '_normalization_transform'):
            raise ValueError('anatomical image has not been normalized')

        coreg_apply_to_file = _apply_perslice_warp(apply_to_file,
                                                   self._coreg_warps,
                                                   .1,
                                                   .1,
                                                   write_dir=self.output_dir,
                                                   caching=self.caching)
        normalized_file = _transform_to_template(
            coreg_apply_to_file, self.template, self.output_dir,
            self._coreg_transform, self._normalization_pretransform,
            self._normalization_transform,
            voxel_size=voxel_size, caching=self.caching)
        return normalized_file

    def inverse_normalize_modality(self, in_file):
        """ Applies inverse normalization from template space to modality space
        """
        if not hasattr(self, '_coref_transform'):
            raise ValueError('anatomical image has not been coregistered to'
                             'the {} space'.format(in_file))
        if not _check_coregistration(in_file, self.coreg_anat):
            raise ValueError('{0} and {1} are not coregistered'.format(
                in_file, self.coreg_anat_))

        inverse_normalized_file = _inverse_transform_to_template(
            in_file, self.template, self.output_dir
            [self._coreg_transform, self._normalization_transform])
        return inverse_normalized_file
