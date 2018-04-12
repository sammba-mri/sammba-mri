import os
from .base import (compute_brain_mask, _bias_correct,
                   _afni_bias_correct, _apply_mask,
                   _apply_perslice_warp,
                   mask_report)
from .perfusion import coregister as coregister_perf
from .func import _realign, _slice_time
from .func import coregister as coregister_func
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
    def __init__(self, brain_volume=None,
                 output_dir=None, caching=False,
                 verbose=True, use_rats_tool=True,
                 mask_clipping_fraction=.2):
        self.brain_volume = brain_volume
        self.output_dir = output_dir
        self.use_rats_tool = use_rats_tool
        self.caching = caching
        self.verbose = verbose
        self.mask_clipping_fraction = mask_clipping_fraction

    def _check_inputs(self):
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
        self.output_dir = os.path.abspath(self.output_dir)

        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)

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

    def check_segmentation(self, in_file):
        """ Quality check mask computation for the chosen
            `mask_clipping_fraction` parameter.
        """
        self._fit()
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

        return mask_report(brain_mask_file, self.brain_volume)

    def _fit(self, y=None):
        self._check_inputs()
        if self.verbose:
            self.terminal_output = 'stream'
        else:
            self.terminal_output = 'none'
        self._set_output_dir()
        return self

    def fit_anat(self, anat_file, brain_mask_file=None):
        self._fit()
        self.anat_ = anat_file
        file_to_mask, brain_file = self.segment(anat_file)
        self._unifized_anat = file_to_mask
        self.anat_brain_ = brain_file
        return self

    def _check_anat_fitted(self):
        if not hasattr(self, 'anat_brain_'):
            raise ValueError(
                'It seems that %s has not been fitted. You must call '
                'fit_anat() before calling fit_modality().'
                % self.__class__.__name__)

    def fit_modality(self, in_file, modality, slice_timing=True,
                     t_r=None, prior_rigid_body_registration=False,
                     brain_mask_file=None):
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
            if brain_mask_file is None:
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
                self._unifized_anat, unbiased_file,
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
        elif modality == 'perf':
            self.perf_brain_ = brain_file
            coregistration = coregister_perf(
                self._unifized_anat, unbiased_file,
                self.output_dir,
                anat_brain_file=self.anat_brain_,
                m0_brain_file=brain_file,
                prior_rigid_body_registration=prior_rigid_body_registration,
                caching=self.caching)
            self.undistorted_perf = coregistration.coreg_m0_
            self._perf_undistort_warps = coregistration.coreg_warps_
            self.anat_in_perf_space = coregistration.coreg_anat_
            self._perf_to_anat_transform = coregistration.coreg_transform_

        return self

    def transform_modality(self, apply_to_file, modality):
        """ Applies modality coregristration to a file in the modality space.
        """
        self._check_anat_fitted()
        modality_coreg_warps = '_coreg_{}_warps'.format(modality)
        if not hasattr(self, modality_coreg_warps):
            raise ValueError('It seems that %s has not been transformed. '
                             'You must call fit_modality() before calling '
                             'transform_modality().' % self.__class__.__name__)

        coreg_apply_to_file = _apply_perslice_warp(
            apply_to_file, self.__getattribute__(modality_coreg_warps), .1, .1,
            write_dir=self.output_dir, caching=self.caching)
        return coreg_apply_to_file

    def fit_transform_modality(self, in_file, apply_to_file,
                               modality, slice_timing=None,
                               t_r=None, prior_rigid_body_registration=None,
                               **fit_params):
        """#Fit to data, then transform it
        Prepare and perform coregistration.

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
        self._check_anat_fitted()
        return self.fit_modality(
            in_file, modality, slice_timing=slice_timing, t_r=t_r,
            prior_rigid_body_registration=
            prior_rigid_body_registration).transform_modality(
            apply_to_file,
            modality=modality)