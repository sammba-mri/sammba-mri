import warnings
from nilearn._utils.exceptions import VisibleDeprecationWarning
from ..segmentation.brain_mask import (compute_histo_brain_mask,
                                       compute_morpho_brain_mask,
                                       _apply_mask)
from ..preprocessing.bias_correction import ants_n4, afni_unifize
from .base import _apply_perslice_warp
from .perfusion import coregister as coregister_perf
from .func import _realign, _slice_time
from .func import coregister as coregister_func
from .base_registrator import BaseRegistrator


class Coregistrator(BaseRegistrator):
    """
    Class for registering anatomical image to perfusion/functional images from
    one animal in native space.

    Parameters
    ----------
    brain_volume : int or None, optional
        Volume of the brain in mm3 used for brain extraction. Typically 400 for mouse
        and 1650 for rat. Used only if prior rigid body registration is needed.

    output_dir : str or None, optional
        Path to the output directory. If None, current directory is used.

    caching : bool, optional
        If True, caching is used for all the registration steps.

    verbose : int, optional
        Verbosity level. Note that caching implies some verbosity in any case.

    use_rats_tool : bool, optional
        If True, brain mask is computed using RATS Mathematical Morphology.
        Otherwise, a histogram-based brain segmentation is used.

    clipping_fraction : float or None, optional
        Clip level fraction is passed to
        nipype.interfaces.afni.Unifize, to tune
        the bias correction step done prior to brain mask segmentation.
        Only values between 0.1 and 0.9 are accepted. Smaller fractions tend to
        make the mask larger. If None, no unifization is done for brain mask
        computation.
    """
    def __init__(self, brain_volume=None,
                 output_dir=None, caching=False,
                 verbose=True, use_rats_tool=True,
                 clipping_fraction=.2):
        self.brain_volume = brain_volume
        self.output_dir = output_dir
        self.use_rats_tool = use_rats_tool
        self.caching = caching
        self.verbose = verbose
        self.clipping_fraction = clipping_fraction

    def _check_inputs(self):
        if self.clipping_fraction is not None:
            if (self.clipping_fraction < .1 or
                    self.clipping_fraction > .9):
                raise ValueError("'clipping_fraction' must be between 0.1"
                                 "and 0.9, you provided {}"
                                 "".format(self.clipping_fraction))

    def fit_anat(self, anat_file, brain_mask_file=None):
        self._fit()
        self.anat_ = anat_file
        self._anat_brain_mask = brain_mask_file
        return self

    def _check_anat_fitted(self):
        if not hasattr(self, '_anat_brain_mask'):
            raise ValueError(
                'It seems that %s has not been anat fitted. You must call '
                'fit_anat() before calling fit_modality().'
                % self.__class__.__name__)

    def fit_modality(self, in_file, modality, slice_timing=True,
                     t_r=None, prior_rigid_body_registration=None,
                     reorient_only=False, brain_mask_file=None):
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
            NOTE: prior_rigid_body_registration is deprecated from 0.1 and
            will be removed in next release. Use `reorient_only` instead.

        reorient_only :  bool, optional
            If True, the rigid-body registration of the anat to the func is
            not performed and only reorientation is done.

        Returns
        -------
        the coregistrator itself
        """
        if prior_rigid_body_registration is not None:
            warn_str = ("The parameter 'prior_rigid_body_registration' is "
                        "deprecated and will be removed in sammba-mri next "
                        "release. Use parameter 'reorient_only' instead.")
            warnings.warn(warn_str, VisibleDeprecationWarning, stacklevel=2)
            reorient_only = not(prior_rigid_body_registration)

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

        unbiased_file = ants_n4(
            to_coregister_file,
            write_dir=self.output_dir,
            terminal_output=self.terminal_output,
            caching=self.caching)

        if reorient_only:
            unifized_anat_file = afni_unifize(
                self.anat_, write_dir=self.output_dir,
                terminal_output=self.terminal_output, caching=self.caching,
                verbose=self.verbose)
            self.anat_brain_ = None
            modality_brain_file = None
        else:
            if brain_mask_file is None or self._anat_brain_mask is None:
                if not isinstance(self.brain_volume, int):
                    raise ValueError('`brain_volume` must be specified to '
                                     'perform rigid-body registration')
            if self.use_rats_tool:
                compute_brain_mask = compute_morpho_brain_mask
            else:
                compute_brain_mask = compute_histo_brain_mask               

            if brain_mask_file is None:
                if self.clipping_fraction:
                    brain_mask_file = compute_brain_mask(
                        to_coregister_file, self.brain_volume,
                        write_dir=self.output_dir,
                        caching=self.caching,
                        terminal_output=self.terminal_output,
                        cl_frac=self.clipping_fraction,
                        verbose=self.verbose)
                else:
                    brain_mask_file = compute_brain_mask(
                        to_coregister_file, self.brain_volume,
                        write_dir=self.output_dir,
                        caching=self.caching,
                        terminal_output=self.terminal_output,
                        unifize=False,
                        verbose=self.verbose)

            modality_brain_file = _apply_mask(
                unbiased_file, brain_mask_file, write_dir=self.output_dir,
                caching=self.caching, terminal_output=self.terminal_output)
            if self._anat_brain_mask is None:
                unifized_anat_file, self.anat_brain_ = self.segment(
                    self.anat_)
            else:
                unifized_anat_file = afni_unifize(
                    self.anat_, write_dir=self.output_dir,
                    terminal_output=self.terminal_output, caching=self.caching,
                    verbose=self.verbose)
                self.anat_brain_ = _apply_mask(
                    unifized_anat_file, self._anat_brain_mask,
                    write_dir=self.output_dir,
                    caching=self.caching,
                    terminal_output=self.terminal_output)

        if modality == 'func':
            self.func_brain_ = modality_brain_file
            coregistration = coregister_func(
                unifized_anat_file, unbiased_file,
                self.output_dir,
                anat_brain_file=self.anat_brain_,
                func_brain_file=modality_brain_file,
                reorient_only=reorient_only,
                caching=self.caching,
                verbose=self.verbose)
            self._func_undistort_warps = coregistration.coreg_warps_
            self.anat_in_func_space_ = coregistration.coreg_anat_
            self._func_to_anat_transform = coregistration.coreg_transform_
            self.undistorted_func_ = _apply_perslice_warp(
                allineated_file, self._func_undistort_warps, .1, .1,
                write_dir=self.output_dir, caching=self.caching)
        elif modality == 'perf':
            self.perf_brain_ = modality_brain_file
            coregistration = coregister_perf(
                unifized_anat_file, unbiased_file,
                self.output_dir,
                anat_brain_file=self.anat_brain_,
                m0_brain_file=modality_brain_file,
                reorient_only=reorient_only,
                caching=self.caching,
                verbose=self.verbose)
            self.undistorted_perf_ = coregistration.coreg_m0_
            self._perf_undistort_warps = coregistration.coreg_warps_
            self.anat_in_perf_space_ = coregistration.coreg_anat_
            self._perf_to_anat_transform = coregistration.coreg_transform_

        return self

    def transform_modality_like(self, apply_to_file, modality):
        """ Applies modality coregristration to a file in the modality space.
        """
        if modality not in ['perf', 'func']:
            raise ValueError("Only 'func' and 'perf' modalities are "
                             "implemented")

        self._check_anat_fitted()
        modality_undistort_warps = '_{}_undistort_warps'.format(modality)
        if not hasattr(self, modality_undistort_warps):
            raise ValueError('It seems that {0} has not been {1} fitted. '
                             'You must call fit_modality() before calling '
                             'transform_modality().'.format(
                                     self.__class__.__name__, modality))

        coreg_apply_to_file = _apply_perslice_warp(
            apply_to_file, self.__getattribute__(modality_undistort_warps),
            .1, .1, write_dir=self.output_dir, caching=self.caching,
            verbose=self.verbose)
        return coreg_apply_to_file
