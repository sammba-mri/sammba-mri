import os
from ..segmentation.brain_mask import (compute_histo_brain_mask,
                                       compute_morpho_brain_mask,
                                       _apply_mask)
from ..preprocessing.bias_correction import ants_n4, afni_unifize
from .base import _apply_perslice_warp, _apply_transforms
from .perfusion import coregister as coregister_perf
from .func import _realign, _slice_time
from .func import coregister as coregister_func
from .struct import anat_to_template
from .base_registrator import BaseRegistrator


class TemplateRegistrator(BaseRegistrator):
    """
    Class for registering anatomical and possibly other modality images from
    one animal to a given head template.

    Parameters
    ----------
    template : str
        Path to the head template image.

    brain_volume : int
        Volume of the brain in mm3 used for brain extraction.
        Typically 400 for mouse and 1650 for rat.

    template_brain_mask : str or None, optional
        Path to the template brain mask image, compliant with the given head
        template.

    dilated_template_mask : str or None, optional
        Path to a dilated head mask, compliant with the given head template.
        If None, the mask is set to the non-background voxels of the head
        template after one dilation.

    output_dir : str, optional
        Path to the output directory. If not specified, current directory is
        used.

    caching : bool, optional
        If True, caching is used for all the registration steps.

    verbose : int, optional
        Verbosity level. Note that caching implies some
        verbosity in any case.

    use_rats_tool : bool, optional
        If True, brain mask is computed using RATS Mathematical Morphology.
        Otherwise, a histogram-based brain segmentation is used.

    clipping_fraction : float or None, optional
        Clip level fraction is passed to
        sammba.externals.nipype.interfaces.afni.Unifize, to tune
        the bias correction step done prior to brain mask segmentation.
        Only values between 0.1 and 0.9 are accepted. Smaller fractions tend to
        make the mask larger.
        If None, no unifization is done for brain mask computation.

    convergence : float, optional
        Convergence limit, passed to
        sammba.externals.nipype.interfaces.afni.Allineate

    registration_kind : one of {'rigid', 'affine', 'nonlinear'}, optional
        The allowed transform kind from the anatomical image to the template.

    Attributes
    ----------
    `template_brain_` : str
        Path to the brain extracted file from the template image

    `anat_brain_` : str
        Path to the brain extracted file from the anatomical image
    """
    def __init__(self, template, brain_volume,
                 template_brain_mask=None,
                 dilated_template_mask=None, output_dir=None, caching=False,
                 verbose=True, use_rats_tool=True,
                 clipping_fraction=.2, convergence=0.005,
                 registration_kind='nonlinear'):
        self.template = template
        self.template_brain_mask = template_brain_mask
        self.dilated_template_mask = dilated_template_mask
        self.brain_volume = brain_volume
        self.output_dir = output_dir
        self.use_rats_tool = use_rats_tool
        self.caching = caching
        self.verbose = verbose
        self.clipping_fraction = clipping_fraction
        self.convergence = convergence
        self.registration_kind = registration_kind

    def _check_inputs(self):
        if not os.path.isfile(self.template):
            raise IOError('`template` must be an existing '
                          'image file, you gave {0}'.format(self.template))

        if not isinstance(self.brain_volume, int):
            raise ValueError('`brain_volume` must be an integer')

        if self.clipping_fraction is not None:
            if (self.clipping_fraction < .1 or
                    self.clipping_fraction > .9):
                raise ValueError("'clipping_fraction' must be between 0.1"
                                 " and 0.9, you provided {}"
                                 "".format(self.clipping_fraction))

    def _check_anat_fitted(self):
        if not hasattr(self, '_normalization_transforms'):
            raise ValueError('It seems that %s has not been anat fitted. '
                             'You must call fit_anat() before calling '
                             'transform_anat() or fit_modality().'
                             % self.__class__.__name__)

    def fit_anat(self, anat_file, brain_mask_file=None):
        """ Estimates registration from anatomical to template space.
        """
        self._fit()
        if self.use_rats_tool:
            compute_brain_mask = compute_morpho_brain_mask
        else:
            compute_brain_mask = compute_histo_brain_mask               

        if not self.template_brain_mask:
            template_brain_mask_file = compute_brain_mask(
                self.template, self.brain_volume,
                write_dir=self.output_dir,
                caching=self.caching,
                terminal_output=self.terminal_output,
                unifize=False,
                verbose=self.verbose)
        else:
            template_brain_mask_file = self.template_brain_mask

        self.template_brain_ = _apply_mask(
            self.template, template_brain_mask_file,
            write_dir=self.output_dir,
            caching=self.caching,
            terminal_output=self.terminal_output)

        self.anat_ = anat_file
        if brain_mask_file is None:
            self._unifized_anat, self.anat_brain_ = self.segment(self.anat_)
        else:
            self._unifized_anat = afni_unifize(
                self.anat_, write_dir=self.output_dir,
                terminal_output=self.terminal_output, caching=self.caching,
                verbose=self.verbose)
            self.anat_brain_ = _apply_mask(
                self._unifized_anat, brain_mask_file,
                write_dir=self.output_dir,
                caching=self.caching,
                terminal_output=self.terminal_output)

        normalization = anat_to_template(
            self._unifized_anat, self.anat_brain_, self.template,
            self.template_brain_, write_dir=self.output_dir,
            dilated_head_mask_filename=self.dilated_template_mask,
            caching=self.caching, verbose=self.verbose, maxlev=None,
            convergence=self.convergence,
            registration_kind=self.registration_kind)

        self.registered_anat_ = normalization.registered
        if self.registration_kind == 'nonlinear':
            self._normalization_transforms = [normalization.transform,
                                              normalization.pretransform]
        else:
            self._normalization_transforms = [normalization.pretransform]

        return self

    def transform_anat_like(self, in_file, interpolation='wsinc5'):
        """ Transforms the given in_file from anatomical space to template
            space.

        Parameters
        ----------
        in_file: str
            Path to the file in the same space as the anatomical image.

        interpolation: one of {'nearestneighbour', 'trilinear', 'tricubic',
                               'triquintic', 'wsinc5'}, optional
            The interpolation method used for the transformed file.

        Retruns
        -------
        transformed_file: str
            Path to the transformed file, in template space.
        """
        self._check_anat_fitted()

        transformed_file = _apply_transforms(
            in_file,
            self.template,
            self.output_dir,
            self._normalization_transforms,
            transforms_kind=self.registration_kind,
            interpolation=interpolation,
            caching=self.caching,
            verbose=self.verbose)
        return transformed_file

    def _check_modality_fitted(self, modality):
        registered_modality = 'registered_{}_'.format(modality)
        if not hasattr(self, registered_modality):
            raise ValueError('It seems that {0} has not been {1} fitted. '
                             'You must call fit_modality() before calling '
                             'transform_modality().'.format(self.__class__.__name__, modality))

    def fit_modality(self, in_file, modality, slice_timing=True, t_r=None,
                     prior_rigid_body_registration=False, voxel_size=None):
        """ Estimates registration from the space of a given modality to
            the template space.

        Parameters
        ----------
        in_file : str
            Path to the modality image. M0 file is expected for perfusion.

        modality : one of {'func', 'perf'}
            Name of the modality.
        """
        self._check_anat_fitted()
        setattr(self, modality + '_', in_file)
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

        unbiased_file = ants_n4(to_coregister_file,
                                      write_dir=self.output_dir,
                                      terminal_output=self.terminal_output,
                                      caching=self.caching)

        if prior_rigid_body_registration:
            if self.use_rats_tool:
                compute_brain_mask = compute_morpho_brain_mask
            else:
                compute_brain_mask = compute_histo_brain_mask               

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
                caching=self.caching,
                verbose=self.verbose)
            self._func_undistort_warps = coregistration.coreg_warps_
            self.anat_in_func_space_ = coregistration.coreg_anat_
            self._func_to_anat_transform = coregistration.coreg_transform_
            self.undistorted_func_ = _apply_perslice_warp(
                allineated_file, self._func_undistort_warps, .1, .1,
                write_dir=self.output_dir, caching=self.caching)

            self.registered_func_ = _apply_transforms(
                self.undistorted_func_, self.template, self.output_dir,
                self._normalization_transforms + [self._func_to_anat_transform],
                transforms_kind=self.registration_kind,
                voxel_size=voxel_size, caching=self.caching)
        elif modality == 'perf':
            self.perf_brain_ = brain_file
            coregistration = coregister_perf(
                self._unifized_anat, unbiased_file,
                self.output_dir,
                anat_brain_file=self.anat_brain_,
                m0_brain_file=brain_file,
                prior_rigid_body_registration=prior_rigid_body_registration,
                caching=self.caching,
                verbose=self.verbose)
            self.undistorted_perf_ = coregistration.coreg_m0_
            self._perf_undistort_warps = coregistration.coreg_warps_
            self.anat_in_perf_space_ = coregistration.coreg_anat_
            self._perf_to_anat_transform = coregistration.coreg_transform_

            self.registered_perf_ = _apply_transforms(
                self.undistorted_perf_, self.template, self.output_dir,
                self._normalization_transforms + [self._perf_to_anat_transform],
                transforms_kind=self.registration_kind,
                voxel_size=voxel_size, caching=self.caching,
                verbose=self.verbose)

        return self

    def transform_modality_like(self, in_file, modality,
                                interpolation='wsinc5', voxel_size=None):
        """Transforms the given file from the space of the given modality to
            the template space. If the given modality has been corrected for
            EPI distorsions, the same correction is applied.

        Parameters
        ----------
        in_file : str
            Path to the file in the same space as the modality image.

        modality : one of {'func', 'perf'}
            Name of the modality.

        interpolation : one of {'nearestneighbour', 'trilinear', 'tricubic',
                               'triquintic', 'wsinc5'}, optional
            The interpolation method used for the transformed file.

        voxel_size : 3-tuple or None, optional
            The target voxels size. If None, the final voxels size will match
            the template.

        Retruns
        -------
        transformed_file : str
            Path to the transformed file, in template space.
        """
        self._check_anat_fitted()
        self._check_modality_fitted(modality)
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
            self._normalization_transforms + [coreg_transform_file],
            transforms_kind=self.registration_kind,
            voxel_size=voxel_size, caching=self.caching, verbose=self.verbose)
        return normalized_file

    def inverse_transform_towards_modality(self, in_file, modality,
                                           interpolation='wsinc5'):
        """ Trasnforms the given file from template space to modality space.

        Parameters
        ----------
        in_file : str
            Path to the file in the same space as the modality image.

        interpolation : one of {'nearestneighbour', 'trilinear', 'tricubic',
                               'triquintic', 'wsinc5'}, optional
            The interpolation method used for the transformed file.

        voxel_size : 3-tuple or None, optional
            The target voxels size. If None, the final voxels size will match
            the template.

        Retruns
        -------
        transformed_file : str
            Path to the transformed file, in template space.
        """
        self._check_anat_fitted()
        self._check_modality_fitted(modality)
        coreg_transform_file = self.__getattribute__(
            '_{}_to_anat_transform'.format(modality))
        modality_file = self.__getattribute__(modality + '_')

        inverted_file = _apply_transforms(in_file, modality_file,
                                          self.output_dir,
                                          self._normalization_transforms +\
                                          [coreg_transform_file],
                                          transforms_kind=self.registration_kind,
                                          inverse=True,
                                          interpolation=interpolation,
                                          caching=self.caching,
                                          verbose=self.verbose)
        return inverted_file
