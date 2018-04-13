import os
import abc
from .base import (compute_brain_mask,
                   _afni_bias_correct, _apply_mask,
                   mask_report)
from sklearn.base import BaseEstimator, TransformerMixin


class BaseRegistrator(BaseEstimator, TransformerMixin):
    """
    Base class for Registrators

    """
    @abc.abstractmethod
    def _check_inputs(self):
        raise NotImplementedError()

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