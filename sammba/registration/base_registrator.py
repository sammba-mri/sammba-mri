import os
import abc
from ..segmentation.brain_mask import (compute_histo_brain_mask,
                                       compute_morpho_brain_mask,
                                       _apply_mask)
from ..preprocessing.bias_correction import afni_unifize
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

    def segment(self, in_file):
        """ Bias field correction and brain extraction
        """
        self._fit()
        unifized_file = afni_unifize(
            in_file, write_dir=self.output_dir,
            terminal_output=self.terminal_output, caching=self.caching)

        if self.use_rats_tool:
            compute_brain_mask = compute_morpho_brain_mask
        else:
            compute_brain_mask = compute_histo_brain_mask               

        if self.mask_clipping_fraction == .2:
            # do not repeat unifization step, as .2 is the default fraction
            brain_mask_file = compute_brain_mask(
                unifized_file, self.brain_volume, write_dir=self.output_dir,
                caching=self.caching,
                terminal_output=self.terminal_output,
                unifize=False)
        elif self.mask_clipping_fraction is None:
            brain_mask_file = compute_brain_mask(
                in_file, self.brain_volume, write_dir=self.output_dir,
                caching=self.caching,
                terminal_output=self.terminal_output,
                use_rats_tool=self.use_rats_tool,
                unifize=False)
        else:
            # custom unifization step, as .2 is the default fraction
            brain_mask_file = compute_brain_mask(
                in_file, self.brain_volume, write_dir=self.output_dir,
                caching=self.caching,
                terminal_output=self.terminal_output,
                cl_frac=self.mask_clipping_fraction)

        brain_file = _apply_mask(unifized_file, brain_mask_file,
                                 write_dir=self.output_dir,
                                 caching=self.caching,
                                 terminal_output=self.terminal_output)
        return unifized_file, brain_file

    def _fit(self, y=None):
        self._check_inputs()
        if self.verbose:
            self.terminal_output = 'stream'
        else:
            self.terminal_output = 'none'
        self._set_output_dir()
        return self