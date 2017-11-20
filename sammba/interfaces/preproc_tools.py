"""Useful preprocessing interfaces

    Change directory to provide relative paths for doctests
    >>> import os
    >>> current_filepath = os.path.realpath(__file__)
    >>> sammba_dir = os.path.dirname(os.path.dirname(current_filepath))
    >>> data_dir = os.path.realpath(os.path.join(sammba_dir, 'testing_data'))
    >>> os.chdir(data_dir)
"""
import os
from sammba.externals.nipype.interfaces.base import (TraitedSpec,
                                                     CommandLineInputSpec,
                                                     CommandLine, traits,
                                                     isdefined)


class RatsMMInputSpec(CommandLineInputSpec):
    in_file = traits.File(
        desc="Input Image",
        exists=True,
        mandatory=True,
        argstr="%s",
        position=0)
    out_file = traits.File(
        desc="Output Image",
        name_template='%s_masked',
        name_source='in_file',
        keep_extension=True,
        argstr="%s",
        position=1)
    volume_threshold = traits.Int(
        desc="Volume threshold (the parameter V in the paper). "
             "[default: 1650]",
        argstr="-v %s")
    intensity_threshold = traits.Int(
        desc="Intensity threshold (the parameter T in the paper). "
             "[default: 500]",
        argstr="-t %s")


class RatsMMOutputSpec(TraitedSpec):
    out_file = traits.File(desc="Brain masked file", exists=True)


class RatsMM(CommandLine):
    """Mathemetical morphology stage for RATS, as described in:
    RATS: Rapid Automatic Tissue Segmentation in rodent brain MRI. Journal
    of neuroscience methods (2014) vol. 221 pp. 175 - 182.

    Author(s): Ipek Oguz, Milan Sonka

    Examples
    ========

    >>> from sammba.interfaces import RatsMM
    >>> rats_masker = RatsMM()
    >>> rats_masker.inputs.in_file = 'structural.nii'
    >>> rats_masker.inputs.intensity_threshold = 1000
    >>> rats_masker.cmdline  # doctest: +IGNORE_UNICODE
    'RATS_MM structural.nii structural_masked.nii -t 1000'
    >>> res = rats_masker.run()  # doctest: +SKIP
    """
    input_spec = RatsMMInputSpec
    output_spec = RatsMMOutputSpec
    _cmd = 'RATS_MM'

    def _list_outputs(self):
        outputs = self.output_spec().get()
        if isdefined(self.inputs.out_file):
            outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        else:
            outputs['out_file'] = os.path.abspath(
                self._filename_from_source('out_file'))
        return outputs
