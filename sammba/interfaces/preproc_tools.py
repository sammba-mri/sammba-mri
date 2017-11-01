from nipype.interfaces.base import (TraitedSpec, CommandLineInputSpec,    
                                    CommandLine, traits)
import os


class RatsMMInputSpec(CommandLineInputSpec):
    in_file = traits.File(
        desc="Input Image",
        exists=True,
        mandatory=True,
        argstr="%s",
        position=0)
    out_file = traits.File(
        desc="Output Image",
        mandatory=True,
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
    >>> from sammba.data_fetchers import fetch_zurich_test_retest
    >>> data = fetch_zurich_test_retest(subjects=[0], verbose=0)
    >>> rats_masker = RatsMM()
    >>> rats_masker.inputs.in_file = data.anat[0]
    >>> rats_masker.inputs.intensity_threshold = 1000
    >>> rats_masker.inputs.out_file = 'masked.nii'
    >>> rats_masker.cmdline  # doctest: +ALLOW_UNICODE
    'RATS_MM /home/travis/nilearn_data/zurich_retest/baseline/1366/3DRARE.nii.gz masked.nii -t 1000'
    >>> res = rats_masker.run()  # doctest: +SKIP
    """
    input_spec = RatsMMInputSpec
    output_spec = RatsMMOutputSpec
    _cmd = 'RATS_MM'

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        return outputs
