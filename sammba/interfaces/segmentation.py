""" Interface for brain mask computation based on nilearn's functions

    Change directory to provide relative paths for doctests
    >>> import os
    >>> current_filepath = os.path.realpath(__file__)
    >>> sammba_dir = os.path.dirname(os.path.dirname(current_filepath))
    >>> data_dir = os.path.realpath(os.path.join(sammba_dir, 'testing_data'))
    >>> os.chdir(data_dir)
"""
import os
import copy
import numpy as np
import nibabel
from nilearn import image, masking
from scipy.ndimage.morphology import (generate_binary_structure,
                                      binary_closing,
                                      binary_fill_holes,
                                      grey_dilation)
from sammba.externals.nipype.interfaces.base import (TraitedSpec,
                                                     BaseInterfaceInputSpec,
                                                     BaseInterface,
                                                     CommandLine,
                                                     CommandLineInputSpec,
                                                     traits,
                                                     isdefined)
from sammba.externals.nipype.utils.filemanip import fname_presuffix
from sammba.externals.nipype import logging

LOGGER = logging.getLogger('interface')
VALID_TERMINAL_OUTPUT = ['stream', 'allatonce', 'file', 'file_split',
                         'file_stdout', 'file_stderr', 'none']


class Info(object):
    """Handle RATS version information.

    version refers to the version of RATS on the system

    """
    @staticmethod
    def version():
        """Check for RATS version on system

        Parameters
        ----------
        None

        Returns
        -------
        version : str
           Version number as string or None if RATS not found

        """
        try:
            clout = CommandLine(command='RATS_MM --version',
                                resource_monitor=False,
                                terminal_output='allatonce').run()
        except IOError:
            # If rats_vcheck is not present, return None
            LOGGER.warn('RATS_MM executable not found.')
            return None

        version_stamp = clout.runtime.stdout.split('\n')[1].split(
            'version ')[0]
        if version_stamp.startswith('RATS_MM'):
            version_stamp = version_stamp.split('version: ')[1]
        else:
            return None

        return version_stamp


class HistogramMaskInputSpec(BaseInterfaceInputSpec):
    in_file = traits.File(
        desc="Input Image",
        exists=True,
        mandatory=True)
    out_file = traits.File(
        name_template='%s_histo_mask',
        name_source='in_file',
        keep_extension=True,
        desc="Output Image")
    volume_threshold = traits.Int(
        1650,
        desc="Volume threshold. [default: 1650]",
        usedefault=True)
    intensity_threshold = traits.Int(
        desc="Intensity threshold. "
             "[default: 500]")
    terminal_output = traits.Enum(
        'stream', 'allatonce', 'file', 'none',
        deprecated='1.0.0',
        desc=('Control terminal output: `stream` - '
              'displays to terminal immediately (default), '
              '`allatonce` - waits till command is '
              'finished to display output, `file` - '
              'writes output to file, `none` - output'
              ' is ignored'),
        nohash=True)
    lower_cutoff = traits.Float(
        .2,
        desc="lower fraction of the histogram to be discarded. In case of "
             "failure, it is usually advisable to increase lower_cutoff "
             "[default: 0.2]",
        usedefault=True)
    upper_cutoff = traits.Float(
        .85,
        desc="upper fraction of the histogram to be discarded."
             "[default: 0.85]",
        usedefault=True)
    connected = traits.Bool(
        True,
        desc="keep only the largest connect component",
        usedefault=True)
    opening = traits.Int(
        5,
        desc="Order of the morphological opening to perform, to keep only "
             "large structures. This step is useful to remove parts of the "
             "skull that might have been included. If the opening order is "
             "`n` > 0, 2`n` closing operations are performed after estimation "
             "of the largest connected constituent, followed by `n` erosions. "
             "This corresponds to 1 opening operation of order `n` followed "
             "by a closing operator of order `n`. Note that turning off "
             "opening (opening=0) will also prevent any smoothing applied to "
             "the image during the mask computation. [default: 2]",
        usedefault=True)
    closing = traits.Int(
        0,
        desc="Number of binary closing iterations to post-process the mask. "
             "[default: 10]",
        usedefault=True)
    dilation_size = traits.Tuple(
        (1, 1, 2),
        desc="Element size for binary dilation if needed",
        usedefault=True)
    verbose = traits.Bool(
        False,
        desc="be very verbose",
        usedefault=True)


class HistogramMaskOutputSpec(TraitedSpec):
    out_file = traits.File(desc="Brain mask file", exists=True)


class HistogramMask(BaseInterface):
    """Interface for nilearn.masking.compute_epi_mask, based on T. Nichols
    heuristics.

    Examples
    ========

    >>> from sammba.interfaces import HistogramMask
    >>> nichols_masker = HistogramMask()
    >>> nichols_masker.inputs.in_file = 'structural.nii'
    >>> nichols_masker.inputs.volume_threshold = 1650
    >>> res = nichols_masker.run()  # doctest: +SKIP
    """
    input_spec = HistogramMaskInputSpec
    output_spec = HistogramMaskOutputSpec
    _terminal_output = 'stream'

    @classmethod
    def set_default_terminal_output(cls, output_type):
        """Set the default terminal output for CommandLine Interfaces.

        This method is used to set default terminal output for
        CommandLine Interfaces.  However, setting this will not
        update the output type for any existing instances.  For these,
        assign the <instance>.terminal_output.
        """

        if output_type in VALID_TERMINAL_OUTPUT:
            cls._terminal_output = output_type
        else:
            raise AttributeError('Invalid terminal output_type: %s' %
                                 output_type)

    def __init__(self, terminal_output=None, **inputs):
        super(HistogramMask, self).__init__(**inputs)
        if terminal_output is not None:
            self.terminal_output = terminal_output

        # Attach terminal_output callback for backwards compatibility
        self.inputs.on_trait_change(self._terminal_output_update,
                                    'terminal_output')

    @property
    def terminal_output(self):
        return self._terminal_output

    @terminal_output.setter
    def terminal_output(self, value):
        if value not in VALID_TERMINAL_OUTPUT:
            raise RuntimeError(
                'Setting invalid value "%s" for terminal_output. Valid values '
                'are %s.' % (value, ', '.join(['"%s"' % v for v in
                                               VALID_TERMINAL_OUTPUT])))
        self._terminal_output = value

    def _terminal_output_update(self):
        self.terminal_output = self.terminal_output

    def _run_interface(self, runtime):
        if isdefined(self.inputs.intensity_threshold):
            threshold = self.inputs.intensity_threshold
            img = image.math_img('img>{0}'.format(threshold),
                                 img=self.inputs.in_file)
        else:
            img = nibabel.load(self.inputs.in_file)

        lower_cutoff = self.inputs.upper_cutoff - 0.05
        mask_img = masking.compute_epi_mask(
            img,
            lower_cutoff=lower_cutoff,
            upper_cutoff=self.inputs.upper_cutoff,
            connected=self.inputs.connected,
            opening=self.inputs.opening)
        mask_data = mask_img.get_data()
        n_voxels_mask = np.sum(mask_data > 0)

        # Find the optimal lower cutoff
        affine_det = np.abs(np.linalg.det(mask_img.affine[:3, :3]))
        n_voxels_min = int(self.inputs.volume_threshold * .9 / affine_det)
        while (n_voxels_mask < n_voxels_min) and (lower_cutoff >
                                                  self.inputs.lower_cutoff):
            lower_cutoff -= .05
            mask_img = masking.compute_epi_mask(
                img,
                lower_cutoff=lower_cutoff,
                upper_cutoff=self.inputs.upper_cutoff,
                connected=self.inputs.connected,
                opening=self.inputs.opening)
            mask_data = mask_img.get_data()
            n_voxels_mask = np.sum(mask_data > 0)
            if self.inputs.verbose:
                print('volume {0}, lower_cutoff {1}'.format(
                    n_voxels_mask * affine_det, lower_cutoff))

        n_voxels_max = int(self.inputs.volume_threshold * 1.1 / affine_det)
        previous_n_voxels_mask = copy.copy(lower_cutoff)
        while (n_voxels_mask >
               n_voxels_max) and (previous_n_voxels_mask >=
                                  n_voxels_mask) and (lower_cutoff + 0.01 <
                                                      self.inputs.upper_cutoff):
            lower_cutoff += .01
            mask_img = masking.compute_epi_mask(
                img,
                lower_cutoff=lower_cutoff,
                upper_cutoff=self.inputs.upper_cutoff,
                connected=self.inputs.connected,
                opening=self.inputs.opening)
            mask_data = mask_img.get_data()
            n_voxels_mask = np.sum(mask_data > 0)
            if self.inputs.verbose:
                print('volume {0}, lower_cutoff {1}'.format(
                    n_voxels_mask * affine_det, lower_cutoff))

        else:
            if n_voxels_mask < n_voxels_min:
                lower_cutoff -= .01
                mask_img = masking.compute_epi_mask(
                    img,
                    lower_cutoff=lower_cutoff,
                    upper_cutoff=self.inputs.upper_cutoff,
                    connected=self.inputs.connected,
                    opening=self.inputs.opening)
                mask_data = mask_img.get_data()
                n_voxels_mask = np.sum(mask_data > 0)
                if self.inputs.verbose:
                    print('volume {0}, lower_cutoff {1}'.format(
                        n_voxels_mask * affine_det, lower_cutoff))

        # Find the optimal opening
        n_voxels_max = int(self.inputs.volume_threshold * 1.5 / affine_det)
        opening = 0
        while n_voxels_mask > n_voxels_max and opening < self.inputs.opening:
            opening += 1
            mask_img = masking.compute_epi_mask(
                img,
                lower_cutoff=lower_cutoff,
                upper_cutoff=self.inputs.upper_cutoff,
                connected=self.inputs.connected,
                opening=opening)
            mask_data = mask_img.get_data()
            n_voxels_mask = np.sum(mask_data > 0)
            if self.inputs.verbose:
                print('volume {0}, lower_cutoff {1}, opening {2}'.format(
                    n_voxels_mask * affine_det, lower_cutoff, opening))

        # Find the optimal closing
        iterations = 0
        n_voxels_min = int(self.inputs.volume_threshold * .8 / affine_det)
        while (n_voxels_mask < n_voxels_min) and (iterations <
                                                  self.inputs.closing):
            iterations += 1
            for structure_size in range(1, 4):
                structure = generate_binary_structure(3, structure_size)
                mask_data = binary_closing(
                    mask_data, structure=structure, iterations=iterations)
                n_voxels_mask = np.sum(mask_data > 0)
                if self.inputs.verbose:
                    print('volume {0}, lower_cutoff {1}, opening {2}, closing '
                          '{3}'.format(n_voxels_mask * affine_det,
                                       lower_cutoff, opening, iterations))
                if n_voxels_mask > n_voxels_min:
                    break

        print('volume {0}, lower_cutoff {1}, opening {2}, closing '
              '{3}'.format(n_voxels_mask * affine_det, lower_cutoff,
                           opening, iterations))

        # Fill holes
        n_voxels_min = int(self.inputs.volume_threshold / affine_det)
        if n_voxels_mask < n_voxels_min:
            for structure_size in range(1, 4):
                structure = generate_binary_structure(3, structure_size)
                mask_data = binary_fill_holes(
                    mask_data, structure=structure)
                n_voxels_mask = np.sum(mask_data > 0)
                if self.inputs.verbose:
                    print('volume {0}, structure_size {1}'.format(
                        n_voxels_mask * affine_det, structure_size))
                if n_voxels_mask > n_voxels_min:
                    break

        # Dilation if needed
        size = self.inputs.dilation_size
        for n in range(3):
            if n_voxels_mask < n_voxels_min * .9:
                previous_n_voxels_mask = copy.copy(n_voxels_mask)
                mask_data = grey_dilation(mask_data, size=size)
                n_voxels_mask = np.sum(mask_data > 0)
                if n_voxels_mask == previous_n_voxels_mask:
                    size2 = (size[0] + 1, size[1] + 1, size[2] + 1)
                    mask_data = grey_dilation(mask_data, size=size2)
                    n_voxels_mask = np.sum(mask_data > 0)

                if self.inputs.verbose:
                    print('volume {0}, grey dilation {1}'.format(
                        n_voxels_mask * affine_det, n))
            else:
                break

        # Fill holes
        for structure_size in range(1, 4):
            structure = generate_binary_structure(3, structure_size)
            mask_data = binary_fill_holes(
                mask_data, structure=structure)
            n_voxels_mask = np.sum(mask_data > 0)

        final_volume = n_voxels_mask * affine_det
        if final_volume > 2 * self.inputs.volume_threshold:
            raise ValueError('Failed brain extraction: final volume is {0} '
                             ''.format(final_volume))
        elif final_volume < .5 * self.inputs.volume_threshold:
            raise ValueError('Failed brain extraction: final volume is {0} '
                             ''.format(final_volume))
        if self.inputs.verbose:
            print('final volume {0}'.format(final_volume))

        mask_img = image.new_img_like(mask_img, mask_data, mask_img.affine)
        if isdefined(self.inputs.out_file):
            mask_img.to_filename(os.path.abspath(self.inputs.out_file))
        else:
            mask_img.to_filename(os.path.abspath(
                fname_presuffix(os.path.basename(self.inputs.in_file),
                                suffix='_histo_mask')))
        return runtime

    def _list_outputs(self):
        outputs = self.output_spec().get()
        if isdefined(self.inputs.out_file):
            outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        else:
            outputs['out_file'] = os.path.abspath(
                fname_presuffix(os.path.basename(self.inputs.in_file),
                                suffix='_histo_mask'))
        return outputs


class MathMorphoMaskInputSpec(CommandLineInputSpec):
    in_file = traits.File(
        desc="Input Image",
        exists=True,
        mandatory=True,
        argstr="%s",
        position=0)
    out_file = traits.File(
        desc="Output Image",
        name_template='%s_morpho_mask',
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


class MathMorphoMaskOutputSpec(TraitedSpec):
    out_file = traits.File(desc="Brain mask file", exists=True)


class MathMorphoMask(CommandLine):
    """Mathemetical morphology stage for RATS, as described in:
    RATS: Rapid Automatic Tissue Segmentation in rodent brain MRI. Journal
    of neuroscience methods (2014) vol. 221 pp. 175 - 182.

    Author(s): Ipek Oguz, Milan Sonka

    Examples
    ========

    >>> from sammba.interfaces import MathMorphoMask
    >>> rats_masker = MathMorphoMask()
    >>> rats_masker.inputs.in_file = 'structural.nii'
    >>> rats_masker.inputs.intensity_threshold = 1000
    >>> rats_masker.cmdline  # doctest: +IGNORE_UNICODE
    'RATS_MM structural.nii structural_morpho_mask.nii -t 1000'
    >>> res = rats_masker.run()  # doctest: +SKIP
    """
    input_spec = MathMorphoMaskInputSpec
    output_spec = MathMorphoMaskOutputSpec
    _cmd = 'RATS_MM'

    def _list_outputs(self):
        outputs = self.output_spec().get()
        if isdefined(self.inputs.out_file):
            outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        else:
            outputs['out_file'] = os.path.abspath(
                self._filename_from_source('out_file'))

        mask_img = nibabel.load(outputs['out_file'])
        affine_det = np.abs(np.linalg.det(mask_img.affine[:3, :3]))
        n_voxels_mask = np.sum(mask_img.get_data() > 0)
        final_volume = n_voxels_mask * affine_det
        if final_volume > 2 * self.inputs.volume_threshold:
            raise ValueError('Failed brain extraction: final volume is {0} '
                             ''.format(final_volume))
        elif final_volume < .5 * self.inputs.volume_threshold:
            raise ValueError('Failed brain extraction: final volume is {0} '
                             ''.format(final_volume))

        return outputs
