.. AUTO-GENERATED FILE -- DO NOT EDIT!

segmentation.interfaces
=======================


.. _sammba.segmentation.interfaces.HistogramMask:


.. index:: HistogramMask

HistogramMask
-------------

`Link to code <None>`__

Interface for nilearn.masking.compute_epi_mask, based on T. Nichols
heuristics.

Examples
~~~~~~~~

>>> from sammba.segmentation import HistogramMask
>>> nichols_masker = HistogramMask()
>>> nichols_masker.inputs.in_file = 'structural.nii'
>>> nichols_masker.inputs.volume_threshold = 1650
>>> res = nichols_masker.run()  # doctest: +SKIP

Inputs::

        [Mandatory]
        in_file: (a file name)
                Input Image

        [Optional]
        closing: (an integer (int or long), nipype default value: 0)
                Number of binary closing iterations to post-process the mask.
                [default: 10]
        connected: (a boolean, nipype default value: True)
                keep only the largest connect component
        dilation_size: (a tuple of the form: (a value of type 'int', a value
                 of type 'int', a value of type 'int'), nipype default value: (1, 1,
                 2))
                Element size for binary dilation if needed
        ignore_exception: (a boolean, nipype default value: False)
                Print an error message instead of throwing an exception in case the
                interface fails to run
        intensity_threshold: (an integer (int or long))
                Intensity threshold. [default: 500]
        lower_cutoff: (a float, nipype default value: 0.2)
                lower fraction of the histogram to be discarded. In case of failure,
                it is usually advisable to increase lower_cutoff [default: 0.2]
        opening: (an integer (int or long), nipype default value: 5)
                Order of the morphological opening to perform, to keep only large
                structures. This step is useful to remove parts of the skull that
                might have been included. If the opening order is `n` > 0, 2`n`
                closing operations are performed after estimation of the largest
                connected constituent, followed by `n` erosions. This corresponds to
                1 opening operation of order `n` followed by a closing operator of
                order `n`. Note that turning off opening (opening=0) will also
                prevent any smoothing applied to the image during the mask
                computation. [default: 2]
        out_file: (a file name)
                Output Image
        upper_cutoff: (a float, nipype default value: 0.85)
                upper fraction of the histogram to be discarded.[default: 0.85]
        verbose: (a boolean, nipype default value: False)
                be very verbose
        volume_threshold: (an integer (int or long), nipype default value:
                 1650)
                Volume threshold. [default: 1650]

Outputs::

        out_file: (a file name)
                Brain mask file

.. _sammba.segmentation.interfaces.MathMorphoMask:


.. index:: MathMorphoMask

MathMorphoMask
--------------

`Link to code <None>`__

Wraps command **RATS_MM**

Mathemetical morphology stage for RATS, as described in:
RATS: Rapid Automatic Tissue Segmentation in rodent brain MRI. Journal
of neuroscience methods (2014) vol. 221 pp. 175 - 182.

Author(s): Ipek Oguz, Milan Sonka

Examples
~~~~~~~~

>>> from sammba.segmentation import MathMorphoMask
>>> rats_masker = MathMorphoMask()
>>> rats_masker.inputs.in_file = 'structural.nii'
>>> rats_masker.inputs.intensity_threshold = 1000
>>> rats_masker.cmdline  # doctest: +IGNORE_UNICODE
'RATS_MM structural.nii structural_morpho_mask.nii -t 1000'
>>> res = rats_masker.run()  # doctest: +SKIP

Inputs::

        [Mandatory]
        in_file: (a file name)
                Input Image
                flag: %s, position: 0

        [Optional]
        args: (a unicode string)
                Additional parameters to the command
                flag: %s
        environ: (a dictionary with keys which are a newbytes or None or a
                 newstr or None and with values which are a newbytes or None or a
                 newstr or None, nipype default value: {})
                Environment variables
        ignore_exception: (a boolean, nipype default value: False)
                Print an error message instead of throwing an exception in case the
                interface fails to run
        intensity_threshold: (an integer (int or long))
                Intensity threshold (the parameter T in the paper). [default: 500]
                flag: -t %s
        out_file: (a file name)
                Output Image
                flag: %s, position: 1
        volume_threshold: (an integer (int or long))
                Volume threshold (the parameter V in the paper). [default: 1650]
                flag: -v %s

Outputs::

        out_file: (a file name)
                Brain mask file
