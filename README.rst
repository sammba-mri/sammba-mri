.. -*- mode: rst -*-

.. image:: https://travis-ci.org/sammba-mri/sammba-mri.svg?branch=master
    :target: https://travis-ci.org/sammba-mri/sammba-mri

.. image:: https://coveralls.io/repos/github/sammba-mri/sammba-mri/badge.svg?branch=master
    :target: https://coveralls.io/github/sammba-mri/sammba-mri?branch=master

.. image:: https://circleci.com/gh/sammba-mri/sammba-mri.svg?style=svg
    :target: https://circleci.com/gh/sammba-mri/sammba-mri

sammba-MRI: small mammals neuroimaging with python
==========

Sammba-MRI provides easy-to-use **pipelines** to process and analyze small mammals brain MRI multimodal images. 
Sammba-MRI will perform automatically several critical steps for MR image analysis.

* Conversion of Bruker DICOM files to NIFTI-1
* Image quality check
* Image registration and creation of a template
* Transformation of individual dataset to the template (or to an atlas)
* Evaluation of cerebral atrophy on the basis of an atlas
* Estimation of cerebral perfusion maps from FAIR EPI images
* Resting state fMRI analysis connectivity  and brain images visualization are straightforward with nilearn once the registration is performed.

Sammba-MRI integrates functionalities from a number of other packages (listed under the dependencies section below)


Dependencies
============

The required dependencies to use the software are the software:

* FSL >= 5.0
* AFNI
* ANTs
* Python >= 3.5
* `RATS <http://www.iibi.uiowa.edu/content/rats-overview/>`_ for brain extraction

as well as the python packages:

* setuptools
* Numpy >= 1.14
* SciPy >= 0.19
* Nibabel >= 2.0.2
* Nilearn >= 0.4.0
* Sklearn >= 0.19
* Nipype >= 1.0.4

If you are running the examples, matplotlib >= 1.5.1 is required.

If you want to run the tests, you need nose >= 1.2.1 and doctest-ignore-unicode.

If you want to convert DICOM files to NIFTI files, you will need the
`DICOM ToolKit (DCMTK) <http://support.dcmtk.org/docs/index.html>`_ package


User guide and gallery of examples are available on
====================================================
https://sammba-mri.github.io/


Installation
============

Sammba-MRI source code can be downloaded with the command::

    git clone https://github.com/sammba-mri/sammba-mri
