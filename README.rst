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

The required dependencies to use the software are the python packages:

* Python 2.7
* setuptools
* Numpy >= 1.6.2
* SciPy >= 0.11
* Nibabel >= 2.0.1
* Nilearn >= 0.1.3
* Matplotlib >= 1.1.1
* Nipype 0.11.0
* NetworkX >= 1.7
* Enthought Traits >= 4.3
* Dateutil >= 1.5

as well as the interfaces:

* FSL >= 4.1.0
* AFNI

User guide and gallery of examples are available on
====================================================
https://sammba-mri.github.io/


Installation
============

Sammba-MRI source code can be downloaded with the command::

    git clone https://github.com/sammba-mri/sammba-mri
