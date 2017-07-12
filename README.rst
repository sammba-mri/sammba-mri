.. -*- mode: rst -*-

.. image:: https://travis-ci.org/small-animal-MRI/sammba-MRI.svg?branch=master
   :target: https://travis-ci.org/small-animal-MRI/sammba-MRI
   :alt: Build Status

.. image:: https://coveralls.io/repos/github/small-animal-MRI/sammba-MRI/badge.svg?branch=master
    :target: https://coveralls.io/github/small-animal-MRI/sammba-MRI?branch=master

sammba-MRI
=====

Mindt builds relevant **pipelines** for processing multimodal small animal brain MRI images.

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

Installation
============

For the moment process-asl is available as a development version. You can download the source code with the command::

    git clone https://github.com/small-animal-MRI/sammba-MRI
