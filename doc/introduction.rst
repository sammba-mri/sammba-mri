=====================================
Introduction: SmAll-maMMals BrAin MRI
=====================================

.. contents:: **Contents**
    :local:
    :depth: 1


What is sammba-MRI: small mammals neuroimaging with python
===========================================================

    sammba-MRI builds **easy-to-use** pipelines to process small mammals brain MRI multimodal images. Examples include (but are not restricted to) mouse models of Alzheimer disease.


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
* Pandas >= 0.12

as well as the interfaces:

* FSL >= 4.1.0
* SPM8/12

If you want to run the tests, you need nose >= 1.2.1

Interfaces configuration
========================
**Configuring FSL**: On an Ubuntu system, FSL is usually installed at :: /usr/share/fsl. You need to add this location to your .bashrc file. Edit this file by running the shell command::

    gedit ~/.bashrc

and add the following lines::

    # FSL
    FSLDIR=/usr/share/fsl
    . ${FSLDIR}/5.0/etc/fslconf/fsl.sh
    PATH=${FSLDIR}/5.0/bin:${PATH}
    export FSLDIR PATH

To test if FSL is correctly installed, open a new terminal and type in the shell command::

    fsl

You should see the FSL GUI with the version number in the header.

**Configuring SPM**: Add the following lines specifying the location of the spm folder to your .bashrc file::

    # SPM8
    export SPM_PATH=/i2bm/local/spm8-standalone/spm8_mcr/spm8

**Using SPM MCR**: If you don't have a matlab licence, specify the location of the Matlab Compiler Runtime and force the
use of the standalone MCR version of spm by appending the following lines to the .bashrc::

    # SPM MCR
    export SPMMCRCMD='/home/salma/Téléchargements/spm8/run_spm8.sh /home/salma/Téléchargements/MCR/v713 script'
    export FORCE_SPMMCR='True'

Installation
============
For the moment sammba-mri is available as a development version. To download the source code, run the shell command::

    git clone https://github.com/sammba-mri/sammba-mri

In the ``sammba-mri`` directory created by the previous step, run
(again, as a shell command)::

    python setup.py install --user

