=====================================
Introduction: Sammba-MRI
=====================================

.. contents:: **Contents**
    :local:
    :depth: 1


What is sammba-MRI: small mammals neuroimaging with python
===========================================================

Sammba-MRI provides **easy-to-use** pipelines to process and analyze small mammals brain MRI multimodal images. 
Sammba-MRI will perform automatically several critical steps for MR image analysis.


° Conversion of Bruker DICOM filles to NIFTI-1

° Image quality check

	° Image registration and creation of a template

	° Transformation of individual dataset to the template (or to an atlas)

	° Evaluation of cerebral atrophy on the basis of an atlas

	° Estimation of cerebral perfusion maps from FAIR EPI images

	° Resting state fMRI analysis connectivity  and brain images visualization are straightforward with nilearn once the registration is performed.

Sammba-MRI integrates functionalities from a number of other packages (listed under the dependencies section below)
Examples include (but are not restricted to) mouse models of Alzheimer disease.


Dependencies
============
The required dependencies to use the package are the software:

* FSL >= 5.0
* AFNI
* ANTS
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


Installation (for Linux environment)
====================================

First install the scientific Python distribution `Anaconda <https://www.anaconda.com/distribution>_` and the version control system `Git <https://git-scm.com>_`. Then install nilearn as explained in `nilearn installation page <http://nilearn.github.io/introduction.html#installing-nilearn/>`_ and `nipype <https://nipype.readthedocs.io/en/latest/users/install.html>_`.

Sammba-mri is available as a development version. To download the source code, run the shell command::

    git clone https://github.com/sammba-mri/sammba-mri.git

In the ``sammba-mri`` directory created by the previous step, run
(again, as a shell command)::

    python setup.py install --user

To check your installation, open IPython by writing "ipython" in the terminal and pressing "Enter" and type in the following line and press "Enter":

In [1]: import sammba

If no error occurs, you have sammba nilearn correctly.


Interfaces configuration 
========================
**Configuring AFNI**: To be able to run AFNI make sure to add the following lines of code to your .bashrc file::

    # AFNI
    export PATH=/usr/lib/afni/bin:$PATH

**Configuring ANTS**: To be able to run ANTS make sure to add the following lines of code to your .bashrc file::

    # ANTs
    export PATH=/usr/local/antsbin/bin:$PATH
