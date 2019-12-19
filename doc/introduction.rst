=====================================
Introduction: Sammba-MRI
=====================================

.. contents:: **Contents**
    :local:
    :depth: 1


What is sammba-MRI: small mammals neuroimaging with Python
===========================================================

Sammba-MRI provides **easy-to-use** pipelines to process and analyze small mammals brain MRI multimodal images. 
Sammba-MRI will perform automatically several critical steps for MRI image analysis.


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

as well as the Python packages:

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


Installation
============

Installing required neuroimaging software
-----------------------------------------
**FSL**: Follow the instructions
from `FSL official installation guide <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation>`_. To be able to run FSL, your system first needs to know where the software is installed at. On a Ubuntu system, this is usually under ``/usr/share/fsl``. Therefore, add the following code to your  ``.bashrc`` file. (To open and edit your .bashrc file on Ubuntu, use the following command: gedit  ``~/.bashrc``)::

    #FSL
    FSLDIR=/usr/share/fsl
    . ${FSLDIR}/5.0/etc/fslconf/fsl.sh
    PATH=${FSLDIR}/5.0/bin:${PATH}
    export FSLDIR PATH

**AFNI**: If you have access to `Neurodebian <http://neuro.debian.net>`_, simply install the `AFNI package <http://neuro.debian.net/pkgs/afni.html>`_ through Neurodebian. Otherwise, go to `AFNI installation page <https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/background_install/main_toc.html>`_. Then add the following lines of code to your  ``.bashrc`` file::

    # AFNI
    export PATH=/usr/lib/afni/bin:$PATH

**ANTs**: The recommended way to install ANTs is to build it from source for your own system. Go to the ANTs section in `Michael Notter's excellent tutorial <http://miykael.github.io/nipype-beginner-s-guide/installation.html>`_ and follow each step. After the installation is complete, add the following lines  to your  ``.bashrc`` file ::

    # ANTs
    export PATH=/usr/local/antsbin/bin:$PATH
    export ANTSPATH=/usr/local/antsbin/bin/


Installing sammba-MRI
---------------------
For the moment, sammba-MRI has been tested only on Linux environment.
First install the scientific Python distribution `Anaconda <https://www.anaconda.com/distribution>`_ and the version control system `Git <https://git-scm.com>`_. Then install the Python modules `nipype <https://nipype.readthedocs.io/en/latest/users/install.html>`_ and `nilearn <http://nilearn.github.io/introduction.html#installing-nilearn/>`_.

Sammba-mri is available as a development version. To download the source code, run the shell command::

    git clone https://github.com/sammba-mri/sammba-mri.git

This will create a directory ``sammba-mri``. You now need to change to this directory and install the package by running (again, as a shell command)::

    cd sammba-mri
    python setup.py install --user

To check your installation, open IPython by writing "ipython" in the terminal and pressing "Enter" and type in the following line and press "Enter"::

    In [1]: import sammba

If no error occurs, you have installed sammba correctly.

