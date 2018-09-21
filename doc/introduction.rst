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
The required dependencies to use the package are 

the softwares:

* AFNI
* ANTS
* `RATS <http://www.iibi.uiowa.edu/content/rats-overview/>`_ for brain extraction

as well as the python packages:

* setuptools
* Numpy >= 1.6.2
* SciPy >= 0.11
* Nibabel >= 2.0.1
* Nilearn >= 0.1.3
* Matplotlib >= 1.1.1
* NetworkX >= 1.7
* Pandas >= 0.12
* packaging
* configparser
* future >= 0.16.0
* traits >= 4.6
* simplejson >= 3.8.0
* prov == 1.5.0
* funcsigs
* click

If you want to run the tests, you need nose >= 1.2.1 and doctest-ignore-unicode

If you want to convert DICOM files to NIFTI files, you will need the
`DICOM ToolKit (DCMTK) <http://support.dcmtk.org/docs/index.html>`_ package


Installation (for Linux environment)
============

First install Anaconda and nilearn as explained in `nilearn installation page <http://nilearn.github.io/introduction.html#installing-nilearn/>`_.

Next install the remaining python dependencies using conda::

    conda install pandas configparser future traits simplejson networkx funcsigs

Some dependencies are only available through pip::

    pip install packaging prov patsy

Finally, you need `graphviz <http://www.graphviz.org/>`_::

    sudo apt-get install graphviz

Sammba-mri is available as a development version. To download the source code, run the shell command::

    git clone https://github.com/sammba-mri/sammba-mri.git

In the ``sammba-mri`` directory created by the previous step, run
(again, as a shell command)::

    python setup.py install --user


Interfaces configuration 
========================
**Configuring AFNI**: To be able to run AFNI make sure to add the following lines of code to your .bashrc file::

    # AFNI
    export PATH=/usr/lib/afni/bin:$PATH

**Configuring AFNI**: To be able to run ANTS make sure to add the following lines of code to your .bashrc file::

    #ANTs
    export PATH=/usr/local/antsbin/bin:$PATH
