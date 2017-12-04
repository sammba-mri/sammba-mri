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
The required dependencies to use the package are 

the softwares:

* AFNI
* FSL >= 5.0
* `RATS <http://www.iibi.uiowa.edu/content/rats-overview/>`_ for brain
  extraction

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


Installation
============

First install Anaconda and nilearn as explained in `nilearn installation page <http://nilearn.github.io/introduction.html#installing-nilearn/>`_.

Next install the remaining python dependencies using conda::

    conda install pandas configparser future traits simplejson networkx packaging funcsigs

Some dependencies are only available through pip::

    pip install packaging prov

Finally, you need `graphviz <http://www.graphviz.org/>`_::

    sudo apt-get install graphviz

For the moment sammba-mri is available as a development version. To download the source code, run the shell command::

    git clone https://github.com/sammba-mri/sammba-mri.git

In the ``sammba-mri`` directory created by the previous step, run
(again, as a shell command)::

    python setup.py install --user


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

**Configuring AFNI**: To be able to run AFNI make sure to add the following lines of code to your .bashrc file::

    # AFNI
    export PATH=/usr/lib/afni/bin:$PATH
