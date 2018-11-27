"""
Preprocessing and analysis tools for small mammals brain MRI data in python
---------------------------------------------------------------------------

Documentation is available in the docstrings and online at
http://sammba-mri.github.io.

Contents
--------
sammba-MRI is a Python module for performant small animals MRI.

Submodules
---------
data_fetchers           --- Utilities to download small mammals brain MRI
                            datasets
io_conversions          --- Convert acquisition/export file formats to simpler
                            ones (usually NIfTI-1) for further processing
registration            --- AFNI-based pipelines to perform registration using
                            nipype interfaces
modality_processors     --- Functions for processing raw data of various MRI
                            modalities that are not BOLD fMRI, such as perfusion
interfaces              --- nipype-like interfaces for small animals tools
"""
import sys  

reload(sys)  
sys.setdefaultencoding('utf8')
import gzip

from .version import _check_module_dependencies, __version__

_check_module_dependencies()

__all__ = ['__version__', 'data_fetchers', 'io_conversions', 'segmentation',
           'preprocessing', 'registration', 'modality_processors']