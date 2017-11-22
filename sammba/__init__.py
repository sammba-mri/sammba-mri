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
registration            --- AFNI-based pipelines to perform fMRI registration,
                            using nipype interfaces
interfaces              --- nipype-like interfaces for small animals tools
"""

import gzip

from .version import _check_module_dependencies, __version__

_check_module_dependencies()

__all__ = ['__version__', 'data_fetchers', 'registration', 'interfaces']