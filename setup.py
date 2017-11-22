#! /usr/bin/env python

descr = """A set of Python modules for functional MRI..."""

import sys
import os

from setuptools import setup, find_packages


def load_version():
    """Executes sammba/version.py in a globals dictionary and return it.
    Note: importing sammba is not an option because there may be
    dependencies like nibabel which are not installed and
    setup.py is supposed to install them.
    """
    # load all vars into globals, otherwise
    #   the later function call using global vars doesn't work.
    globals_dict = {}
    with open(os.path.join('sammba', 'version.py')) as fp:
        exec(fp.read(), globals_dict)

    return globals_dict


def is_installing():
    # Allow command-lines such as "python setup.py build install"
    install_commands = set(['install', 'develop'])
    return install_commands.intersection(set(sys.argv))

# Make sources available using relative paths from this file's directory.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

_VERSION_GLOBALS = load_version()
DISTNAME = 'sammba-mri'
DESCRIPTION = 'preprocessing and analysis tools for small mammals brain MRI data in Python'
LONG_DESCRIPTION = open('README.rst').read()
MAINTAINER = 'Nachiket Nadkarni and Salma Bougacha'
MAINTAINER_EMAIL = 'salmabougacha@hotmail.com'
URL = 'http://sammba-mri.github.io'
LICENSE = 'CeCILL-B'
DOWNLOAD_URL = 'http://sammba-mri.github.io'
VERSION = _VERSION_GLOBALS['__version__']


if __name__ == "__main__":
    if is_installing():
        module_check_fn = _VERSION_GLOBALS['_check_module_dependencies']
        module_check_fn(is_sammba_installing=True)

    install_requires = \
        ['%s>=%s' % (mod, meta['min_version'])
            for mod, meta in _VERSION_GLOBALS['REQUIRED_MODULE_METADATA']
            if not meta['required_at_installation']]

    setup(name=DISTNAME,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          license=LICENSE,
          url=URL,
          version=VERSION,
          download_url=DOWNLOAD_URL,
          long_description=LONG_DESCRIPTION,
          zip_safe=False,  # the package can run out of an .egg file
          classifiers=[
              'Intended Audience :: Science/Research',
              'Intended Audience :: Developers',
              'License :: OSI Approved',
              'Programming Language :: C',
              'Programming Language :: Python',
              'Topic :: Software Development',
              'Topic :: Scientific/Engineering',
              'Operating System :: Microsoft :: Windows',
              'Operating System :: POSIX',
              'Operating System :: Unix',
              'Operating System :: MacOS',
              'Programming Language :: Python :: 2',
              'Programming Language :: Python :: 2.7',
              'Programming Language :: Python :: 3.3',
              'Programming Language :: Python :: 3.4',
          ],
          packages=find_packages(),
          package_data={'sammba.tests': ['*.nii.gz', '*.npz'],
                        #'sammba.description': ['*.rst'],
                        },
          install_requires=install_requires,)
