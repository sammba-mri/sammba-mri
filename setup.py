import os
from setuptools import setup


# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="sammba",
    version="alpha",
    maintainer="Nachiket Nadkarni",
    maintainer_email="nadkarni@fastmail.fm",
    description=("MRI Rodents processing in python."),
    license="CeCILL-B",
    keywords="rodents registration dicom",
    url="https://sammba-mri.github.io",
    packages=['sammba', ],
    long_description=read('README.rst'),
    classifiers=[
        "Topic :: Scientific/Engineering",
        "License :: CeCILL-B",
    ],
)
