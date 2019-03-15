# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
_version_major = 0
_version_minor = 1
_version_micro = ''  # use '' for first of series, number for 1 and above
_version_extra = 'dev'
# _version_extra = ''  # Uncomment this for full releases

# Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro:
    _ver.append(_version_micro)
if _version_extra:
    _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

CLASSIFIERS = ["Development Status :: 3 - Alpha",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: MIT License",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

# Description should be a one-liner:
description = "SPEEDCOM: A program to predict absorption and emission \
                spectra of small organic molecules"
# Long description will go up on the pypi page
long_description = """
SPEEDCOM (Spectra Prediction for the Excitation and Emission of Dyes and
Other Conjugated Organic Molecules)
========
SPEEDCOM is an open-source python package that predicts the fluorescence
excitation and emission spectra of small conjugated organic molecules.

Users can input the SMILES sting for a molecule via a GUI and be returned
the molecule's spectra, as well as some physical characteristics, such as
the molar extinction coefficient and quantum yield.

SPEEDCOM uses a deep-learning model trained on data from the PhotoChemCAD
database. Upon external use, the model is be pre-trained so can quickly
predict the spectra of an un-seen molecule.
License
=======
``SPEEDCOM`` is licensed under the MIT license. See the "LICENSE" file for
information on the history of this software, terms & conditions for usage,
and a DISCLAIMER OF ALL WARRANTIES.
All trademarks referenced herein are property of their respective holders.
Copyright (c) 2018 -- Emissible (Joe Abbott, Ryan Beck, Hang Hu, Yang Lui, Lixin
Li) at the The University of Washington.
"""

NAME = "speedcom"
MAINTAINER = "Joe Abbott"
MAINTAINER_EMAIL = "jwa7@uw.edu"
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = "https://github.com/emissible/SPEEDCOM"
DOWNLOAD_URL = ""
LICENSE = "MIT"
AUTHOR = "Joe Abbott, Ryan Beck, Hang Hu, Yang Lui, Lixin Li"
AUTHOR_EMAIL = "jwa7@uw.edu, rbeck4@uw.esu, hanghu@uw.edu, yliu92uw@uw.edu, \
                lixinlu@uw.edu"
PLATFORMS = "OS Independent"
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
VERSION = __version__
PACKAGES = ["speedcom"]
# PACKAGE_DATA = {'nets': ['speedcom/nets/*.h5']}
REQUIRES = []
