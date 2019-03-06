import scipy.signal

import ../core


"""
The data cleaning function.  Will take in a list of data from the path and
will output a numpy array containing the SMILES string with the spectral
features.
"""

#DATA_DIR global variable containing the path to the data
def get_emission_files():
  """
  Get the total files containing emission spectra.  Expects the files to be of
  type *.txt and have the naming convention *.ems.txt will return a list of
  strings to be parsed later.
  """
  files = []
  return files

def get_absorption_files():
  """
  Get the total files containing absorption spectra.  Expects the files to be
  of type *.txt and have the naming convention *.abs.txt will return a list of
  strings to be parsed later.
  """
  files = []
  return files

def get_molucule_name(file_name):
  
  return molucule_name

def get_spectra(file_name):
  spectra = []

  return spectra

def get_peaks(spectra):
  peaks = scipy.signal.find_peaks(spectra)
  return peaks
