import cirpy
import numpy as np
import os
#import pubchempy as pcp
import scipy.signal

import ../core


"""
The data cleaning functions.  Will take in a list of data from the path and
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
  return [fn for fn in os.listdir(DATA_DIR)
              if any(fn.endswith(.ems.txt))]

def get_absorption_files():
  """
  Get the total files containing absorption spectra.  Expects the files to be
  of type *.txt and have the naming convention *.abs.txt will return a list of
  strings to be parsed later.
  """
  return [fn for fn in os.listdir(DATA_DIR)
              if any(fn.endswith(.abs.txt))]

def get_molecule_cas(file_name):
  """
  Take a file_name, remove the ending and prefix and return the CAS number.
  If the file_name does not have the cas number, try the name of the system,
  else return None so as to disreguard the molecule from training.
  """
  file = file_name.split("_")
  #try to get the cas # if not then try name.
  try:
    cas = file[1].split('-')
    cas = int(cas[0] + cas[1] + cas[2])
    cid = 

  #file doesn't have a cas number with it, try name
  except:
    pass

  try:
    cid = pcp.get_compounds[file[1], 'name'][0].cid

  #Couldn't determine name.
  except:
    cid = None

  return cas

def get_molecular_weight(cid):
  properties = []
  molecule = pcp.Compound.from_cid(cid)
  return properties

def get_spectra(file_name):
  """
  Iterates through the file and pulls out the normalized spectrum contained
  within.  The X is in wavelentgh, and the intensity is normalized to the
  most intense peak.
  """
  spectra = []
  with open (file_name,'r') as file:
    next(file)
    for line in file:
      tmp = line.split()
      spectra.append([float(tmp[0]), float(tmp[1])])
  return np.asarray(spectra)

def get_peaks(spectra):
  peaks = []
  """
  This will take in the spectra from the get_spectra function and will return
  a list of [x,Y] coordinates for the peaks.
  """
  peak_locations = scipy.signal.find_peaks(spectra[:,1])[0]
  for i in peak_locations:
    peaks.append(spectra[i, :])
  return np.asarry(peaks)
