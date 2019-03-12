"""
Unit tests for the data_extract.py file
"""

import unittest

from SPEEDCOM.data_clean import data_clean 

#Variables for testing
testing_dir = 'DATA_CLEAN_TEST_DIR'
cid_test = "ZZZ_64-17-5_ethanol.txt"
name_test = "ZZZ_ethanol_ZZZ.txt"

# Unit tests
class get_emission_files(unittest.TestCase):
  def test_nonexistant_directory(self):
    """
    Tests that a non-existant directory returns an error.  Should be handeld
    by built in python os package.
    """
    DATA_DIR = "NONEXISTANTDIRECTORY"
    self.assertRaises(Exception, lambda:data_extract.get_emission_files())

  def test_emission_files(self):
    """
    Tests the proper number of ems.txt files are returned.
    """
    DATA_DIR = testing_dir
    assert len(data_extract.get_emission_files()) == 3, \
      'get_emission_files gets improper number of files'

class get_absorption_files(unittest.TestCase):
  def test_nonexistant_directory(self):
    """
    Tests that a non-existant directory returns an error.  Should be handeld
    by built in python os package.
    """
    DATA_DIR = "NONEXISTANTDIRECTORY"
    self.assertRaises(Exception, lambda:data_extract.get_absorption_files())

  def test_emission_files(self):
    """
    Tests the proper number of abs.txt files are returned.
    """
    DATA_DIR = testing_dir
    assert len(data_extract.get_absorption_files()) == 3, \
      'get_emission_files gets improper number of files'

class get_molecule_cid(unittest.TestCase):
  def test_cid_from_cas(self):
    """
    Tests able to get CID from the cas number
    """
    assert data_extract.get_molecule_cid(cid_test) == 702, \
      'get_molecule_cid unable to ID from cid'
    return

  def test_cid_from_name(self):
    """
    Tests able to get CID from name
    """
    assert data_extract.get_molecule_cid(name_test) == 702, \
      'get_molecule_cid unable to ID from name'
    return

class get_molecular_weight(unittest.TestCase):
  def test_proper_weight(self):
    """
    Ensures able to get moleular weight from CID = 702
    """
    assert data_extract.get_molecular_weight(702) == 46.07, \
      'get_molecular_weight returns improper weight'
    return

class get_molecular_smiles(unittest.TestCase):
  def test_proper_smiles(self):
    """
    Ensures able to get isomeric SMILES from CID
    """
    assert data_extract.get_molecular_smiles(702) == 'CCO', \
      'get_molecular_smiles returns improper canonical SMILES'
    return

class get_spectra(unittest.TestCase):
  def test_proper_length(self):
    """
    Ensure returned spectrum is of the proper length
    """
    assert len(data_extract.get_spectra(testing_dir + "1_118-96-7_X.abs.txt"))\ 
      == 20, "get_spectra returns spectum of improper size"
    return

class get_peaks(unittest.TestCase):
  def test_proper_length(self):
    """
    Ensure returned peak list is of appropriate size
    """
    spectra = data_extract.get_spectra(testing_dir + "1_118-96-7_X.abs.txt")
    assert len(data_extract.get_peaks(spectra) == 2, \
      'get_peaks retuns improper number of peaks'
    return

class initiate_molecules(unittest.TestCase):
  def test_proper_size(self):
    """
    Ensure the appropriate number of molecule objects are returned
    upon initialization.
    """
    DATA_DIR = testing_dir
    assert len(data_extract.initiate_molecules()) == 3, \
      "initiate_molecules unable to find all molecules"
    return
