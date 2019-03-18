"""
Unit tests for the context.data_extract.py file
"""

import unittest
import speedcom.tests.context

#Variables for testing
testing_dir = 'DATA_CLEAN_TEST_DIR'
cid_test = "ZZZ_64-17-5_ethanol.txt"
name_test = "ZZZ_ethanol_ZZZ.txt"

# Unit tests
class get_emission_files(unittest.TestCase):
    def test_nonexistant_directory(self):
        """
        Tests that a non-existant directory returns an error.  Should be
        handeled by built in python os package.
        """
        DATA_DIR = "NONEXISTANTDIRECTORY"
        self.assertRaises(Exception, lambda:context.data_extract.\
            get_emission_files())

    def test_emission_files(self):
        """
        Tests the proper number of ems.txt files are returned.
        """
        assert len(context.data_extract.get_emission_files(testing_dir)) == 3,\
            'get_emission_files gets improper number of files'

class get_absorption_files(unittest.TestCase):
    def test_nonexistant_directory(self):
        """
        Tests that a non-existant directory returns an error.  Should be 
        handeled by built in python os package.
        """
        self.assertRaises(Exception, lambda:context.data_extract.\
        get_absorption_files("NONXISTANTDIRECTORY"))

    def test_emission_files(self):
        """
        Tests the proper number of abs.txt files are returned.
        """
        assert len(context.data_extract.get_absorption_files(testing_dir)) == \
        3, 'get_emission_files gets improper number of files'

class get_molecule_cid(unittest.TestCase):
    def test_cid_from_cas(self):
        """
        Tests able to get CID from the cas number
        """
        assert context.data_extract.get_molecule_cid(cid_test) == 702, \
            'get_molecule_cid unable to ID from cid'
        return

    def test_cid_from_name(self):
        """
        Tests able to get CID from name
        """
        assert context.data_extract.get_molecule_cid(name_test) == 702, \
            'get_molecule_cid unable to ID from name'
        return

class get_molecular_weight(unittest.TestCase):
    def test_proper_weight(self):
        """
        Ensures able to get moleular weight from CID = 702
        """
        assert context.data_extract.get_molecular_weight(702) == 46.07, \
            'get_molecular_weight returns improper weight'
        return

class get_molecular_smiles(unittest.TestCase):
    def test_proper_smiles(self):
        """
        Ensures able to get isomeric SMILES from CID
        """
        assert context.data_extract.get_molecular_smiles(702) == 'CCO', \
            'get_molecular_smiles returns improper canonical SMILES'
        return

class get_spectra(unittest.TestCase):
    def test_proper_length(self):
        """
        Ensure returned spectrum is of the proper length
        """
        assert len(context.data_extract.get_spectra(testing_dir + \
        "1_118-96-7_X.abs.txt")) == 20, \
        "get_spectra returns spectum of improper size"
        return

class get_peaks(unittest.TestCase):
    def test_proper_length(self):
        """
        Ensure returned peak list is of appropriate size
        """
        spectra = context.data_extract.get_spectra(testing_dir + \
            "1_118-96-7_X.abs.txt")
        assert len(context.data_extract.get_peaks(spectra)) == 2, \
            'get_peaks retuns improper number of peaks'
        return

class initiate_molecules(unittest.TestCase):
    def test_proper_size(self):
        """
        Ensure the appropriate number of molecule objects are returned
        upon initialization.
        """
        assert len(context.data_extract.initiate_molecules(testing_dir)) == 3,\
            "initiate_molecules unable to find all molecules"
        return

class electgrostatic_potentials(unittest.TestCase):
    def test_valid_return_with_known(self):
        assert context.data_extract.electrostatic_potentials('acetic acid') \
            == 6.20, "electrostatic_potentils returns incorrect value"
        return

    def test_valid_return_unknown(self):
        assert context.data_extract.electrostatic_potentials("NOT MOLECULES")\
            == 1.0, "electrostatic_potentials unable to return for unknowns"
        return
