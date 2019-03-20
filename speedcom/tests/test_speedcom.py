"""
Unit tests for the speedcom.py file
"""

import unittest
import speedcom.tests.context as context

testing_dir = './speedcom/tests/DATA_CLEAN_TEST_DIR/'

class initiate_molecules(unittest.TestCase):
    def test_proper_size(self):
        """
        Ensure the appropriate number of molecule objects are returned
        upon initialization.
        """
        mol_list_len = len(context.core.initiate_molecules(testing_dir, \
            test=True))
        assert mol_list_len == 5, "initiate_molecules unable to find all \
            molecules, len = %f" %(mol_list_len)
        return


