"""
Unit tests for the speedcom.py file
"""

import unittest
import speedcom.tests.context as context

class initiate_molecules(context.speedcom.TestCase):
    def test_proper_size(self):
        """
        Ensure the appropriate number of molecule objects are returned
        upon initialization.
        """
        assert len(speedcom.initiate_molecules(testing_dir)) == 3,\
            "initiate_molecules unable to find all molecules"
        return


