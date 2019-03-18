

class initiate_molecules(unittest.TestCase):
    def test_proper_size(self):
        """
        Ensure the appropriate number of molecule objects are returned
        upon initialization.
        """
        assert len(context.speedcom.initiate_molecules(testing_dir)) == 3,\
            "initiate_molecules unable to find all molecules"
        return


