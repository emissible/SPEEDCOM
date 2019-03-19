"""
Unit tests for the dataUtils.py file
"""

import speedcom.tests.context as context
from context.dataUtils import DataUtils

test_filename = 'random_files'
if os.path.exists(test_filename):
    os.romove




def test_readData():
    """
    Test function for readData
    """
    # test raise exception for random names

    try:
        DataUtils.readData('ARandomFile')
    except Exception as e:
        assert isinstance(e, AssertionError)
