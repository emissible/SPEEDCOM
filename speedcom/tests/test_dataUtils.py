"""
Unit tests for the dataUtils.py file
"""
import json
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

        
def test_load_wordmap_from_json():
    """
    test function for load word map
    """
    with open('temp.json', 'w') as fp:
        json.dump({1:'#',2:'!'})
    try:
        word_map = DataUtils.load_wordmap_from_json('temp.json')
    except AssertionError as e:
        
        