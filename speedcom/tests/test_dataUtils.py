"""
Unit tests for the dataUtils.py file
"""
import json
import pandas as pd
import numpy as np
import os
import speedcom.tests.context as context
#from context import dataUtils.DataUtils as DataUtils

DataUtils =  context.dataUtils.DataUtils()

def test_readData():
    """
    Test function for readData
    """
    test_filename = 'test_readData_temp.tsv'
    if os.path.exists(test_filename):
        os.remove(test_filename)

    # test raise exception for random names
    try:
        DataUtils.readData(test_filename)
    except Exception as e:
        assert isinstance(e, AssertionError)

    test_df = pd.DataFrame([1,2])
    test_df.to_csv(test_filename)

    # test outout type
    df_read = DataUtils.readData(test_filename)
    assert isinstance(df_read, np.ndarray)
    return

def test_get_xy():
    """
    Test function for get_xy
    """

    data
    # test assertion
    try:
        DataUtils.get_xy(test_filename)
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
    return
