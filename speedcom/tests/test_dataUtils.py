"""
Unit tests for the dataUtils.py file
"""
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

    data = 1
    # test assertion
    try:
        DataUtils.get_xy(data, 1, 1)
    except Exception as e:
        assert isinstance(e, AssertionError)

    return

def test_splitData():
    """
    Test function for splitData
    """
    x = pd.DataFrame([1])
    y = pd.DataFrame([1, 2])
    # test assertion
    try:
        DataUtils.splitData(x,y)
    except Exception as e:
        assert isinstance(e, AssertionError)

    return

def test_get_wordmap():
    """
    Test function for splitData
    """
    SMILES = 'CC'
    # test output file type
    result = DataUtils.get_wordmap(SMILES)
    assert isinstance(result, dict)

    return
