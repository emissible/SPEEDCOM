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
    test function for load word map from an exsting json
    """
    with open('temp.json', 'w') as fp:
        json.dump({1:'#',2:'!'})
    try:
        word_map = DataUtils.load_wordmap_from_json('temp.json')
    except Exception as e:
        assert isinstance(e, AssertionError)
    assert isinstance(word_map, dict), \
        'loaded json does not create a dict'
    assert len(word_map) == 2
    os.remove('temp.json')
    
    return
    

def test_save_wordmap_json():
    """
    test funtion for saving wordmap into json file 
    """
    wordmap = {1:'#',2:'!'}
    try:
        DataUtils.save_wordmap_json(wordmap, 'temp_saved.json')
    except Exception as e:
        assert isinstance(e, TypeError)
    assert os.path.isfile('temp_saved.json'), \
        'did not successfuly save to file'
    os.remove('temp.json')
    
    return
    
def test_decode_num_smiles():
    """
    test if convert a numeric input back into strings
    """
    test_numeric_smiles = np.array([0,1,1,2,3,4])
    test_rev_wordmap = {0:'!', 1: 'C', 2: '#', 3: '=', 4:'E'}
    output = DataUtils.decode_num_smiles(test_numeric_smiles, test_rev_wordmap)
    assert len(output[0]) == 4, 'not getting correct smiles decoding'
    return
    
def test_reverse_wordmap():
    """
    test if convert a word map to reverse word map
    """
    wordmap = {'C':0,'#':1,'!':2}
    rev_wordmap = DataUtils.reverse_wordmap(reverse_wordmap(word_map))
    assert len(rev_wordmap) == 3
    return 
    
    
def test_numeric_encoding():
    """
    test if getting correct numeric encoding
    """
    smiles_list = np.array(['CC=C#C'])
    test_rev_wordmap = {0:'!', 1: 'C', 2: '#', 3: '=', 4:'E'}
    uniform_length = 8
    numeric = DataUtils.numeric_encoding(x_list, uniform_length, word_map)
    assert len(numeric[0]) == 8
    return
    
    
    