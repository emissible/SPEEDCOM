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
    
    
def test_load_wordmap_from_json():
    """
    test function for loading dictionary from json file
    """
    content = {'C':0, '#':1, '!':2, 'E':3}
    with open('temp.json', 'w') as fp:
        json.dump(content, fp)
    try:
        word_map = DataUtils.load_wordmap_from_json('temp.json')
    except Exception as e:
        assert isinstance(e, AssertionError)
    assert isinstance(word_map, dict), \
        'loaded json does not create a dict'
    assert len(word_map) == 4
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
        assert isinstance(e,FileNotFoundError)
    assert os.path.isfile('temp_saved.json'), \
        'did not successfuly save to file'
    os.remove('temp_saved.json')
    
    return
    
def test_decode_num_smiles():
    """
    test if convert a numeric input back into strings
    """
    test_numeric_smiles = np.array([[0,1,1,2,3,4]])
    test_rev_wordmap = {0:'!', 1: 'C', 2: '#', 3: '=', 4:'E'}
    output = DataUtils.decode_num_smiles(test_numeric_smiles, test_rev_wordmap)
    assert len(output[0]) == 4, 'not getting correct smiles decoding'
    return
    
def test_reverse_wordmap():
    """
    test if convert a word map to reverse word map
    """
    wordmap = {'C':0,'#':1,'!':2}
    rev_wordmap = DataUtils.reverse_wordmap(wordmap)
    assert len(rev_wordmap) == 3
    return 
    
    
def test_numeric_encoding():
    """
    test if getting correct numeric encoding
    """
    smiles_list = np.array(['CC=C#C'])
    test_wordmap = {'C':0,'#':1,'!':2, '=':3, 'E':4}
    uniform_length = 10
    numeric = DataUtils.numeric_encoding(smiles_list, uniform_length, test_wordmap)
    assert len(numeric[0]) == 10

    return

def test_get_wordmap():
    """
    Test function for get_wordmap
    """
    SMILES = np.array(['CC'])
    # test output file type
    result = DataUtils.get_wordmap(SMILES)
    assert isinstance(result, dict)

    return

def test_get_max_len():
    """
    Test function for get_max_len
    """
    SMILES = np.array(['CCC', 'CC'])
    max_len = DataUtils.get_max_len(SMILES)
    assert max_len == 3
    return

def test_onehot_encoding():
    """
    Test function for onehot_encoding
    """
    wordmap = {'C':0, '#':1, '!':2, 'E':3}
    x = np.array(['#C!'])
    result = DataUtils.onehot_encoding(x,5,wordmap)
    assert result.shape == (1,5,4)

    return
    
    
    