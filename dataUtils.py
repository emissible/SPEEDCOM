import numpy as np 
import pandas as pd
from sklearn.model_selection import train_test_split

class DataUtils:

  
  @staticmethod
  def readData(fname,sep='\t'):
    """
    Read data from a file, print first couple lines
    return nd array of the data    
    """
    dset = pd.read_csv(fname, sep=sep)
    assert len(dset) >0
    print(dset.head())
    return dset.values
  
  @staticmethod
  def get_xy(dataset, x_col_index, y_col_index):
    """
    select x, y from the nd array dataset, return two arrays
    x and y
    """
    return dataset[:,x_col_index], dataset[:,y_col_index]
  
  @staticmethod
  def splitData(x, y, ratio=0.1, random_state=42):
    """
    Split data for training and testing
    """
    return train_test_split(x, y, test_size=ratio, random_state=random_state)

  
  @staticmethod
  def get_wordmap(x_smiles):
    """
    get integer map for character and int
    return dict and max length of all str
    """
    charset = set(''.join(x_smiles.flatten()) + "!E")
    char_to_int = dict((c,i) for i,c in enumerate(charset))
    return char_to_int
  
  @staticmethod  
  def get_max_len(x_smiles):
    """
    find the max length of smiles of all cases
    """
    return max(map(lambda x: len(x), x_smiles.flatten()))
  
  @staticmethod
  def onehot_encoding(x_list, uniform_length, word_map):
    """
    Encode the input with one-hot encoding
    """
    one_hot =  np.zeros((x_list.shape[0], uniform_length , len(word_map)),dtype=np.int8)
    for i, item in enumerate(data_list):
        #encode the startchar
        one_hot[i,0,word_map["!"]] = 1
        #encode the rest of the chars
        for j,c in enumerate(item):
            one_hot[i,j+1,word_map[c]] = 1
           #Encode endchar
        one_hot[i,len(item)+1:,word_map["E"]] = 1
       #Return one-hot matrix
    return one_hot

  
  @staticmethod
  def numeric_encoding(x_list, uniform_length, word_map):
    """
    Encode the input data with numerical encoding
    """
    ret = np.ndarray((x_list.shape[0], uniform_length), dtype=np.int8)
    for i, item in enumerate(x_list):
        ret[i, 0] = word_map["!"]
        for j, c in enumerate(item):
            ret[i, j + 1] = word_map[c]
        ret[i, len(item)+1:] = word_map['E']
    return ret
