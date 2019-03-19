import numpy as np
import pandas as pd
import json
import os
from sklearn.model_selection import train_test_split

class DataUtils:


    @staticmethod
    def readData(fname,sep='\t'):
        """
        Read data from a file, print first couple lines
        and return the nd array of the data
        Args:
        -----
        fname: str

        return:
        -----
        np.array
        """
        assert os.path.exists(fname)
        dset = pd.read_csv(fname, sep=sep)
        assert len(dset) >0
        print(dset.head())
        return dset.values

    @staticmethod
    def get_xy(dataset, x_col_index, y_col_index):
        """
        select x, y from the nd array dataset, return two arrays
        x and y
        Args:
        -----
        dataset: np.array
        x_col_index: int
        y_col_index: int

        return:
        -----
        x: np.array
        y: np.array
        """
        assert isinstance(dataset, np.ndarray)
        return dataset[:,x_col_index], dataset[:,y_col_index]

    @staticmethod
    def splitData(x, y, ratio=0.1, random_state=42):
        """
        Split data for training and testing
        Args:
        -----
        x: np.array
        y: np.array
        ratio: float (0 < ratio < 1)
        random_state: int

        return:
        -----
        X_train: np.array
        X_test: np.array
        y_train: np.array
        y_test: np.array

        """
        assert x.shape[0] == y.shape[0]
        return train_test_split(x, y,test_size=ratio,random_state=random_state)


    @staticmethod
    def get_wordmap(x_smiles):
        """
        get numeric map for characters
        return dict of this map
        start and end as "!" and "E"
        need to sort it to keep consistency of wordmap
        as trained models are sensitive to numeric encoding
        Args:
        -----
        x_smiles: np.array of str
        return:
        -----
        char_to_int: dict {char : int}
        """
        charset = sorted(list(set(''.join(x_smiles.flatten()) + "!E")))
        char_to_int = dict((c,i) for i,c in enumerate(charset))
        return char_to_int

    @staticmethod
    def get_max_len(x_smiles):
        """
        find the max length of smiles of all cases
        Args:
        -----
        x_smiles:

        return:
        -----
        max_len: int
        """
        return max(map(lambda x: len(x), x_smiles.flatten()))

    @staticmethod
    def onehot_encoding(x_list, uniform_length, word_map):
        """
        Encode the input with one-hot encoding and pad it to
        uniform length each row
        Args:
        -----
        x_list: np.array of str
        uniform_length: int
        word_map: dict

        return:
        -----
        np.array
        """
        one_hot = np.zeros((x_list.shape[0],uniform_length,len(word_map)),\
                dtype=np.int8)
        for i, item in enumerate(x_list):
            #encode the startchar
            one_hot[i,0,word_map["!"]] = 1
            #encode the rest of the chars
            for j,c in enumerate(item):
                one_hot[i,j+1,word_map[c]] = 1
               #Encode endchar
            one_hot[i,len(item)+1:,word_map["E"]] = 1

        return one_hot


    @staticmethod
    def numeric_encoding(x_list, uniform_length, word_map):
        """
        Encode the input smiles with numerical encoding
        according to the word map and pad them to
        uniform length
        Args:
        -----
        x_list: np.array
        uniform_length: int
        word_map: dict

        return:
        -----
        ret: np.array
        """
        ret = np.ndarray((x_list.shape[0], uniform_length), dtype=np.int8)
        for i, item in enumerate(x_list):
            ret[i, 0] = word_map["!"]
            for j, c in enumerate(item):
                ret[i, j + 1] = word_map[c]
            ret[i, len(item)+1:] = word_map['E']
        return ret

    @staticmethod
    def reverse_wordmap(word_map):
        """
        based on char-to-int to get the map of
        int-to-char
        Args:
        -----
        word_map: {char: int}

        return:
        -----
        rev_wordmap: {int: char}
        """
        rev_wordmap={}
        for key, val in word_map.items():
            rev_wordmap[val] = key
        return rev_wordmap


    @staticmethod
    def decode_num_smiles(numeric_smiles_list, rev_wordmap):
        """
        given numeric encoded smiles arrays, decode the numeric
        back into its corresponding char and reconstruct the smiles
        Args:
        -----
        numeric_smiles_list: np.array (dtypes: int)
        rev_wordmap: dict (dtypes {int: char})

        return:
        -----
        np.array (dtype: str)
        """
        smiles_list = []
        for numeric_smiles in numeric_smiles_list:
            smiles = ''
            for num in numeric_smiles:
                char = rev_wordmap[num]
                if (char !='!' and char != 'E'):
                    smiles += char
            smiles_list.append(smiles)
        return np.array(smiles_list)

    @staticmethod
    def save_wordmap_json(wordmap, json_fname):
        """
        save wordmap dictionary to json file
        Args:
        -----
        word_map: dict

        return:
        -----
        json_fname: str
        """
        with open(json_fname, 'w') as fp:
            json.dump(wordmap, fp)

    @staticmethod
    def load_wordmap_from_json(json_fname):
        """
        load saved word map from json
        Args:
        -----
        json_fname: str

        return:
        -----
        dict
        """
        assert os.path.exists(json_fname), \
            'json_fname does not exist'
        with open(json_fname, 'r') as f:
            word_map_loaded = json.load(f)
        return word_map_loaded
