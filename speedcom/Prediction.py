import numpy as np
import pandas as pd

from scipy.stats import boxcox
from dataUtils import DataUtils
from keras.models import model_from_json
from scipy.special import inv_boxcox

class Models:
    """
    A wrapper class using keras pre-trained models for 
    specrtra prediciton.
    Initilized with loaded model for absorption, emission,
        epsilon and quantumn yield; loaded preset SMILES word map
    
    Attributes:
    -----
        Models
        _load_model_weight
        UNIFORM_LENGTH

    Methods:
    -----
        predict_abs
        predict_ems
        predict_quantum_yield
        predict_epsilon
        predict_all
        save_table_file
        save_visual
        _predict
        _reshapeX
        _load_model_weight

    """
    
    # set uniform length for all SMILES string
    UNIFORM_LENGTH = 279
    
    def __init__(self):
        """
        model constructor with preloaded features
        """
        self.model_abs = self._load_model_weight('model_smiles_cnn.json',\
            'weights.best.hdf5')
        self.model_ems = self._load_model_weight('model_smiles_ems_dropna.json',\
            'ems_dropna.best1.hdf5')
        self.model_epsilon = self._load_model_weight('model_smiles_epsilon_lstm.json',\
            'epsilon.best_-1.hdf5')
        self.model_qy = self._load_model_weight('model_lstm_qy.json', \
            'weights.qy_best.hdf5')
        self.word_map = DataUtils.load_wordmap_from_json('smiles_wordmap.json')

    def predict_abs(self, x):
        """
        predict absorption wavelength
        Args:
        x -- SMILES(str)
        """
        return self._predict(self.model_abs, x)

    def predict_ems(self, x):
        """
        predict emission wavelength
        Args:
        x -- SMILES(str)
        """
        return self._predict(self.model_ems, x)
        
    def predict_quantum_yield(self, x):
        """
        predict quantum yield
        Args:
        x -- SMILES(str)
        """
        return self._predict(self.model_qy, x)

    def predict_epsilon(self, x, boxcox_lambda=0.2):
        """
        predict epsilon
        Args:
        x -- SMILES(str)
        boxcox_lambda -- int/float(preset to 0.2 for preloaded model)
        """
        y = self._predict(self.model_epsilon, x)
        return inv_boxcox(y, boxcox_lambda)

        
    def predict_all(self, x, boxcox_lambda=0.2):

        """
        predict abs, ems, epsilon, quantum yield
        Args:
        -----
        x -- SMILES(str)
        boxcox_lambda --int/float
        
        return:
        -----
        table -- np.array: data to fill the frontend table
        data_plot -- np.array: data for ploting
        """
        abs_l = self.predict_abs(x)
        ems_l = self.predict_ems(x)
        epsilon = self.predict_epsilon(x, boxcox_lambda)
        qy = self.predict_quantum_yield(x)
        table = np.column_stack([abs_l, ems_l, epsilon, qy])
        data_plot = np.column_stack([np.array(x), abs_l[0], [1.0], ems_l[0],[1.0]])
        return table, data_plot
    
    def save_table_file(self, x, table):
        """
        save table to file as txt with smiles name
        Args:
        -----
        table: np.array

        """
        table_df = pd.DataFrame(table, columns=['lambda_abs', 'lambda_ems','epsilon', 'quantum_yield'])
        table_df.to_csv('%s_example_table.txt' %x, sep = '\t', index=False)
        return
    
    def save_visual(self, x, data_for_plot):
        """
        save data_for_plot to smiles name
        Args:
        -----
        data_for_plot: np.array

        """
        plot_df = pd.DataFrame(data_for_plot, columns=['SMILES', 'abs_wl_max',\
           'abs_intensity','ems_wl_max','ems_intensity'])
        plot_df.to_csv('%s_example_plot_data.txt' % x, sep = '\t', index=False)
        return
    
    def _predict(self, model, x):
        """
        prediction
        Args:
        -----
        model: keras.model
        x: str 
        return:
        np.array(float/int)
        """
        x_array = self._reshapeX(x)
        x_numeric = DataUtils.numeric_encoding(x_array, self.UNIFORM_LENGTH, self.word_map) 
        return model.predict(x_numeric)

    def _reshapeX(self, x):
        """
        expand dimension
        """
        return np.array([x])

    def _load_model_weight(self, model_file, weight_file):
        """
        load model and its weight
        Args:
        model_file: json filename
        weight_file: hdf5 filename
        return:
        loaded_model: keras.model

        """
        json_file = open(model_file, 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        loaded_model = model_from_json(loaded_model_json)
        loaded_model.load_weights(weight_file)
        return loaded_model

if __name__ == '__main__':
    models = Models()
    input_smiles = 'C1=CC=CC=C1'
    table, visual_data = models.predict_all(input_smiles)
    models.save_table_file(input_smiles, table)
    models.save_visual(input_smiles, visual_data)
    print(table)
    print(visual_data)
