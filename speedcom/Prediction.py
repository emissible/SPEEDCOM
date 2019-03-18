import numpy as np
import pandas as pd

from scipy.stats import boxcox
from dataUtils import DataUtils
from keras.models import model_from_json
from scipy.special import inv_boxcox

class Models:
    UNIFORM_LENGTH = 279
    
    def __init__(self):
        """
        
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
        return self._predict(self.model_abs, x)

    def predict_ems(self, x):
        return self._predict(self.model_ems, x)
        
    def predict_quantum_yield(self, x):
        return self._predict(self.model_qy, x)

    def predict_epsilon(self, x, boxcox_lambda=0.2):
        y = self._predict(self.model_epsilon, x)
        return inv_boxcox(y, boxcox_lambda)

        
    def predict_all(self, x, boxcox_lambda=0.2):
        abs_l = self.predict_abs(x)
        ems_l = self.predict_ems(x)
        epsilon = self.predict_epsilon(x, boxcox_lambda)
        qy = self.predict_quantum_yield(x)
        table = np.column_stack([abs_l, ems_l, epsilon, qy])
        data_plot = np.column_stack([np.array(x), abs_l[0], [1.0], ems_l[0],[1.0]])
        return table, data_plot
    
    def save_table_file(self, x, table):
        table_df = pd.DataFrame(table, columns=['lambda_abs', 'lambda_ems','epsilon', 'quantum_yield'])
        table_df.to_csv('%s example_table.txt' %x, sep = '\t', index=False)    
    
    def save_visual(self, x, data_for_plot):
        plot_df = pd.DataFrame(data_for_plot, columns=['SMILES', 'abs_wl_max',\
           'abs_intensity','ems_wl_max','ems_intensity'])
        plot_df.to_csv('%s example_plot_data.txt' % x, sep = '\t', index=False)

    
    def _predict(self, model, x):
        x_array = self._reshapeX(x)
        x_numeric = DataUtils.numeric_encoding(x_array, self.UNIFORM_LENGTH, self.word_map) 
        return model.predict(x_numeric)

    def _reshapeX(self, x):
        return np.array([x])

    def _load_model_weight(self, model_file, weight_file):
        json_file = open(model_file, 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        loaded_model = model_from_json(loaded_model_json)
        loaded_model.load_weights(weight_file)
        return loaded_model

if __name__ == '__main__':
    models = Models()
    input = 'C1=CC=CC=C1'
    table, visual_data = models.predict_all(input)
    models.save_table_file(input, table)
    models.save_visual(input, visual_data)
    print(table)
    print(visual_data)
