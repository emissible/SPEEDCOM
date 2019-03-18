import numpy as np 
import pandas as pd
import matplotlib as mplt
import matplotlib.pyplot as plt

class ModelUtils:
  
    @staticmethod
    def coeff_determination(y_true, y_pred):
        from keras import backend as K
        SS_res =  K.sum(K.square( y_true-y_pred ))
        SS_tot = K.sum(K.square( y_true - K.mean(y_true) ) )
        return ( 1 - SS_res/(SS_tot + K.epsilon()) )
    
    
  
    @staticmethod
    def plot_model_error(x_train, x_test, y_train, y_test, model, save_fig_fname, label):
        """
        input model and actual data, plot train, test predicted data and error
        """
        y_train = y_train.reshape(-1,1)
        y_test = y_test.reshape(-1,1)
    
        plt.figure(figsize=(8, 10), dpi=100)
        plt.subplot(211)
        plt.scatter(y_train, model.predict(X_train), color = 'r', label = 'train')
        plt.scatter(y_test, model.predict(X_test), color = 'blue', label = 'test')
        plt.xlabel('Actual %s' % label, fontsize=12)
        plt.ylabel('Predicted %s' %label, fontsize=12)
        plt.legend(loc='upper center')

        plt.subplot(212)
        plt.scatter(y_train, model.predict(X_train)-y_train, color = 'r', label = 'train', marker= 'x')
        plt.scatter(y_test, model.predict(X_test)-y_test, color = 'blue', label = 'test', marker = 'x')
        plt.axhline(0, ls='--')
        plt.xlabel('Actual %s' % label, fontsize=12)
        plt.ylabel('Prediction error %s' % label, fontsize=12)
        plt.legend(loc='upper center')
        
        if save_fig_fname is not None:
            plt.savefig(save_fig_fname)
        else:
            plt.show()
    
    @staticmethod
    def combine_columns(nd_arrays):
        """
        combine nd array as columns
        nd_arrays: tuple of nd arrays of same dimension
        return: np.array
        """
        return np.column_stack((nd_arrays))
    
    @staticmethod
    def get_y_category(y_actual,min_edge, bin_width):
        """
        get assing Y into classes based on y actual value

        y_actual -- np.array
        min_edge -- int/float
        bin_width --int/foat

        return Y_category -- np.array:  class variable of Y
        """
        Y_category = (y_actual - min_edge) // bin_width
        return Y_category
    
  
    @staticmethod
    def get_class_count(Y_category):
        """
        get counts of each class and return as dictionary

        Y_category -- np.array: categorical y
        return cls_counts -- dict: {class(key): counts(value)}

        """
        cls, counts = np.unique(Y_category, return_counts = True)
        cls_counts = dict(zip(cls, counts))
        return cls_counts
    
    
  
    @staticmethod
    def subsampling(dataset, class_column_index, class_max_count, class_dict):
        """
        subsample each class to make sure every class has closer counts

        dataset -- nd.array: includes x column, y(category) etc
        class_column_index -- int: index of the category data in dataset
        class_max_count -- int: maximum number of counts allowed in each class 
        class_dict -- dict: {class: count}

        return subsampled data -- nd array
        """
        out = []
        for row in dataset:
            cls = row[class_column_index]
            rInt = np.random.randint(0, class_dict[cls])
            if rInt <= class_max_count:
                out.append(row)
        return np.array(out)
    
    @staticmethod
    def onehot_encode_y(y, num_class):
        """
        one-hot encoding of categorical y
        """
        one_hot = np.zeros((y.shape[0],num_class),dtype=np.int8)
        for index, cls in enumerate(y):
          one_hot[index, int(cls)] = 1
        return one_hot

    @staticmethod  
    def coeff_determination(y_true, y_pred):
        """
        calculates R square from y actual and y predict
        """
        from keras import backend as K
        SS_res =  K.sum(K.square( y_true-y_pred ))
        SS_tot = K.sum(K.square( y_true - K.mean(y_true) ) )
        return ( 1 - SS_res/(SS_tot + K.epsilon()) )
    
    @staticmethod
    def get_lr_metric(optimizer):
        def lr(y_true, y_pred):
            return optimizer.lr
        return lr
