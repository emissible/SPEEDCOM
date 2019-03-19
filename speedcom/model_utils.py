import numpy as np
import pandas as pd
import matplotlib as mplt
import matplotlib.pyplot as plt

class ModelUtils:

    @staticmethod
    def coeff_determination(y_true, y_pred):
        """
        Calculates the R squared of true value and predicted value
            with keras backend.

        Args:
        -----
            y_true (np.ndarray) -- an array of actual y-values
            y_pred (np.ndarray) -- an array of predicted y-values

        Returns:
        -----
            R_squared (float) -- the R squared value between the
                input arrays.
        """
        # Assertions
        # assert y_true.shape[0] = y_pred.shape[0], \
        #     'Lengths of input arrays must be the same.'
        # assert isinstance(y_true, np.ndarray), \
        #     'Input y value not a numpy ndarray.'
        # assert isinstance(y_pred, np.ndarray), \
        #     'Input y pred not a numpy ndarray.'
        # Functionality
        from keras import backend as K
        SS_res =  K.sum(K.square( y_true-y_pred ))
        SS_tot = K.sum(K.square( y_true - K.mean(y_true) ) )
        R_squared = (1 - SS_res/(SS_tot + K.epsilon()))

        return R_squared

    @staticmethod
    def plot_model_error(x_train, x_test, y_train, y_test, model, \
                         save_fig_fname, label):
        """
        Plots the training and test data, and the associated error
            produced by a certain model.

        Args:
        -----
            x_train (np.ndarray) -- the training x values
            x_test (np.ndarray) -- the testing x values
            y_train (np.ndarray) -- the training y values
            y_test (np.ndarray) -- the testing y values
            model (object i.e. Keras.model) -- the model that will
                be used to generate the predicted values.
            save_fig_fname (str) -- the desired output filename. Must
                include the .png extension. If None, will just show
                the plot rather than saving it to file.
            label (str) -- the context of the data being plotted
                i.e. 'Absorption' or 'Emission'

        Returns:
        --------
            Either an output .png file containing the plot if
                save_fig_fname is specified, else no return (just
                shows the plot).
        """
        # Assertions
        assert isinstance(x_train, np.ndarray)
        assert isinstance(x_test, np.ndarray)
        assert isinstance(y_train, np.ndarray)
        assert isinstance(y_test, np.ndarray)
        assert isinstance(save_fig_fname, (str, type(None))), \
            'Wrong Type: desired output file name must be a string'
        if type(save_fig_fname) is str:
            assert save_fig_fname.endswith('.png'),\
                'output_file string must include the .png extension'

        # Functionality
        y_train = y_train.reshape(-1,1)
        y_test = y_test.reshape(-1,1)
        plt.figure(figsize=(8, 10), dpi=100)
        plt.subplot(211)
        plt.scatter(y_train, model.predict(x_train), color='r', label='train')
        plt.scatter(y_test, model.predict(x_test), color='blue', label='test')
        plt.xlabel('Actual %s' % label, fontsize=12)
        plt.ylabel('Predicted %s' %label, fontsize=12)
        plt.legend(loc='upper center')

        plt.subplot(212)
        plt.scatter(y_train, model.predict(x_train)-y_train, color = 'r', \
                label = 'train', marker= 'x')
        plt.scatter(y_test, model.predict(x_test)-y_test, color = 'blue', \
                label = 'test', marker = 'x')
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
        Combine the columns from an array of arrays into arrays
            themselves.

        Args:
        -----
            nd_arrays (np.ndarray) -- an array of arrays of rows
                to be combined into columns.

        Returns:
        -----
            columns_combined (np.array) -- an array of arrays of
                the combined columns.
        """
        # Assertions
        assert isinstance(nd_arrays, list), \
           'Input must be a list'
        for i in nd_arrays:
            assert isinstance(i, np.ndarray), \
                'Each element in list must be a numpy ndarray'
        l_prev = len(nd_arrays[0])
        for i in range(1, len(nd_arrays)):
            l = len(nd_arrays[i])
            assert l == l_prev, \
                'Array at index ' + str(i) + ' is the wrong size. ' \
                'All arrays in input array must be the same length'
            l_prev = l
        # Functionality
        columns_combined = np.column_stack((nd_arrays))

        return columns_combined

    @staticmethod
    def get_y_category(y_actual, min_edge, bin_width):
        """
        Divide continuous y_actual into discrete classes/bins.

        Args:
        -----
            y_actual (np.ndarray) -- the y values to be discretized
            min_edge (int/float)  -- the value at which to start the
                discretization.
            bin_width (int/foat)  -- the width of the bins.

        Returns:
        --------
            Y_category (np.ndarray) -- the actual y values divided
                discrete classes.

        """
        # Assertions
        assert max(y_actual) > min_edge, \
            'min_edge has to be smaller than maximum y_actual'
        assert isinstance(y_actual, np.ndarray), \
            'Input array must be a numpy ndarray.'
        assert isinstance(min_edge, (int, float)), \
            'min_edge should be a float or an int.'
        assert isinstance(bin_width, (int, float)), \
            'bin_width should be a float or an int.'
        # Functionality
        Y_category = (y_actual - min_edge) // bin_width

        return Y_category

    @staticmethod
    def get_class_count(Y_category):
        """
        Get counts of each class and returns it as dictionary.

        Args:
        -----
            Y_category (np.array) -- the actual y values divided
                into categrories.

        Returns:
        --------
            cls_counts (dict) -- {class(key): counts(value)}
        """
        # Assertions
        assert isinstance(Y_category, np.ndarray), \
            'Input must be a numpy ndarray.'
        cls, counts = np.unique(Y_category, return_counts = True)
        cls_counts = dict(zip(cls, counts))

        return cls_counts

    @staticmethod
    def subsampling(dataset, class_column_index, class_max_count, class_dict):
        """
        Subsamples each class to ensure even distribution amongst
            classes.

        Args:
        -----
            dataset (nd.array) -- includes x column, y(category) etc
            class_column_index (int) -- index of the category data
                in dataset.
            class_max_count (int) -- maximum number of counts
                allowed in each class.
            class_dict (dict) -- a dictionary of the form
                {class: count}.

        Returns:
        --------
            ss_data (np.ndarray) -- an array containing the subsampled
                data.
        """
        out = []
        for row in dataset:
            cls = row[class_column_index]
            rInt = np.random.randint(0, class_dict[cls])
            if rInt <= class_max_count:
                out.append(row)
        ss_data = np.array(out)

        return ss_data

    @staticmethod
    def onehot_encode_y(y, num_class):
        """
        One-hot encoding of the categories y

        Args:
        -----
            y (np.ndarray)  --  the y values
            num_class (int) --  the number of classes y is divided into

        return:
            one_hot (np.array) -- an array that contains the encoded y
                data.
        -----
        """
        # Assertions
        assert isinstance(y, np.ndarray), \
            'y must be a numpy ndarray'
        assert isinstance(num_class, int), \
            'num_class must be an int'
        # Functionality
        one_hot = np.zeros((y.shape[0],num_class),dtype=np.int8)
        for index, cls in enumerate(y):
          one_hot[index, int(cls)] = 1

        return one_hot

    @staticmethod
    def get_lr_metric(optimizer):
        """
        Returns the learning rate from optimizer.
        Args:
        -----
            optimizer (keras.optimizer) -- helps to optimize the
                compiling of the keras model.
        return:
        -----
            lr (float) -- the learning rate of the optimizer
        """
        def lr(y_true, y_pred):
            return optimizer.lr

        return lr
