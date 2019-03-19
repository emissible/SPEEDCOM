import json
from keras.models import model_from_json
import numpy as np
from sklearn.linear_model import LinearRegression
import speedcom.tests.context as context
# import context.Models as Models

mod_util = context.ModUtil()

def test_coeff_determination():
    """
    Tests the function that calculates the R squared value
        between a true and predicted value.
    """
    # ytrue_test = np.array([5, 6, 7])
    # ypred_test = np.array([5, 6, 7])
    #
    # try:
    #     R_sq_test = coeff_determination(ytrue_test, ypred_test)
    # except (AssertionError):
    #     pass

    return

def test_plot_model_error():
    """
    Tests the function that plots the training and test data
        and the associated error.
    """
    # xtrain_test = np.array([1, 2, 3])
    # xtest_test  = np.array([4, 5, 6])
    # ytrain_test = np.array([7, 8, 9])
    # ytest_test  = np.array([10, 11, 12])
    # filename = 'model_error_test_plot.png'
    # # X = np.array([[1, 1], [1, 2], [2, 2], [2, 3]])
    # # y = np.dot(X, np.array([1, 2])) + 3
    # model = LinearRegression().fit(xtrain_test.reshape(-1,1), \
    #                                ytrain_test.reshape(-1,1))
    #
    #
    # try:
    #     plot_model_error(xtrain_test, x_test_test, ytrain_test, ytest_test)
    #
    #
    # assert os.path.isfile(filename), \
    #     "Function hasn't written feature results, or has not done so to the \
    #     correct directory."
    # os.remove(filename)

    return

def test_combine_columns():
    """
    Tests the function that combines an array of arrays of rows into
        columns.
    """
    nd_arrays = [ np.array([5,6,7]), np.array([1,2,3]) ]
    try:
       columns_combined = mod_util.combine_columns(nd_arrays)
    except Exception as e:
        assert isinstance(e, TypeError)
    assert isinstance(columns_combined, np.ndarray), \
        'Function not producing the correct output type.'
    col_dim = len(nd_arrays[0])
    assert len(columns_combined) == col_dim, \
        'Output array has a dimensional error.'

    return

def test_get_y_category():
    """
    Tests the function that discretizes the actual Y values into
        classes/bins.
    """
    y_actual = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

    try:
        disc = mod_util.get_y_category(y_actual, 1, 10)
    except Exception as e:
        assert isinstance(e, TypeError)
    assert len(disc) == len(y_actual), \
        'Lengths of input and ouput arrays must be equal.'

    return

def test_get_class_count():
    """
    Tests the function that counts the elements in each class of a
        discretized array and returns them as a dictionary.
    """
    categories = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4])
    try:
        counts = mod_util.get_class_count(categories)
    except Exception as e:
        assert isinstance(e, TypeError)
    assert isinstance(counts, dict), \
        'Function not properly outputting a dictionary'
    assert len(counts) <= len(categories), \
        'Length of class count cannot be greater than length of input array'

    return

def test_subsampling():
    """
    Tests the function that subsamples an array of data to ensure
        even distribution amongst classes.
    """
    dataset = np.array([[1,0],[2,0],[3,0],[3,0],[3,0],[3,0],[3,0],[4,1],\
                    [5,1],[6,1],[7,2],[8,2],[9,2],[10,3]])
    counts = {0: 7, 1: 3, 2: 3, 3: 1}
    try:
        ss_data = mod_util.subsampling(dataset, 1, 3, counts)
    except Exception as e:
        assert isinstance(e, TypeError)
    assert len(ss_data) <= len(dataset), \
        'Function should output a data set smaller than or equal to input data'
    assert len(ss_data[0]) == len(dataset[0]), \
        'Number of columns in input and output data should be equal'

    return

def test_onehot_encode():
    """
    Tests the function that encodes the y categories using the
        'one-hot' method.
    """
    y = np.array([0,1,3,2])
    num_class = 4

    try:
        one_hot = mod_util.onehot_encode_y(y, num_class)
    except Exception as e:
        assert isinstance(e, TypeError)
    assert isinstance(one_hot, np.ndarray), \
        'Function not outputting a numpy ndarray'

    return
