import speedcom.tests.context as context

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
    xtrain_test = np.array([1, 2, 3])
    xtest_test  = np.array([4, 5, 6])
    ytrain_test = np.array([7, 8, 9])
    ytest_test  = np.array([10, 11, 12])



    assert os.path.isfile(filename), \
        "Function hasn't written feature results, or has not done so to the \
        correct directory."
    os.remove(filename)
