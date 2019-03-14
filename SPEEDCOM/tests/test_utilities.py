import numpy as np



def test_remove_deliminators():
    """
    Tests the function that removes delimiters from
        a list of numbers.
    """

    tester = ['6,230', '4,890']

    try:
        my_array = remove_deliminators(tester)
    except Exception as e:
        assert isinstance(e, TypeError)

    assert isinstance(my_array, np.array), \
        'Wrong Type: function has not returned a numpy array.'
    for i in my_array:
        assert isinstance(float, float64), \
            'value at index' + str()



    return
