import numpy as np
from numpy.testing import assert_array_equal
#import speedcom.tests.context as context
import context.Prediction.Models as Models


test_SMILES = 'C1=CC=CC=C1'
test_model = Models()

def test__init__():
    try:
        model = context.Prediction.Models()
    except:
        raise RuntimeError('Error in constructing model')
    assert isinstance(model, context.Prediction.Models), \
            'molecule object constructed as wrong type.'

def test_predict_abs():
    """
    test the output type and shape of the predicted absorption 
        wavelength 
    """
    output = test_model.predict_abs(test_SMILES)
    assert isinstance(output, np.ndarray),\
        'output has to be array'
    assert_array_equal(test_shape, (1,1), 'output shape wrong')
    return
    
    
def test_predict_ems():
    """
    test the output type and shape of the predicted 
        emission wavelength 
    """
    output = test_model.predict_ems(test_SMILES)
    assert isinstance(output, np.ndarray),\
        'output has to be array'
    assert_array_equal(test_shape, (1,1), 'output shape wrong')
    return

def test_predict_quantum_yield():
    """
    test the output type and shape of the predicted 
        emission wavelength 
    """
    return
    