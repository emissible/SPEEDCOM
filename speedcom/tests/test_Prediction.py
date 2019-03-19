import numpy as np
from numpy.testing import assert_array_equal
#import speedcom.tests.context as context
import speedcom.tests.context as context
#import ontext.Prediction.Models as Models


test_SMILES = 'C1=CC=CC=C1'
test_model = context.Prediction.Models('..')

def test__init__():
    try:
        model = context.Prediction.Models('..')
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
        quantum yield
    """
    output = test_model.predict_qy(test_SMILES)
    assert isinstance(output, np.ndarray),\
        'output has to be array'
    assert_array_equal(test_shape, (1,1), 'output shape wrong')
    return
    
def test_predict_epsilon():
    """
    test the output type and shape of the predicted 
        epsilon
    """
    output = test_model.predict_epsion(test_SMILES)
    assert isinstance(output, np.ndarray),\
        'output has to be array'
    assert_array_equal(test_shape, (1,1), 'output shape wrong')
    return

def test_predict_all():
    """
    test the output shape and each output shape
    expected length of output is 2
    first output and second output length is 4 and 5
    """
    result = test_model.predict_all(test_SMILES)
    assert len(result) == 2
    assert len(result[0]) == 4
    assert len(result[1]) == 5
    return

def test_save_table_file():
    """
    test if table actually saves out
    """
    table, visual_data = test_models.predict_all(test_SMILES)
    target_fn = 'tempdata'
    try:
        test_model.save_table_file(target_fn, table)
    except Exception as e:
        assert isinstance(e, TypeError)
    assert os.path.isfile(target_fn), \
        "Function hasn't written feature results, or has not done so to the \
        correct directory."
    os.remove(target_fn)

    return
    
def test_save_visual():
    """
    test if visualization data successfully saves out    
    """
    table, visual_data = test_models.predict_all(test_SMILES)
    test_name = test_SMILES
    try:
        models.save_visual(best_name, visual_data)
    except Exception as e:
        assert isinstance(e, TypeError)
    assert os.path.isfile('%s_example_plot_data.txt' %s), \
        "Function hasn't written feature results, or has not done so to the \
        correct directory."
    os.remove('%s_example_plot_data.txt' %s)
    return