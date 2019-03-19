import numpy as np
import speedcom.tests.context as context

# Initialize the class with a test moleucle.
test_SMILES = 'C1=CC=CC=C1'
molecule = context.NNModels.Descriptors(test_SMILES)

def test__init__():
    """
    Testing for correct construction of the Descriptor class.
    """
    try:
        mol = context.NNModels.Descriptors()
    except:
        raise RuntimeError('Error in constructing molecule')
    # Checking that the model object is the correct type
    assert isinstance(mol, context.NNModels.Descriptors), \
        'molecule object constructed as wrong type.'

    return

def test_set_molecule():
    """
    Testing the function that reassigns the molecule
        object to a new SMILES string.
    """
    test_SMILES2 = 'C[C@@H](CC1=CC=CC=C1)NC'
    try:
        molecule.set_molecule(test_SMILES2)
    except Exception as e:
        assert isinstance(e, TypeError)

    return

def test_get_properties():
    """
    Tests the function that returns a molecule's properties.
    """
    molecule.set_molecule(test_SMILES)
    try:
        test_props = molecule.get_properties()
        test_indiv = molecule.get_properties(property_name='exactmw')
        test_list = molecule.get_properties(property_name=['exactmw','NumRings'])
    except Exception as e:
        assert isinstance(e, TypeError)
    # Checking that each element in the list is a float
    index = 0
    for prop_val in test_props.values():
        assert isinstance(prop_val, float), \
            "Error in output dictionary; value at index " \
            + str(index) + " isn't a float."
        index += 1
    # Checking that the individual property returned is correct
    assert isinstance(test_indiv, float), \
        'exactmw must be a float'
    assert test_indiv == 78.046950192 , \
        "exactmw doesn't match that of benzene"
    assert isinstance(test_list, dict), \
        'Ouput of function upon imputting a list must be a dict'
    assert test_list['exactmw'] == 78.046950192, \
        'Output of function for these exactmw not correct'
    assert test_list['NumRings'] == 1.0, \
        'Output of function for these exactmw not correct'

    return

def test_get_features():
    """
    Tests the function that returns a dictionary of features for a
        given molecule.
    """
    try:
        feats = molecule.get_features()
    except Exception as e:
        assert isinstance(e, TypeError)
    assert isinstance(feats, dict), \
        'Function should be outputting a dict'

    return

def test_get_Morgan_fingerprint():
    """
    Tests the function that generates a unique Morgan
        fingerprint for a molecule.
    """
    test_usefeat = True
    try:
        test_fp = molecule.get_Morgan_fingerprint(use_feat=test_usefeat)
    except Exception as e:
        assert isinstance(e, TypeError)
    # Checking that each element in the list is a float
    index = 0
    for i in test_fp:
        assert isinstance(i, int), \
            'Error in output: value at index ' + str(index) + \
            ' isnt an integer'
        index += 1

    return

def test_get_coulomb_matrix():
    """
    Tests the function that returns the Coulomb matrix of a molecule.
    """
    try:
        out_cm = molecule.get_coulomb_matrix()
        out_eig = molecule.get_coulomb_matrix(output_eigval=True)
    except Exception as e:
        assert isinstance(e, TypeError)
    # Checking outputs
    assert isinstance(out_cm, np.ndarray), \
        'Function not outputting a numpy array'
    assert isinstance(out_eig, np.ndarray), \
        'Function not outputting a numpy array'

    return
