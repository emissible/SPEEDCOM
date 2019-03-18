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
    assert isinstance(mol, NNModel.Descriptors), \
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
    try:
        test_props = molecule.get_properties()
    except Exception as e:
        assert isinstance(e, TypeError)
    index = 0
    for prop_val in test_props.values():
        assert isinstance(prop_val, float), \
            "Error in output dictionary; value at index " \
            + str(index) + " isn't a float."
        index += 1

    return

def test_get_features():
    """
    Tests the function that returns a dictionary of features for a
        given molecule.
    """
    try:
        molecule.get_features()
    except Exception as e:
        assert isinstance(e, TypeError)

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
    assert isinstance(out_cm, np.ndarray), \
        'Function not outputting a numpy array'
    assert isinstance(out_eig, np.ndarray), \
        'Function not outputting a numpy array'

    return

# def test___get_charges_coords():
#     """
#     Tests the function that generates a pandas data frame of the
#         charges and coordinates of each atom in a molecule.
#     """
#     try:
#         out_charge_coords = molecule.__get_charges_coords()
#     except Exception as e:
#         assert isinstance(e, TypeError)
#     assert isinstance(out_charge_coords, pd.DataFrame), \
#         'Function not outputting a pandas DataFrame'
#
#     return
