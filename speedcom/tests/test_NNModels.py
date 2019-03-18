import speedcom.tests.context as context

def test__init__():
    """
    Testing for correct construction of the Descriptor class.
    """
    try:
        mol_obj = context.NNModels.Descriptors()
    except:
        raise RuntimeError('Error in constructing molecule')
    assert isinstance(mol_obj, NNModel.Descriptors), \
        'molecule object constructed as wrong type.'

    return

def test_set_molecule():
    """
    Testing the function that reassigns the molecule
        object to a new SMILES string.
    """
    test_SMILES = 'C1=CC=CC=C1'
    test_SMILES2 = 'C[C@@H](CC1=CC=CC=C1)NC'
    molecule = context.NNModels.Descriptors(test_SMILES)
    try:
        mol_obj.set_molecule(test_SMILES2)
    except Exception as e:
        assert isinstance(e, TypeError)


    return

def test_get_properties():
    """
    Tests the function that returns a molecule's properties.
    """
    test_SMILES = 'C1=CC=CC=C1'
    test_mol = context.NNModels.Descriptors(test_SMILES)

    try:
        test_props = test_mol.get_properties()
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
    test_SMILES = 'C1=CC=CC=C1'
    test_mol = context.NNModels.Descriptors(test_SMILES)

    try:
        test_mol.get_features()
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
       test_fp = get_Morgan_fingerprint(test_usefeat)
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
        get_coulomb_matrix(smile)
    except Exception as e:
        pass
    return
