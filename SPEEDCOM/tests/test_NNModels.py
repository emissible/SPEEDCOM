from SPEEDCOM.NNModels import Descriptors

# Defining some test SMILES strings
SMILES = 'C1=CC=CC=C1'       # benzene
SMILES2 = 'CC1=CC=CC=C1	'    # toluene

def test__init__():
    """
    Testing for proper construction of the class
    """


    return

def test_set_molecule():
    """
    Testing the function that reassigns the molecule
        object to a new SMILES string.
    """

    molecule = Descriptors(SMILES)
    try:
        molecule.set_molecule(SMILES2)
    except Exception as e:
        assert isinstance(e, TypeError)

    molecule = Descriptors(SMILES2)
    assert isinstance(molecule, NNModels.Descriptors), \
        'Wrong Type: should be of type NNModels.Descriptors'

    return

def test_get_properties():

    try:


    return
