import numpy as np
import os
import pandas as pd
import speedcom.tests.context as context
import rdkit.Chem.Draw as draw


def test_remove_deliminators():
    """
    Tests the function that removes delimiters from
        a list of numbers.
    """
    test_list = ['6,230', '4,890']
    try:
        my_array = context.utilities.remove_deliminators(test_list)
    except Exception as e:
        assert isinstance(e, TypeError)

    assert isinstance(my_array, np.ndarray), \
        'Wrong Type: function has not returned a numpy array.'
    index = 0
    for i in my_array:
        assert isinstance(i, float), \
            'value at index ' + str(index) + ' is not a float'
        index += 1

    assert context.utilities.remove_deliminators(['1,000,000']) == [1000000], \
        'Not correctly removing 2 delimiters from a number.'
    assert context.utilities.remove_deliminators(['0.1']) == [0.1], \
        'Incorrectly removing the decimal point from a float.'
    assert context.utilities.remove_deliminators(['1234']) == [1234], \
        'Incorrecly changing value of number passed'
    assert context.utilities.remove_deliminators(['1,234.5']) == [1234.5], \
        'Incorrecly changing value of number passed'

    return

def test_remove_cations():
    """
    Tests the function that removes counterions from a molecule's
        SMILES string.
    """
    # Shouldn't modify this SMILES...
    test_SMILES = 'C[C@@H](CC1=CC=CC=C1)NC'
    # ... or this one ...
    test_SMILES2 = 'CCN(CC)C1=CC2=C(C=C1)N=C3C=CC(=[N+](CC)CC)C=C3O2'
    # ... but should modify this one...
    test_SMILES3 = 'C1=CC2=C(C=C1N)SC3=CC(=[NH2+])C=CC3=N2.[Cl-]'
    try:
        removed1 = context.utilities.remove_cations(test_SMILES)
    except Exception as e:
        assert isinstance(e, TypeError)

    removed2 = context.utilities.remove_cations(test_SMILES2)
    removed3 = context.utilities.remove_cations(test_SMILES3)

    if not removed1 == test_SMILES:
        raise RuntimeError("Error: function removing \
            counterions that don't exist!")
    elif not removed2 == test_SMILES2:
        raise RuntimeError("Error: function removing \
            counterions that it shouldn't!")
    elif removed3 == test_SMILES3:
        raise RuntimeError("Error: function not removing \
            counterions that it should!")
    else:
        pass

    return

def test_draw_molecule():
    """
    Tests the function that draws a molecule from its SMILES string
        and outputs the drawing to a .png file.
    """
    test_SMILES = 'C[C@@H](CC1=CC=CC=C1)NC'
    test_filename = 'testDrawMolecule.png'
    try:
        context.utilities.draw_molecule(test_SMILES, test_filename)
    except Exception as e:
        assert isinstance(e, TypeError)
    assert os.path.isfile(test_filename), \
        "Function hasn't written file with molecule drawing, or has not \
        done so to the correct directory."
    os.remove(test_filename)

    return

def test_get_l_max():
    """
    Tests the function that returns the wavelength of maximum
        intensity from a list of wavelengths and their intensities.
    """
    test_wl_int = np.array([[400, 0.5], [500, 1]])

    try:
        test_l_max = context.utilities.get_l_max(test_wl_int)
    except Exception as e:
        assert isinstance(e, TypeError)
    assert isinstance(test_l_max, float), \
        'Wrong Type: function should output a float.'

    return

def test_get_em_max():
    """
    Tests the function that retrieves the max values from existing
        files in a list of emission file names.
    """

    return

def test_pad_ndarrays():
    """
    Tests the function that pads out all arrays in a dictionary
        with zeros so that they are all of the same size; the
        size of the largest original input array in the dictionary.
    """
    test_dict = {0:[1, 2, 3], 1:[1, 2, 3, 4, 5]}
    try:
        finger_dict = context.utilities.pad_ndarrays(test_dict)
    except Exception as e:
        assert isinstance(e, TypeError)
    assert isinstance(finger_dict, dict), \
        'Wrong Type: function not returning a dict'
    assert len(finger_dict) == len(test_dict), \
        'Function inserting keys! input and output dictionaries \
        should be of the same size'

    return

def test_compute_fingerprints():
    """
    Tests the function that computes the fingerprints of molecules
        in a dataframe and returns them either in a dictionary or
        to a csv json txt file.
    """
    test_df = pd.DataFrame({'#':[0,1,2],\
    'SMILES':['C1=CC=CC=C1','CC1=CC=CC=C1','CC1=CC=CC=C1C']})
    test_keyname = '#'
    filename = 'fingerprint_test.txt'

    try:
        finger_dict = context.utilities.compute_fingerprints(test_df, \
            key_name=test_keyname)
        finger_dict_file = context.utilities.compute_features(test_df, \
            key_name=test_keyname, output_file=filename)
    except Exception as e:
        assert isinstance(e, TypeError)
    assert isinstance(finger_dict, dict), \
        'Function not outputting a dictionary'
    assert len(finger_dict) == len(test_df), \
        'Function not evaluating fingerprints for all SMILES in input df.'
    assert os.path.isfile(filename), \
        "Function hasn't written property results, or has not done so to the \
        correct directory."
    os.remove(filename)

    return

def test_compute_coulomb_matrixes():
    """
    Tests the function that computes the coulomb matrices of all
        the molecules in an input dataframe.
    """
    test_df = pd.DataFrame({'#':[0,1,2],\
    'SMILES':['C1=CC=CC=C1','CC1=CC=CC=C1','CC1=CC=CC=C1C']})
    test_keyname = '#'
    filename = 'coulomb_test.txt'

    try:
        out_cm = context.utilities.compute_coulumb_matrixes(test_df, \
            key_name=test_keyname)
        out_eig = context.utilities.compute_coulumb_matrixes(test_df, \
            key_name=test_keyname, use_eigval=True)
        out_cm_file = context.utilities.compute_coulumb_matrixes(test_df, \
            key_name=test_keyname, output_file=filename)
    except Exception as e:
        assert isinstance(e, TypeError)
    assert isinstance(out_cm, dict), \
        'Function not outputting a dictionary'
    assert isinstance(out_eig, dict), \
        'Function not outputting a dictionary'
    assert len(out_eig.values()) == len(test_df), \
        'Function not evaluating Coulomb Matrix for all molecules in df'
    assert os.path.isfile(filename), \
        "Function hasn't written property results, or has not done so to the \
        correct directory."
    os.remove(filename)

    return

def test_compute_properties():
    """
    Tests the function that returns the properties of the
        molecules in a dataframe from their SMILES strings.
    """
    test_df = pd.DataFrame({'#':[0,1,2],\
    'SMILES':['C1=CC=CC=C1','CC1=CC=CC=C1','CC1=CC=CC=C1C']})
    test_idxname = '#'
    filename = 'properties_test.txt'

    try:
        out_prop = context.utilities.compute_properties(test_df, \
            index_name=test_idxname)
        out_prop_file = context.utilities.compute_properties(test_df, \
            index_name=test_idxname, output_file=filename)
    except Exception as e:
        assert isinstance(e, TypeError)
    assert isinstance(out_prop, pd.DataFrame), \
        'Function not outputting a pandas dataframe.'
    assert len(out_prop) == len(test_df), \
        'Function not evaluating properties for all molecules in df'
    assert os.path.isfile(filename), \
        "Function hasn't written property results, or has not done so to the \
        correct directory."
    os.remove(filename)

    return

def test_compute_features():
    """
    Tests the function that computes the features of the molecules
        in an input data frame from their SMILES strings.
    """
    test_df = pd.DataFrame({'#':[0,1,2],\
    'SMILES':['C1=CC=CC=C1','CC1=CC=CC=C1','CC1=CC=CC=C1C']})
    test_keyname = '#'
    filename = 'features_test.txt'

    try:
        out_feat = context.utilities.compute_features(test_df, \
            key_name=test_keyname)
        out_feat_file = context.utilities.compute_features(test_df, \
            key_name=test_keyname, output_file=filename)
    except Exception as e:
        assert isinstance(e, TypeError)
    assert isinstance(out_feat, dict), \
        'Function not outputting a dictionary properly.'
    assert os.path.isfile(filename), \
        "Function hasn't written feature results, or has not done so to the \
        correct directory."
    os.remove(filename)

    return

def test_broaden_spectrum():
    """
    Tests the function responsible for broadening the peaks in the
        spectra.
    """
    test_spect = [250, 1]
    sigma = 5.0
    try:
        coords = context.utilities.broaden_spectrum(test_spect, sigma)
    except Exception as e:
        assert isinstance(e, TypeError)
    assert min(coords[0]) == test_spect[0] - 50.0, \
        'wavelength range not correct'
    assert max(coords[0]) == test_spect[0] + 50.0, \
        'wavelength range not correct'
    assert np.isclose(max(coords[1]), 1.0), \
        'maximum normalized intensity not correct.'

    return

def test_visualize():
    """
    Tests the function that plots the predicted spectra and saves
        to a file.
    """
    test_SMILES = 'CCC'
    data_test = [test_SMILES, 200., 1., 250., 1.]

    try:
        context.utilities.visualize(data=data_test, save_dir="/tmp/")
    except Exception as e:
        assert isinstance(e, TypeError)
    base_dir = os.path.dirname(os.path.realpath(__file__))
    visual_file = '/tmp/' + test_SMILES + ".png"
    assert os.path.isfile(visual_file), \
        "Figure hasn't been written to file, or hasn't been saved in \
        the correct directory."
    os.remove(visual_file)

    return
