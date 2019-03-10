import numpy as np
import pandas as pd

import json
from NNModels import spDescriptors
#import rdkit
#from rdkit import Chem
#from rdkit.Chem import AllChem
# from rdkit.ForceField.rdForceField import MMFFMolProperties as properties
# import rdkit.Chem.Draw as draw

# def get_largest_num_atoms(train_df):
#     """
#     Returns the number of atoms in the largest molecule of the
#         training data. Intended for determining the size of the
#         Coulomb Matrices.
#
#     Args:
#     -----
#         train_df (pandas.DataFrame) -- the input training data set
#             as retrieved from the PhotoChemCAD database.
#     Returns:
#     --------
#         largest_num_atoms (int)  -- the number of atoms of the
#             largest molecule in the training data set, with
#             Hydrogens included.
#     """
#     # Assertions
#
#     # Generating a list containing the number of atoms in each
#     # molecule in the training DF
#     num_atoms = []
#     for indexi, rowi in df.iterrows():
#         smiles = df.name_smiles[indexi]
#         molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
#         num_atoms.append(molecule.GetNumAtoms())
#     # Adding a new column to the DF that contains the atom count
#     df['num_atoms'] = num_atoms
#     df = df.set_index('#')
#     # Calculating the max number of atoms
#     largest_num_atoms = max(df.num_atoms)
#
#     return largest_num_atoms

def pad_ndarrays(input_dict):
    """
    pad the arrays in the input dictionary as the the maxium length
    among them, with zeros
    """
    lens_of_arrays = []
    for array_i in input_dict.values():
        lens_of_arrays.append(len(array_i))

    max_len = max(lens_of_arrays)
    for array_i_key in input_dict.keys():
        array_i_len = len(input_dict[array_i_key])
        if(array_i_len < max_len):
            input_dict[array_i_key] = \
                np.pad(input_dict[array_i_key], (0, max_len-array_i_len),
                       'constant').tolist()

    return input_dict

def compute_fingerprints(df,SMILES_column='SMILES',key_name=None,radius=2,
                         nBits=2048, use_features=False, padding=True,
                         output_file=None):
    """
    Compute the fingerprints for an input dataframe with all the SIMLES, and
    output the results as an dictionary with json txt format

    Args:
    -----
        df (pandas.DataFrame)  -- an input dataframe with SMILES info
        SMILES_column (str)    -- the column name of SMILES
        key_name (str)         -- the column name for output dict key
        radius (int)           --
        nBits (int)            -- maxium number of bits for fingerprints
                                  computation
        use_features (boolean) -- If True (default as False), use featrues to
                                  compute figureprints
        padding (boolean)      -- If True (default), pad all the fingerprints
                                  to the maxium length in the dictionary
                                  with zeros
        output_file (str)      -- If None, return a dict,
                                  Otherwise output an json txt file.
    Returns:
    --------
        fps_dict (dict)   --  an dictionary contains the fingerprints
                              key    -- name or index of the molecules
                              values -- a list of int

    """
    #Assertions

    #initilized the descriptor engine:
    spD_engine = spDescriptors()

    fps_dict = {}
    for rowi_idx, rowi in df.iterrows():
        spD_engine.set_molecule(rowi[SMILES_column])
        rowi_fp = spD_engine.get_Morgan_fingerprint(radius,nBits,use_features)
        if(key_name is None):
            fps_dict[rowi_idx] = rowi_fp
        else:
            fps_dict[rowi[key_name]] = rowi_fp

    if(padding):
        pad_ndarrays(fps_dict)

    if(output_file is not None):
        with open(output_file, 'w') as f:
            f.write(json.dumps(fps_dict))
    else:
        return fps_dict

def compute_coulumb_matrixes(df,SMILES_column='SMILES', key_name=None,
                             eig_sort=True, padding=True, output_file=None):
    """
    Compute the fingerprints for an input dataframe with all the SIMLES, and
    output the results as an dictionary with json txt format

    Args:
    -----
        df (pandas.DataFrame) -- an input dataframe with SMILES info
        SMILES_column (str)   -- the column name of SMILES
        key_name (str)        -- the column name for output dict key
        eig_sort (boolean)    -- If True (default), sort the coulomb matrixes
                                 with their eigenvalues
        padding (boolean)     -- If True (default), pad all the coulomb matrixes
                                 to the maxium length in the dictionary
                                 with zeros
        output_file (str)     -- If None, return a dict
                                 Otherwise output an json txt file.

    Returns:
    --------
        fps_dict (dict)   --  an dictionary contains the fingerprints
                              key    -- name or index of the molecules
                              values -- a list of list of floats
    """
    #Assertions

    #initilized the descriptor engine:
    spD_engine = spDescriptors()

    CMs_dict = {}
    for rowi_idx, rowi in df.iterrows():
        spD_engine.set_molecule(rowi[SMILES_column])
        rowi_CM = spD_engine.get_coulomb_matrix()
        if(key_name is None):
            CMs_dict[rowi_idx] = rowi_CM
        else:
            CMs_dict[rowi[key_name]] = rowi_CM

    if(padding):
        pad_ndarrays(CMs_dict)

    if(output_file is not None):
        with open(output_file, 'w') as f:
            f.write(json.dumps(CMs_dict))
    else:
        return CMs_dict

def compute_properties(df,SMILES_column='SMILES',index_name=None,
                       output_file=None):
    """
    Compute the fingerprints for an input dataframe with all the SIMLES, and
    output the results as a csv txt file (exported by pandas)

    Args:
    -----
        df (pandas.DataFrame) -- an input dataframe with SMILES info
        SMILES_column (str)   -- the column name of SMILES
        index_name (str)      -- the index name for output DataFrame index
        output_file (str)     -- If None, return a dataframe
                                 Otherwise output an csv txt file.
    Returns:
    --------

    """

def compute_features(df,SMILES_column='SMILES',key_name=None, output_file=None):
    """
    Compute the fingerprints for an input dataframe with all the SIMLES, and
    output the results as a csv txt file (exported by pandas)

    Args:
    -----
        df (pandas.DataFrame) -- an input dataframe with SMILES info
        SMILES_column (str)   -- the column name of SMILES
        key_name (str)        -- the column name for output dict key
        output_file (str)     -- If None, return a dict
                                 Otherwise output an json txt file.
    Returns:
    --------
    """
