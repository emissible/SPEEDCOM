import math
import numpy as np
import pandas as pd
from itertools import zip_longest

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
# from rdkit.ForceField.rdForceField import MMFFMolProperties as properties
# import rdkit.Chem.Draw as draw


def get_largest_num_atoms(train_df):
    """
    Returns the number of atoms in the largest molecule of the
        training data. Intended for determining the size of the
        Coulomb Matrices.

    Args:
    -----
        train_df (pandas.DataFrame) -- the input training data set
            as retrieved from the PhotoChemCAD database.
    Returns:
    --------
        largest_num_atoms (int)  -- the number of atoms of the
            largest molecule in the training data set, with
            Hydrogens included.
    """
    # Assertions

    # Generating a list containing the number of atoms in each
    # molecule in the training DF
    num_atoms = []
    for indexi, rowi in df.iterrows():
        smiles = df.name_smiles[indexi]
        molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
        num_atoms.append(molecule.GetNumAtoms())
    # Adding a new column to the DF that contains the atom count
    df['num_atoms'] = num_atoms
    df = df.set_index('#')
    # Calculating the max number of atoms
    largest_num_atoms = max(df.num_atoms)

    return largest_num_atoms

def grouper(n, iterable, fillvalue=None):
    """
    A function to aggregates items in a list into groups of
        assigned size.

    Args:
    -----
        n (int)         -- the number of items per group.
        iterable (list) -- the list to iterate over, whose items will
            be grouped.
        fillvalue       -- the value to fill in empty compartments of a
            group if the number of items in the iterator isn't a multiple
            of the group size (n).

    Returns:
    --------
        grouper (itertools.zip_longest) -- the iterator for the
            aggregation of items in the input list.
    """
    args = [iter(iterable)] * n
    grouper = zip_longest(fillvalue=fillvalue, *args)

    return grouper

def get_padded_coulomb_matrix(coulomb_matrix, M):
    """
    Takes an existing NxN coulomb matrix for a molecule with N number
        of atoms and, if N is less than the max size of the matrix M,
        pads out the coulomb matrix with zeros so that it is of size
        MxM.

    Args:
    -----
        coulomb_matirx (numpy.ndarray) -- the coulomb matrix for a
            given molecule.
        M (int) -- the side length of the matrix to attain an MxM
            matrix.

    Returns:
    --------
        padded_coulomb (numpy.ndarray) -- an MxM coulomb matrix.

    """
    N = len(coulomb_matrix)
    padded = np.zeros(shape=(M,M))
    padded[:N,:N] = coulomb

    return padded_coulomb
