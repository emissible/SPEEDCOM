from itertools import zip_longest
import math
import numpy as np
import os
import pandas as pd

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import Properties
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig

class Descriptors:
    """
    A wrapper class using rdkit to generate the different
        descpritors for specrtra prediciton.

    Initilized with SMILES of a molecule

    Attributes:
    -----------
        Molecule       -- an object of rdkit.Chem.rdchem.Mol
        __feat_factory -- an object to caluate the features
            for molecules, from
            rdkit.Chem.rdMolChemicalFeatures.MolChemicalFeatureFactory

    Methods:
    --------
        set_molecule
        get_features
        get_properties
        get_coulomb_matrix
        get_Morgan_fingerprint
        __config_feature_factory
        __get_charges_coords

    """

    def __init__(self, SMILES = None):
        """
        Descriptor class constructor method.
        """
        if(SMILES is not None):
            self.set_molecule(SMILES)
        else:
            self.Molecule = None

        #functions extend from rdkit package
        self.__feat_factory  = None
#        self.Features    = None
#        self.Fingerprint = None


    def set_molecule(self, SMILES):
        """
        Assigns a molecule to the initialized Descriptor class.

        Args:
        -----
            SMILES (str) -- the string representation of a molecule
        """
        # Assertions
        assert isinstance(SMILES, str),\
            'the SMILES must be a string'
        # Functionality
        self.Molecule = Chem.MolFromSmiles(SMILES)

        return

    def get_properties(self, property_name=None):
        """
        Returns the properties of a molecule object using
            the rdkit package.

        Args:
        -----
            property_name (str) -- the name of a property to be
                returned. If None, function returns all properties.

        Returns:
        --------
            f_dict2 (dict) -- a dictionary containing all the
                molecule's properties as keys and their float
                values as the dict values.
        """
        # Assertions
        assert type(self.Molecule) == Chem.rdchem.Mol
        # Functionality
        f_dict = dict(zip(Properties().GetPropertyNames(),\
                    Properties().ComputeProperties(self.Molecule)))

        if(property_name is None):
            return f_dict
        else:
            if type(property_name) == str:
                return f_dict[property_name]
            elif type(property_name) == list:
                f_dict2 = {}
                for i in range(len(property_name)):
                    f_dict2[property_name[i]] = f_dict[property_name[i]]
                return f_dict2

    def get_features(self):
        """
        Gets the features of the molecule currently constructed in
            the Descriptor class.

        Args:
        -----
        Returns:
        --------
            features_dict (dict) -- a dictionary conatining the
                feature names as keys and their values as the
                values of the dictionary.
        """
        if(self.__feat_factory is None):
            self.__config_feature_factory()

        assert type(self.Molecule) == Chem.rdchem.Mol

        features = self.__feat_factory.GetFeaturesForMol(self.Molecule)
        features_dict = {}
        for i in range(len(features)):
            f_type = features[i].GetType()
            f_ids  = features[i].GetAtomIds()
            if (f_type in features_dict.keys()):
                features_dict[f_type].append(f_ids)
            else:
                features_dict[f_type] = [f_ids]

        return features_dict

    def get_Morgan_fingerprint(self, radius=2, nBits=2048, use_feat=False):
        """
        Returns a list of integers that as a whole represents the
            Morgan fingerprint of a molecule.

        Args:
        -----
            radius (int)    -- the radius parameter to be passed to
                the
            nBits (int)     -- the number of bits the molecule
                is characterized by.
            use_feat (bool) -- (default False) If True, uses the
                molecule's features when calculating the finger-
                print using RDKit.

        Returns:
        --------
            out_fp (list) -- a list containing the bits that
                comprise the Morgan fingerprint of the molecule.
        """
        # Assertions
        assert type(self.Molecule) == Chem.rdchem.Mol
        # Functionality
        fp = AllChem.GetMorganFingerprintAsBitVect(self.Molecule, radius,
                                                    nBits=nBits,
                                                    useFeatures=use_feat)
        out_fp = list(fp.ToBinary())

        return out_fp

    def get_coulomb_matrix(self, eig_sort=True, output_eigval=False):
        """
        Generates the coulomb matrix for a given molecule from its
            SMILES string, of size MxM, where M is the number of
            atoms in the molecule. in the training data set.

        Args:
        -----
            eig_sort (bool) -- (Default True) If False, doesn't sort
                the coulomb matrix by eignevalue.
            output_eigval (bool) -- (Default False) If True, outputs
                the eigenvalues of the Coulomb Matrix instead of the
                Coulomb matrix itself.

        Returns:
        --------
            coulomb_matrix (numpy.ndarray) -- the coulomb matrix for
                a given molecule's nuclear geometry in the form of a
                2D aray.

            eig (np.array) -- Outputs if the output_eigval argument
                is True as a 1D array of the eigenvalues of the Coulomb
                matrix.
        """
        # Assertions
        assert type(self.Molecule) == Chem.rdchem.Mol
        assert isinstance(eig_sort, bool), 'eig_sort must be a bool'
        assert isinstance(output_eigval, bool), 'output_eigval must be a bool'

        # Generating the coulomb matrix
        molecule_df = self.__get_charges_coords()
        num_atoms = len(molecule_df)
        coulomb_matrix = np.zeros(shape=(num_atoms,num_atoms))
        for indexi, rowi in molecule_df.iterrows():
            for indexj, rowj in molecule_df.iterrows():
                Zi = rowi.charge
                xi = rowi.x
                yi = rowi.y
                zi = rowi.z
                Zj = rowj.charge
                xj = rowj.x
                yj = rowj.y
                zj = rowj.z
                if indexi == indexj:
                    element = 0.5 * math.pow(Zi, 2.4)
                else:
                    norm_diff = math.sqrt(math.pow((xi-xj),2) + math.pow((yi-yj),2) + math.pow((zi-zj),2))
                    element = Zi * Zj / norm_diff
                coulomb_matrix[indexi][indexj] = element

        if output_eigval:
            eig = np.linalg.eig(coulomb_matrix)[0]
            eig.sort()

            return eig

        if eig_sort:
            eig = np.linalg.eig(coulomb_matrix)[0]
            eig_idx = np.argsort(eig)
            temp_matrix = np.zeros(shape=(num_atoms,num_atoms))
            for i in range(num_atoms):
                temp_matrix[i] = coulomb_matrix[eig_idx[i]]
            temp_matrix = temp_matrix.transpose()
            for i in range(num_atoms):
                coulomb_matrix[i] = temp_matrix[eig_idx[i]]

        return coulomb_matrix

    def __config_feature_factory(self):
        """
        Initialize the 'feature factory' rdkit module with the
            current molecule.
        """
        fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
        self.__feat_factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

        return

    def __get_charges_coords(self):
        """
        Generates a pandas dataframe containing the charges and cartesian
            coordinates of each atom in a molecule.

        Args:
        -----
            SMILES (str) -- the SMILES string representation of the
                molecule.

        Returns:
        --------
            molecule_df (pandas.DataFrame) -- contains the charge and
                cartesian coordinate information for each atom within
                a molecule.
        """
        # Assertions

        # Building the benzene molecule and ADDING HYDROGENS
        molecule = Chem.AddHs(self.Molecule)
        # 'Embedding' the molecular coordinates, optimising structure
        AllChem.EmbedMolecule(molecule)
        AllChem.MMFFOptimizeMolecule(molecule)
        # Generating universal force field model
        ff = AllChem.UFFGetMoleculeForceField(molecule)
        # Getting the positions of nuclei; returned as a tuple
        # of the form (x1, y1, z1, x2, y2, z2, x3, ...)
        positions = ff.Positions()

        # Creating a list of the atomic numbers
        atomic_nums = []
        for atom in molecule.GetAtoms():
            atomic_nums.append(atom.GetAtomicNum())

        # Creating lists of the cartesian coordinates of the atoms
        x = []
        y = []
        z = []
        for item1, item2, item3 in self.__grouper(3, positions):
            x.append(item1)
            y.append(item2)
            z.append(item3)

        # Building a DF with predictors
        molecule_df = pd.DataFrame()
        molecule_df['charge'] = atomic_nums
        molecule_df['x'] = x
        molecule_df['y'] = y
        molecule_df['z'] = z

        return molecule_df

    def __grouper(self, n, iterable, fillvalue=None):
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
