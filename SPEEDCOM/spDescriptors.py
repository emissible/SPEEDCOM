import os
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import Properties
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig

class spDescriptors:
    """
    A Class to generate the descpritors for specrtra prediciton
    """

    def __init__(self, SMILEs = None):
        """
        spDescriptor Constructor
        """
        if(SMILEs is not None):
            self.set_molecule(SMILEs)
        else:
            self.Molecule = None

        #functions extend from rdkit package
        self.feat_factory  = None
#        self.Features    = None
#        self.Fingerprint = None


    def set_molecule(self, SMILEs):
        """ set molecule of the spDecriptor"""
        self.Molecule = Chem.MolFromSmiles(SMILEs)
        return

    def get_properties(self, feature_name = None):
        """  """

        assert type(self.Molecule) == Chem.rdchem.Mol
        f_dict = dict(zip(Properties().GetPropertyNames(),\
                    Properties().ComputeProperties(self.Molecule)))

        if(feature_name is None):
            return f_dict
        else:
            if type(feature_name) == str:
                return f_dict[feature_name]
            elif type(feature_name) == list:
                f_dict2 = {}
                for i in range(len(feature_name)):
                    f_dict2[feature_name[i]] = f_dict[feature_name[i]]
                return f_dict2

    def get_features(self):
        if(self.feat_factory is None):
            self.config_feature_factory()

        assert type(self.Molecule) == Chem.rdchem.Mol

        features = self.feat_factory.GetFeaturesForMol(self.Molecule)
        features_dict = {}
        for i in range(len(features)):
            f_type = features[i].GetType()
            f_ids  = features[i].GetAtomIds()
            if (f_type in features_dict.keys()):
                features_dict[f_type].append(f_ids)
            else:
                features_dict[f_type] = [f_ids]

        return features_dict

    def get_Morgan_fingerprint(self, radius=2, nBits=2048, use_features=False):
        """  """
        assert type(self.Molecule) == Chem.rdchem.Mol
        fp = AllChem.GetMorganFingerprintAsBitVect(self.Molecule, radius, nBits=nBits,
                                                   useFeatures=use_features)
        return list(fp.ToBinary())

    def get_coulomb_matrix(self):
        """ """
        

    def config_feature_factory(self):
        fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
        self.feat_factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
        return
