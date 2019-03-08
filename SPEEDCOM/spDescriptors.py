import os
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import Properties
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig

class spDescriptor:
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

    def get_properties(self):
        """  """

        assert type(self.Molecule) == Chem.rdchem.Mol
        return dict(zip(Properties().GetPropertyNames(),\
                    Properties().ComputeProperties(self.Molecule)))

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
        fp = AllChem.GetMorganFingerprintAsBitVect(self.Molecule, raduis, nBits=nBits,
                                                   useFeatures=use_features)
        return list(fp)

    def config_feature_factory(self):
        fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
        self.feat_factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
        return
