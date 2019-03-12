# Core.py below

#from models import utilities
from data_clean import data_clean
#from data_clean import data_extract

class Molecule:
  """
  The molecule class holds the information about the different molecules:
    * absorption_peaks: The absorption peaks (2D numpy array).
    * cid:              The PUBCHEM CID (float).
    * columb:           The columb matrix for the molecule.
    * emission_peaks:   The emission peaks (2D numpy array).
    * smiles:           The canonical SMILES string (string).
    * weight:           The molecular weight (float).
  """
  def __init__(self, absorp, cid, columb, emiss, smiles, weight):
    self.absorption_peaks = absorp
    self.cid = cid
    self.coulomb = columb
    self.emission_peaks = emiss
    self.smiles = data_clean.remove_cations(smiles)
    self.weight = weight
    #self.Descriptors = models.Descriptors(smiles)


