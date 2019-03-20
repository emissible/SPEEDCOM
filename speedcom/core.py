import os
import speedcom.utilities as utilities
import speedcom.NNModels as NNModels
import speedcom.data_extract as data_extract

#import speedcom

def initiate_molecules(data_dir):
    """
    Takes in the directory containing all of the data files and returns a list
    of populated molecule objects.
    """
    my_molecules = []
    my_added_molecules = []
    #For all applicable emission spectra put properties into molecule objects,
    #and place objects into list for molecules.
    for ems_file in data_extract.get_emission_files(data_dir):
        pubchem_molecule = data_extract.get_molecule_cid(ems_file) 
        if pubchem_molecule:
            cid = pubchem_molecule.cid
            spectra = data_extract.get_spectra(data_dir + "/" + ems_file)
            absorp = None
            columb = None
            emiss = data_extract.get_peaks(spectra)
            smiles = pubchem_molecule.isomeric_smiles
            weight = pubchem_molecule.molecular_weight
            file_name = ems_file.split(".")[0]
            my_molecules.append(Molecule(absorp, cid, columb, emiss, \
                smiles, weight, file_name))
            my_added_molecules.append(cid)
        else:
           pass
       #For all applicable absorption spectra put properties into molecule 
       #objects, and place objects into list for molecules.
        for abs_file in data_extract.get_absorption_files(data_dir):
            found = 0
            pubchem_molecule = data_extract.get_molecule_cid(abs_file)
            #Check to see if the molecule already exists in the list, and if 
            #it does, does it already have valid absorption peaks:
            if pubchem_molecule:
                cid = pubchem_molecule.cid
                if cid in my_added_molecules:
                    found = 1
                    mol_index = my_added_molecules.index(cid)
                    spectra = data_extract.get_spectra(data_dir + "/" + \
                        abs_file)
                    my_molecules[mol_index].absorption_peaks = \
                        data_extract.get_peaks(spectra)
                    break
                else:
                  pass
                if not found:
                    spectra = data_extract.get_spectra(data_dir + "/" + \
                        abs_file)
                    absorp = data_extract.get_peaks(spectra)
                    columb = None
                    emiss = None
                    smiles = pubchem_molecule.isomeric_smiles
                    weight = pubchem_molecule.molecular_weight
                    file_name = abs_file.split(".")[0]
                    my_molecules.append(Molecule(absorp, cid, \
                        columb, emiss, smiles, weight, file_name))
    return my_molecules

class Molecule:
  """
  The molecule class holds the information about the different molecules:
    * absorption_peaks: The absorption peaks (2D numpy array).
    * cid:              The PUBCHEM CID (float).
    * coulomb:           The columb matrix for the molecule.
    * emission_peaks:   The emission peaks (2D numpy array).
    * smiles:           The canonical SMILES string (string).
    * weight:           The molecular weight (float).
  """
  def __init__(self, absorp, cid, coulomb, emiss, smiles, weight, flnm):
    self.absorption_peaks = absorp
    self.cid = cid
    self.coulomb = coulomb
    self.emission_peaks = emiss
    self.smiles = utilities.remove_cations(smiles)
    self.weight = weight
    self.Descriptors = NNModels.Descriptors(smiles)
    self.file_name = flnm

if __name__ =='__main__':
#    data_dir = os.path.join(os.path.dirname(__file__),'YOUR_SPECTRA_FILE_DIR')
#    molecule_list = initiate_molecules(data_dir)
#    print(len(molecule_list))
    pass
