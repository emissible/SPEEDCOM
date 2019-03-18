import numpy as np
import os
import pubchempy as pcp
import scipy.signal

# import the speedcom.Molecule class
#import speedcom

"""
The data obtaining functions.  Please note that many of these functions
should only need to be run when training a new model.
"""

#data_dir global variable containing the path to the data
def get_emission_files(data_dir):
    """
    Get the total files containing emission spectra.  Expects the files to 
    be of type *.ems.txt and have the naming convention *.ems.txt will return
    a list of strings to be parsed later.
     """
    return [fn for fn in os.listdir(data_dir) if fn.endswith(".ems.txt")]

def get_absorption_files(data_dir):
    """
    Get the total files containing absorption spectra.  Expects the files to
    be of type *.txt and have the naming convention *.abs.txt will return a 
    list of strings to be parsed later.
    """
    return [fn for fn in os.listdir(data_dir) if fn.endswith(".abs.txt")]

def get_molecule_cid(file_name):
    """ 
    Take a file_name, remove the ending and prefix and return the pubchem
    molecule object.  If the file_name does not have the cas number, try the 
    name of the system,   else return None so as to disreguard the molecule
    from training.
    """
    my_file = file_name.split("_")
    #try to get the cas # if not then try name.
    cas = my_file[1]
    try:
        cid = pcp.get_compounds(cas, 'name')[0]
    except:
        #file doesn't have a cas number with it, try name
        cas = my_file[-1]
        cas.replace("_", ' ')
        my_name = cas.split(".")
        try:
            cid = pcp.get_compounds(my_name[0], 'name')[0]
        except:
            cid = None
    finally:
        return cid

def get_spectra(file_name):
    """
    Iterates through the file and pulls out the normalized spectrum contained
    within.  The X is in wavelentgh, and the intensity is normalized to the
    most intense peak.
    """
    spectra = []
    with open (file_name,'r') as file:
        next(file)
        for line in file:
            tmp = line.split()
            if len(tmp) > 1:
                spectra.append([float(tmp[0]), float(tmp[1])])
            else:
                continue
    return np.asarray(spectra)

def get_peaks(spectra):
    peaks = []
    """
    This will take in the spectra from the get_spectra function and will return
    a list of [x,Y] coordinates for the peaks.
    """
    peak_locations = scipy.signal.find_peaks(spectra[:,1])[0]
    for i in peak_locations:
        peaks.append(spectra[i, :])
    return np.asarray(peaks)

def electrostatic_potentials(file_name):
    """
    Returns the dielectic constant for the various solvents within the
    database.  If a solvent is not in the list will return 1. Takes a 
    string and returns a double
  
    Sources: 
    acetic acid to water (pH 7): 
        https://www.organicdivision.org/wp-content/uploads/2016/12/\
            organic_solvents.html
    Glycerol - Sulfuric Acid: 
        http://www.appliedmc.com/content/images/Dielectric_Constants.pdf
    Chloronaphthalene - http://www.stenutz.eu/chem/solv6.php?\
        name=1-chloronaphthalene
    3-methylpentante - https://journals.aps.org/pra/pdf/10.1103/\
        PhysRevA.42.4735
    PBS - http://www.materials-talks.com/blog/2015/03/09/how-to-get-the\
        -latest-software-for-the-zetasizer-nano/
    Borate Buffer - https://www.honeywellprocess.com/library/marketing/\
        tech-specs/Dielectric%20Constant%20Table.pdf
    Sodium Hydroxide - http://www.microkat.gr/msdspd90-99/Sodium%20\
        hydroxide.html
    """
    return{
        "none":                                     1.0,
        "acetic acid":                              6.20,
        "acetonitrile":                             36.64,
        "benzene":                                  2.28,
        "chlorobenzene":                            5.69,
        "chloroform":                               4.81,
        "cyclohexane":                              2.02,
        "dichloromethane":                          9.08,
        "diethyl ether":                            4.267,
        "dimethylformamide":                        38.25,
        "dioxane":                                  2.21,
        "DMSO":                                     47.0,
        "30% tris buffered (in DMSO)":              47.0,
        "ethanol":                                  24.6, 
        "ethanol ":                                 24.6,
        "ethanol(basic)":                           24.6,
        "ethanol (basic)":                          24.6,
        "hexane":                                   1.89,
        "isopropanol":                              18.3,
        "methanol":                                 32.6,
        "propanol":                                 20.1,
        "pyridine":                                 12.3,
        "tetrahydrofuran":                          7.52,
        "toluene":                                  2.38,
        "water":                                    78.54,
        "water (pH 5 to 9)":                        78.54,
        "water (pH 5 to 9 )":                       78.54,
        "water (pH 4)":                             78.54,
        "water (pH 7)":                             78.54,
        "water (pH 10)":                            78.54,
        "glycerol":                                 13.2,
        "0.5 M H2SO4":                              31.0,
        "H2SO4 aq (1 N)":                           31.0,
        "chloronaphthalene":                        5.0,
        "3-methylpentane":                          9.7,
        "PBS":                                      79.0, 
        "PBS ":                                     79.0,
        "phosphate buffer (pH 7, 0.1 M)":           79.0,
        "borate buffer (pH 10)":                    8.2,
        "NaOH aq":                                  57.5,
        "0.1 N NaOH aq":                            57.5
        }.get(file_name, 1.0)
