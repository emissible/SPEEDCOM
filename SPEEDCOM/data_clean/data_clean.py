"""
Functions for assisting in the cleaning of the obtained information.
"""

import numpy as np

def remove_deliminators(my_strings):
  """
  Remove deliminators from numbers (comma) so as
  to be able to process numbers as int or float types in place of
  strings.  Takes in 'my_array', an array of strings, and will return
  an array of floats.
  """
  my_array = []
  for i in my_strings:
    number = i
    if ',' in i:
      tmp = i.split(",")
      number = tmp[0] + tmp[1]
    try:
      my_array.append(float(number))  
    except:
      print("String" + i + " not able to be cast to float, characters \
        other than ',' or '.'?")
  return np.asarray(my_array)

def remove_catioins(smiles):
  """
  Remove the 1st and 7th row cat/anions from the smiles strings.
  """
  split_smiles = smiles.split(".")
  ion_list = ['[Li+]', '[Na+]', '[K+]', '[Rb+]', '[Cs+]', '[Fr+]', '[F-]', \
    '[Cl-]', '[Br-]', '[I-]', '[At]']
  smiles = [i for i in split_smiles if i not in ion_list]
  smiles = '.'.join(smiles)
  return smiles
