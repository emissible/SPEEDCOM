# Core.py below

#from models import utilities
from data_clean import data_clean
from data_clean import data_extract

data_dir = "../DATA/PhotochemCAD3/PCAD3 Compd Database 2018"

molecule_list = data_extract.initiate_molecules(data_dir)

print(len(molecule_list))
