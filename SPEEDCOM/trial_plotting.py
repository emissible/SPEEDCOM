#Python file to do the trial plotting.  Should be discarded as is unnecessary
#for master.

import pandas as pd
import numpy as np

import utilities

data_file = 'models/trial_abs_ems_data.tsv'

df = pd.read_csv(data_file, sep='\t')

utilities.visualize(df.iloc[0].values)
