from NNModels import Descriptors
import pandas as pd

def get_properties_df(df):
    """
    Computes and returns a data frame containing all the properties
        of each molecule in an input data frame, as returned by the
        get_properties function in the NNModels class.

    Args:
    -----
        df (pandas.DataFrame) -- the input data frame

    Returns:
    --------
        properties_df(pandas.DataFrame) -- the output data frame
            containing all property descriptors for all molecules.

    """
    # Assertions
    assert isinstance(df, pd.DataFrame), 'Input not a pandas DF'
    # Generating the output DataFrame
    properties = dict.fromkeys(df.index)

    for indexi, rowi in df.iterrows():
        molecule = Descriptors(rowi.name_smiles)
        properties[indexi] = molecule.get_properties()

    properties_df = pd.DataFrame(properties).transpose()
    properties_df[['Wavelength', 'Epsilon', 'Quantum Yield']] \
        = df[['Wavelength', 'Epsilon', 'Quantum Yield']]

    return properties_df
