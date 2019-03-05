# Use Cases

## Background



## Front End

### GUI
User interactions are handled through our Graphical User Interface (GUI). Here the user inputs a SMILES string for the molecule for which they would like to generate a spectrum. This string is passed to the data-cleaning portion of our program in the back end.

## Back End

## Data Cleaning

The SMILES string for the molecule is received from the front end GUI. After being converted into an approximate 3D geometry, the nuclear coordinates and atomic charges are used to generate a _Coulomb Matrix_, which 

### Neural Network

Because of the complexity and quantity of the descriptors attributed to each molecule, our program involves the utilization of a neural network to build the predictors. This is implemented with [keras](https://github.com/keras-team/keras), an open-source high-level neural-network programming interface. 