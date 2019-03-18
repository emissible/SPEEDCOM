[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.com/emissible/SPEEDCOM.svg?branch=master)](https://travis-ci.com/emissible/SPEEDCOM)
[![Coverage Status](https://coveralls.io/repos/github/emissible/SPEEDCOM/badge.svg?branch=master)](https://coveralls.io/github/emissible/SPEEDCOM?branch=master)

<p align="center"><img src="doc/source/_static/logo.png" alt="SPEEDCOM" title="SPEEDCOM"/></p>

Authors: **Joe Abbott**, **Ryan Beck**, **Hang Hu**, **Yang Liu**, **Lixin Lu**.

<span style="color:red"> _**SPEEDCOM** is currently under development. We should be up and running soon._ </span>

## Overview

_SPEEDCOM_ is an open source python package that aims to predict the fluorescence emission and absorption spectra of small conjugated organic molecules. These features are predicted using a neural network, implemented with [keras](https://github.com/keras-team/keras), and are trained on data from the [PhotochemCAD database](http://www.photochemcad.com/PhotochemCAD.html). The software has a graphical-user-interface (GUI) where users can input the [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) string for a given molecule and be returned its predicted spectra and associated characteristic quantities. For further details on the background science, and the operations of our program, please see our [use cases](https://github.com/emissible/SPEEDCOM/blob/master/use_cases.md).


## Configuration

### Pre-requirements:

* Python version 3.6.7 or later
* conda version 4.6.8 or later

### Installation
From your computer's terminal application, enter the following command:

``conda install speedcom``

## Contributions

Any contributions to the project are warmly welcomed! If you discover any bugs, please report them in the [issues section](https://github.com/emissible/SPEEDCOM/issues) of this repository and we'll work to sort them out as soon as possible. If you have data that you think will be good to train our model on, please contact one of the authors. 


## License

SPEEDCOM is licensed under the [MIT license](https://github.com/emissible/SPEEDCOM/blob/master/LICENSE).
