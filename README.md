# cardiac_RF_shimming
### Extended RF Shimming Code
Contains code to compute the extended constrained RF shimming optimisation as outlined in a full paper published in 2017 (available [here](https://doi.org/10.1002/nbm.3701)). This repo also includes an example data set to compute RF shims and optimised pulse amplitude and duration that obey all the relvant hardware and safety constraints for an in vivo cardiac experiment.

### Usage
To run this code, it is necessary to first download and install the [CVX version 2.1](http://cvxr.com/cvx/) package from [here](http://cvxr.com/cvx/download/). The standard version contains a conflict when using Hermitian Q-matrices generated in MATLAB so after install, the `mtimes.m` file located in the CVX directory at cvx/builtins/@cvx must be replaced with the version in the `lib/` directory in this repository.

The example script can then be run in MATLAB after CVX install by running `main.m`. All necessary in-vivo and EM simulation data including VOPs and the global Q-matrix are stored in the `bin/` directory.

### Releases and Data

A citable version of release v2.2 is available: [![DOI](https://zenodo.org/badge/41297878.svg)](https://zenodo.org/badge/latestdoi/41297878)

All B<sub>1</sub><sup>+</sup> and SSFP imaging data from the publication is available at: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.231228.svg)](https://doi.org/10.5281/zenodo.231228)
