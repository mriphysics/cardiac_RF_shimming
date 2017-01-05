# cardiac_RF_shimming
### Constrained RF Shimming Code
Cotains code to compute SAR constrained RF shimming optimisation as outlined in a future publication. Currently this repo includes an example data set to compute shims for an in vivo cardiac experiment for a fixed RF pulse amplitude. Soon the full computation for variable RF pulse amplitudes will be added.

### Usage
To run this code, it is necessary to first download and install the [CVX version 2.1](http://cvxr.com/cvx/) package from [here](http://cvxr.com/cvx/download/). The standard version contains a conflict when using Hermitian Q-matrices generated in MATLAB so after install, the `mtimes.m` file located in the CVX directory at cvx/builtins/@cvx must be replaced with the version in the `lib/` directory in this repository.

The example script can then be run in MATLAB after CVX install by running `main.m`. All necessary in-vivo and EM simulation data including VOPs and the global Q-matrix are stored in the `bin/` directory.

A citable version of release v2.2 is available: [![DOI](https://zenodo.org/badge/41297878.svg)](https://zenodo.org/badge/latestdoi/41297878)

All B<sub>1</sub><sup>+</sup> and SSFP imaging data from the publication is available at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.231228.svg)](https://doi.org/10.5281/zenodo.231228)
