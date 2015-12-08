# cardiac_RF_shimming
### Constrained RF Shimming Code
Cotains code to compute SAR constrained RF shimming optimisation as outlined in a future publication. Currently this repo includes an example data set to compute shims for an in vivo cardiac experiment for a fixed RF pulse amplitude. Soon the full computation for variable RF pulse amplitudes will be added.

### Usage
To run this code, it is necessary to first download and install the [CVX version 2.1](http://cvxr.com/cvx/) package from [here](http://cvxr.com/cvx/download/). The standard version contains a conflict when using Hermitian Q-matrices generated in MATLAB so after install, the `mtimes.m` file located in the CVX directory at cvx/builtins/@cvx must be replaced with the version in this repository.
