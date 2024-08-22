Tidy, Scientific Full-Factorial Numerical Simulations in Matlab
===============================================================

[![DOI](https://zenodo.org/badge/105531812.svg)](https://zenodo.org/badge/latestdoi/105531812)

SimRunner enables you to easily run any numerical simulation main file over a 
full-factorial design range of paramter choices, including numeric and 
cell-array of string values. Intelligent estimates of end-time are computed 
from completed experiments using ensemble regression trees. Final results are 
packaged and saved sequentially in case of any problem.

To aid scientific computing, random seeds are controlled by default, and the 
final results package contains the specific parameters used for each and every 
experiment.

Support is also provided for running experiments in parallel (requiring the 
parallel computing toolbox).

To get started, type 'help SimRunner' at the MATLAB prompt. A test main file 
and runfile are provided.
