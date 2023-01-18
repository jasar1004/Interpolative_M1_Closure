# Interpolative first-order maximum entropy (M1) Closure

This directory contains C++ code for performing interpolation of maximum entropy solutions.

The maximum entropy solutions are obtained by solving the associated entropy optimization problem

for sets of moments up to first order.

Description of the different subdirectories:

- M1_Optimization: contains the code for solving the entropy optimization for given sets of angular 
  moments up to first-order. It calls the optimization algorithm from the nlopt open source code 
  for non-linear optimization.

- M1_Closure_interp: contains the code that performs the interpolation of the maximum entropy solutions.
