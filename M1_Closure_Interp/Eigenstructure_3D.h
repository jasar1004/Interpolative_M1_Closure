#ifndef _EIGENSTRUCTURE_3D_H_INCLUDED
#define _EIGENSTRUCTURE_3D_H_INCLUDED

// #include "../../Packages/Alglib/src/alglibinternal.h"
// #include "../../Packages/Alglib/src/linalg.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <limits>
#include <mpi.h>
#include <unistd.h>
#include <complex>
#include <cmath>

#ifndef _M1_Model_1D_H_INCLUDED
#include "../M1_Optimization/M1_Optimization.h"
#endif // _M1_Model_1D_H_INCLUDED

using namespace std;
// using namespace alglib;
// using namespace alglib_impl;

#define TOLER_HYPER         1.0e-4

long double Chi_2_Lower_bound_1D(const long double &N1);

long double Chi_2_Upper_bound_1D(const long double &N1);

long double xi(const long double &N1);

long double Chi_2_Gray(const long double &N1);

void Check_Realizability(const long double &N1, const long double &Chi2);

long double unit_vec_N1(const int &i, const long double &N1, const long double &phi);

long double I2ij(const int &i, const int &j, const long double &I0, const long double &N1, const long double &phi, const long double &Chi2);

long double dI2ij_dI0(const int &i, const int &j, const long double &r_E, const long double &N1, const long double &phi, const long double &Chi2, const long double &E_dChi_dE);
  
long double dI2ij_dIi(const int &i, const int &j, const int &index_N1, const long double &r_E, const long double &N1, const long double &phi, const long double &Chi2, const long double &dChi_dN1);

void Cubic_solver(const long double &aa, const long double &bb, const long double &cc, const long double &dd, long double &x1, long double &x2, long double &x3, bool flag_Imaginary_Part = false);

void Evaluate_eigenvalues_1D(const long double &r, const long double &N1, const long double &Chi2, const long double &E_dChi2_dE, const long double &dChi2_dN1, long double &lambda_1, long double &lambda_2);

void Evaluate_Chi2_derivatives_Gray(const long double &r, const long double &N1, long double &E_dChi2_dE, long double &dChi2_df);

long double dChi2_dN1_Gray(const long double &N1);

long double d2Chi2_dN1_2_Gray(const long double &N1);

long double d3Chi2_dN1_3_Gray(const long double &N1);

void Test_Finite_Difference();

void Eigenvalues_Gray_M1_Closure_2D(long double &lambda_1, long double &lambda_2, long double &lambda_3, const long double &N1, const long double &phi);

void Evaluate_Lambdas_x(long double &lambda_1, long double &lambda_2, long double &lambda_3, const long double &r_E, const long double &N1, const long double &phi, const long double &Chi2, const long double &E_dChi2_dE, const long double &dChi2_df, bool flag_Imaginary_Part = false);

void Evaluate_eigenvalues_2D(const long double &r, const long double &N1, const long double &phi, const long double &Chi2, const long double &E_dChi2_dE, const long double &dChi2_dN1, long double &lambda_1, long double &lambda_2, long double &lambda_3);

void Setup_lc(long double lc_vec[][3], const int &index_U, const long double &r_E, const long double &N1, const long double &phi, const long double &Chi2, const long double &E_dChi2_dE, const long double &dChi2_df);

void Setup_rc(long double rc_vec[][3], const int &index_U, const long double &r_E, const long double &N1, const long double &phi, const long double &Chi2, const long double &E_dChi2_dE, const long double &dChi2_df);

void Eigenstructure_Check(const long double &r_E, const long double &N1, const long double &phi, const long double &Chi2, const long double &E_dChi2_dE, const long double &dChi2_df, const int &Problem_Type);

#endif // _EIGENSTRUCTURE_3D_H_INCLUDED
