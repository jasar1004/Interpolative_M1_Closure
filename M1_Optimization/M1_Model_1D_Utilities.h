#ifndef _M1_Model_1D_UTILITIES_H_INCLUDED
#define _M1_Model_1D_UTILITIES_H_INCLUDED

/*******************************************************************
  File: M1_Model_1D_Utilities.h

  Description:  ...  

  Author:  Joachim A.R. Sarr

  Date:    May 05th, 2020
*******************************************************************/

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

#include "M1_Optimization.h"

#define DIFF_E      1000
#define DIFF_N1     1001

using namespace std;

long double binomialCoefficients(const long double &n, const long double &k);

long double forward_finite_difference(const long double *u_y, const int &index_e, const int &index_f, const int &N_points_E, const int &N_points_f, const int &N_points_Theta, const int &N_points_Phi, const long double &h, const int &order, const int &VAR_NUM);
long double backward_finite_difference(const long double *u_y, const int &index_e, const int &index_f, const int &N_points_E, const int &N_points_f, const int &N_points_Theta, const int &N_points_Phi, const long double &h, const int &order, const int &VAR_NUM);
long double central_finite_difference(const long double *u_y, const int &index_e, const int &index_f, const int &N_points_E, const int &N_points_f, const int &N_points_Theta, const int &N_points_Phi, const long double &h, const int &order);
long double Lebedev_Quadrature_Orthog_SH(const int &index_SH_1, const int &index_SH_2, const int &order);
long double Lebedev_Quadrature_Matrix_SH_Temp(const long double *Matrix, const int &index_SH, const int &degree_SH, const int &order);
void Orthogonality_Test(const int &order);
int grid_E(const long double &r_E, const long double *rE_uniform, const int &Npts_E, const int &Npts_f, const int &index_N1);

void LU_Solve_Back(const int &n, long double *x, const long double *A, const long double *b);
void LU_Solve_Forward(const int &n, long double *x, const long double *A, const long double *b);
void LUdecomposition(const long double *A, long double *L, long double *U, const int &n);
void Solve_A_x_b(const long double *VanderMonde_Matrix, long double *Coeff_Vand_Matrix, const long double *VanderMonde_Vector, const int &n, const int &index_Cheby);
void Check_A_x_b(const long double *VanderMonde_Matrix, const long double *Coeff_Vand_Matrix, const int &n);
int SH_Linear_Index(const int &Order_SH, const int &degree_SH);
long double Precompute_SH_Basis_Coeffs(const long double &l, const long double &m, const long double &i, const long double &j, const long double &k);
long double Cartesian_Spherical_harmonics(const long double &x, const long double &y, const long double &z, const int &l, const int &m);
long double Test_Spherical_Harmonics(const long double *Vals, const int &index_Cheby, const int &N_points_Cheby, const long double &x, const long double &y, const long double &z, const int &Order_SH);
void Duffy_Transformation(long double &zeta_p, long double &eta_p, const long double &zeta, const long double &eta);
long double VanderMonde_Matrix_1_var(const long double &var_1, const int &Order_poly_1);
long double VanderMonde_Matrix_2_vars(const long double &var_1, const long double &var_2, const int &Order_poly_1, const int &Order_poly_2);
long double VanderMonde_Vector_N_vars(const int &Index_Entry, const int &Index_Point);
void Vandermonde_Interpolation_1var(long double *Coefficient_Matrix_Fit_In_plus, const long double *f_In_plus_NON_GRAY_Cheby, const int &N_Coeffs_x, const int &Lebedev_Order, const long double *x_NON_GRAY_Cheby, const int &var_num);

long double Precompute_First_Kind_Chebyshev_Basis_Coeffs(const int &n, const int &k);

void Chebyshev_First_Kind_to_Monomial_Basis_ratio_E_No_SH(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_mu);

void Chebyshev_First_Kind_to_Monomial_Basis_Norm_f_No_SH(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_mu);

void Chebyshev_First_Kind_to_Monomial_Basis_mu_No_SH(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_mu);

void Chebyshev_First_Kind_to_Monomial_Basis_ratio_E(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_SH);

void Chebyshev_First_Kind_to_Monomial_Basis_Norm_f(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_SH);

void Spherical_Harmonics_to_Monomial_Basis(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_SH, const int &Order_SH, const int *Array_l_SH, const int *Array_m_SH);

long double Lebedev_Quadrature_Orthog_SH(const int &order_SH1, const int &order_SH2, const int &m1, const int &m2, const int &order);

void Orthogonality_Test(const int &Quad_Rule);

#endif // _M1_Model_1D_UTILITIES_H_INCLUDED
