#ifndef _QUADRULE_HPP_INCLUDED
#define _QUADRULE_HPP_INCLUDED

# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

void cdgqf ( int nt, int kind, long double alpha, long double beta, long double t[], 
  long double wts[] );
void cgqf ( int nt, int kind, long double alpha, long double beta, long double a, long double b, 
  int lo, long double t[], long double wts[] );
void chebyshev_set ( int n, long double x[], long double w[] );
void chebyshev1_compute ( int n, long double xtab[], long double weight[] );
long double chebyshev1_integral ( int expon );
void chebyshev1_set ( int n, long double x[], long double w[] );
void chebyshev2_compute ( int n, long double xtab[], long double weight[] );
long double chebyshev2_integral ( int expon );
void chebyshev2_set ( int n, long double x[], long double w[] );
void chebyshev3_compute ( int n, long double xtab[], long double weight[] );
long double chebyshev3_integral ( int expon );
void chebyshev3_set ( int n, long double x[], long double w[] );
long double class_matrix ( int kind, int m, long double alpha, long double beta, long double aj[], 
  long double bj[] );
void clenshaw_curtis_compute ( int n, long double x[], long double w[] );
void clenshaw_curtis_set ( int n, long double xtab[], long double weight[] );
void fejer1_compute ( int n, long double x[], long double w[] );
void fejer1_set ( int n, long double xtab[], long double weight[] );
void fejer2_compute ( int n, long double x[], long double w[] );
void fejer2_set ( int n, long double xtab[], long double weight[] );
long double gegenbauer_integral ( int expon, long double alpha );
void gegenbauer_ek_compute ( int n, long double alpha, long double a, long double b,
  long double xtab[], long double weight[] );
void gegenbauer_ss_compute ( int n, long double alpha, long double xtab[], 
  long double weight[] );
void gegenbauer_ss_recur ( long double &p2, long double &dp2, long double &p1, long double x,
  int order, long double alpha, long double c[] );
void gegenbauer_ss_root ( long double &x, int order, long double alpha, long double &dp2, 
  long double &p1, long double c[] );
void gen_hermite_dr_compute ( int order, long double alpha, long double x[], long double w[] );
void gen_hermite_ek_compute ( int order, long double alpha, long double x[], long double w[] );
long double gen_hermite_integral ( int expon, long double alpha );
void gen_laguerre_ek_compute ( int order, long double alpha, long double xtab[], 
  long double weight[] );
long double gen_laguerre_integral ( int expon, long double alpha );
void gen_laguerre_ss_compute ( int order, long double alpha, long double xtab[], 
  long double weight[] );
void gen_laguerre_ss_recur ( long double *p2, long double *dp2, long double *p1, long double x, 
  int order, long double alpha, long double b[], long double c[] );
void gen_laguerre_ss_root ( long double *x, int order, long double alpha, long double *dp2, 
  long double *p1, long double b[], long double c[] );
void hermite_ek_compute ( int n, long double x[], long double w[] );
void hermite_gk16_set ( int n, long double x[], long double w[] );
void hermite_gk18_set ( int n, long double x[], long double w[] );
void hermite_gk22_set ( int n, long double x[], long double w[] );
void hermite_gk24_set ( int n, long double x[], long double w[] );
long double hermite_integral ( int n );
long double hermite_integral2 ( long double a );
void hermite_probabilist_set ( int n, long double x[], long double w[] );
void hermite_set ( int order, long double xtab[], long double weight[] );
void hermite_1_set ( int order, long double xtab[], long double weight[] );
void hermite_ss_compute ( int order, long double xtab[], long double weight[] );
void hermite_ss_recur ( long double *p2, long double *dp2, long double *p1, long double x,
  int order );
void hermite_ss_root ( long double *x, int order, long double *dp2, long double *p1 );
int i4_factorial2 ( int n );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_power ( int i, int j );
void imtqlx ( int n, long double d[], long double e[], long double z[] );
void jacobi_ek_compute ( int order, long double alpha, long double beta, long double xtab[], 
  long double weight[] );
long double jacobi_integral ( int expon, long double alpha, long double beta );
void jacobi_ss_compute ( int order, long double alpha, long double beta, long double xtab[], 
  long double weight[] );
void jacobi_ss_root ( long double *x, int order, long double alpha, long double beta, 
  long double *dp2, long double *p1, long double b[], long double c[] );
void jacobi_ss_recur ( long double *p2, long double *dp2, long double *p1, long double x, int order, 
  long double alpha, long double beta, long double b[], long double c[] );
void kronrod_set ( int order, long double xtab[], long double weight[] );
void laguerre_ek_compute ( int order, long double xtab[], long double weight[] );
long double laguerre_integral ( int expon );
void laguerre_set ( int order, long double xtab[], long double weight[] );
void laguerre_1_set ( int order, long double xtab[], long double weight[] );
void laguerre_ss_compute ( int order, long double xtab[], long double weight[] );
void laguerre_ss_recur ( long double *p2, long double *dp2, long double *p1, long double x, 
  int order, long double b[], long double c[] );
void laguerre_ss_root ( long double *x, int order, long double *dp2, 
  long double *p1, long double b[], long double c[] );
long double laguerre_sum ( long double func ( long double x ), long double a, int order, 
  long double xtab[], long double weight[] );
void legendre_dr_compute ( int order, long double xtab[], long double weight[] );
void legendre_ek_compute ( int n, long double x[], long double w[] );
long double legendre_integral ( int expon );
void legendre_recur ( long double *p2, long double *dp2, long double *p1, long double x, 
  int order );
void legendre_set ( int order, long double xtab[], long double weight[] );
void lobatto_compute ( int n, long double x[], long double w[] );
void lobatto_set ( int order, long double xtab[], long double weight[] );
void nc_compute_weights ( int n, long double a, long double b, long double x[], long double w[] );
void ncc_compute ( int order, long double xtab[], long double weight[] );
void ncc_set ( int order, long double xtab[], long double weight[] );
void nco_compute ( int n, long double x[], long double w[] );
void nco_set ( int order, long double xtab[], long double weight[] );
void ncoh_compute ( int n, long double x[], long double w[] );
void ncoh_set ( int order, long double xtab[], long double weight[] );
void parchk ( int kind, int m, long double alpha, long double beta );
void patterson_set ( int order, long double xtab[], long double weight[] );
void psi_values ( int *n_data, long double *x, long double *fx );
long double r8_abs ( long double x );
long double r8_epsilon ( );
long double r8_factorial ( int n );
long double r8_factorial2 ( int n );
long double r8_gamma ( long double x );
long double r8_gamma_log ( long double x );
long double r8_huge ( );
long double r8_hyper_2f1 ( long double a, long double b, long double c, long double x );
long double r8_max ( long double x, long double y );
long double r8_psi ( long double xx );
long double r8_sign ( long double x );
void r8vec_copy ( int n, long double a1[], long double a2[] );
long double r8vec_dot_product ( int n, long double a1[], long double a2[] );
void r8vec_linspace ( int n, long double a_first, long double a_last, long double a[] );
long double *r8vec_linspace_new ( int n, long double a_first, long double a_last );
void r8vec_print ( int n, long double a[], string title );
void r8vec_reverse ( int n, long double x[] );
void radau_compute ( int n, long double x[], long double w[] );
void radau_set ( int order, long double xtab[], long double weight[] );
void scqf ( int nt, long double t[], int mlt[], long double wts[], int nwts, int ndx[], 
  long double swts[], long double st[], int kind, long double alpha, long double beta, long double a, 
  long double b );
void sgqf ( int nt, long double aj[], long double bj[], long double zemu, long double t[], 
  long double wts[] );
void timestamp ( );

#endif //_QUADRULE_HPP_INCLUDED
