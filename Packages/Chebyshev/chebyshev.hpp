#ifndef _CHEBYSHEV_HPP_INCLUDED
#define _CHEBYSHEV_HPP_INCLUDED

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

#define UNIFORM_DISTRIBUTION                  10000
#define CHEBYSHEV_FIRST_KIND_DISTRIBUTION     10001
#define CHEBYSHEV_SECOND_KIND_DISTRIBUTION    10002
#define LOBATTO_DISTRIBUTION                  10003

#define SPHERICAL_HARMONIC_DISTRIBUTION       10004

#define PI_CHEBY 3.141592653589793238462643383279502884

using namespace std;

long double Chebyshev_Polynomial_Basis(const long double &x, const int &Index);
long double Chebyshev_Second_Kind_Polynomial_Basis(const long double &x, const int &N);
long double chebyshev_points_first_kind ( int index_node, int n );
long double chebyshev_points_second_kind ( int index_node, int n );
long double zeros_shifted ( int index_node, int n, long double a, long double b, int node_distribution);
void chebyshev_quadrature ( long double &weight, long double &x, int index_node, int n, long double a, long double b, int node_distribution);
void chebyshev_quadrature ( long double *weight, long double *x, int n, long double a, long double b, int node_distribution);
long double Uniform_Distribution(const int &i, const int &Np, const long double &val_min, const long double &val_max);
long double Uniform_Distribution_No_Endpoint(const int &i, const int &Np, const long double &val_min, const long double &val_max);

#endif // _CHEBYSHEV_HPP_INCLUDED
