#ifndef _CHEBYSHEV_HPP_INCLUDED
# include "chebyshev.hpp"
#endif // _CHEBYSHEV_HPP_INCLUDED

//****************************************************************************
// This routine computes the Chebyshev polynomial of the first kind of 
// degree N
//****************************************************************************
long double Chebyshev_Polynomial_Basis(const long double &x, const int &N) {
    long double Cheby, Cheby_n_m_1, Cheby_n;
    Cheby_n_m_1 = 1.0;
    Cheby_n = x;
    
    if (N == 0) {
        Cheby = Cheby_n_m_1;
    } else if (N == 1) {
        Cheby = Cheby_n;
    } else {
        for (int i = 0; i < N-1; i++) {
            Cheby = 2.0*x*Cheby_n - Cheby_n_m_1;
            Cheby_n_m_1 = Cheby_n;
            Cheby_n = Cheby;
        }
    }
    
    return Cheby;
}

//****************************************************************************
// This routine computes the Chebyshev polynomial of the second kind of 
// degree N
//****************************************************************************
long double Chebyshev_Second_Kind_Polynomial_Basis(const long double &x, const int &N) {
    long double Cheby, Cheby_n_m_1, Cheby_n;
    Cheby_n_m_1 = 1.0;
    Cheby_n = 2.0*x;
    
    if (N == 0) {
        Cheby = Cheby_n_m_1;
    } else if (N == 1) {
        Cheby = Cheby_n;
    } else {
        for (int i = 0; i < N-1; i++) {
            Cheby = 2.0*x*Cheby_n - Cheby_n_m_1;
            Cheby_n_m_1 = Cheby_n;
            Cheby_n = Cheby;
        }
    }
    
    return Cheby;
}

//****************************************************************************
// This routine returns roots of the Chebyshev polynomial of the first kind 
// also known as Chebyshev points of the first kind or Chebyshev nodes, or, 
// more formally, ChebyshevGauss points
//****************************************************************************
long double chebyshev_points_first_kind ( int index_node, int n) {
    // n represents the degree of the polynomial whose 
    // zeros are being computed
    long double angle;
    long double x;
    int k;
    long double pi_val = PI_CHEBY;
    
    k = (n - 1) - index_node; // k = 0 ... n - 1
    
    angle = ( long double) ( 2 * k + 1 ) * pi_val / ( long double ) ( 2 * n );
    x = cos ( angle );
    
    return x;
}

//****************************************************************************
// This routine returns extrema of the Chebyshev polynomial of the first kind 
// also known as Chebyshev extreme points, or Chebyshev-Lobatto points
//****************************************************************************
long double chebyshev_points_second_kind ( int index_node, int n ) {
    // n represents the degree of the polynomial whose 
    // extremas are being computed, which is a Lobatto polynomial
    // of the form Lo(n,x) = (1 - x^2) * U(n-1,x) = (1 - x^2) * dT(n,x)/dx, where T
    // is the Chebyshev polynomial of the first kind
    long double angle;
    int i;
    long double x;
    int k;
    long double pi_val = PI_CHEBY;
    
    k = (n - 1) - index_node; // k = 0 ... n - 1
    
    angle = ( long double) ( k ) * pi_val / ( long double ) ( n - 1 );
    x = cos ( angle );
    
    return x;
}

long double zeros_shifted ( int index_node, int n, long double a, long double b, int node_distribution) {
    long double x;
    if (node_distribution == CHEBYSHEV_FIRST_KIND_DISTRIBUTION) {
        x = chebyshev_points_first_kind ( index_node, n);
    } else if (node_distribution == CHEBYSHEV_SECOND_KIND_DISTRIBUTION) {
        x = chebyshev_points_second_kind ( index_node, n);
    } else {
        cout << "Nodal distribution not specified" << endl;
    }
    
    x = 0.5 * ( a + b ) + x * 0.5 * ( b - a );
    return x;
}

void chebyshev_quadrature ( long double &weight, long double &x, int index_node, int n, long double a, long double b, int node_distribution) {
    long double pi_val = PI_CHEBY;
    if (node_distribution == CHEBYSHEV_FIRST_KIND_DISTRIBUTION) {
        x = chebyshev_points_first_kind ( index_node, n);
        weight = sqrt(1.0 - x*x)*pi_val/(n - 1);
    } else if (node_distribution == CHEBYSHEV_SECOND_KIND_DISTRIBUTION) {
        x = chebyshev_points_second_kind ( index_node, n);
        weight = pi_val/(n - 1);
        if (index_node == 0 || index_node == n - 1) {
           weight *= 2.0; 
        }
        weight *= sqrt(1.0 - x*x);
    } else {
        cout << "Nodal distribution not specified" << endl;
    }
    
    x = 0.5 * ( a + b ) + x * 0.5 * ( b - a );
    weight *= 0.5 * ( b - a );
}

void chebyshev_quadrature ( long double *weight, long double *x, int n, long double a, long double b, int node_distribution) {
    long double weight_temp, x_temp;
    for (int i = 0; i < n; i++) {
        chebyshev_quadrature ( weight_temp, x_temp, i, n, a, b, node_distribution);
        weight[i] = weight_temp;
        x[i] = x_temp;
    }
    
}

long double Uniform_Distribution(const int &i, const int &Np, const long double &val_min, const long double &val_max) {
    long double val;
     
    val = val_min + i*(val_max - val_min)/Np;
    
    if (Np == 0) {
        val = val_min;
    }
    
    return val;
}

long double Uniform_Distribution_No_Endpoint(const int &i, const int &Np, const long double &val_min, const long double &val_max) {
    long double val;
     
    val = val_min + (i + 1)*(val_max - val_min)/(Np + 2);
    
    return val;
}
