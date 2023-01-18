/*******************************************************************
  File: M1_Model_1D_Utilities.cc

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

#include "M1_Model_1D_Utilities.h"

long double binomialCoefficients(const long double &n, const long double &k) {
    long double c;
    if (k < 0) {
        return 0;
    }
    
    if (k == 0 || k == n) {
        return 1;
    }
    
    c = 1.0;
    for (int i = 0; i < k; i++) {
        c = c * (n - i) / (i + 1);
    }
    return c;
}

int SH_Linear_Index(const int &Order_SH, const int &degree_SH) {
   int index_SH;
   index_SH = Order_SH*Order_SH + Order_SH + degree_SH;
   return index_SH;
}

int factorial(const int &n) { 
    if (n == 0) 
        return 1; 
    return n * factorial(n - 1); 
} 

int factorial(const int &n, const int &a) { 
    if (n == a) 
        return 1; 
    return n * factorial(n - 1, a); 
} 

long double factorial_inv(const int &n, const int &a) { 
    if (n == a) 
        return 1.0; 
    return (long double) (1.0/n) * factorial_inv(n - 1, a); 
} 


long double factorial_ratios(const int &num, const int &den) { 
    long double ratio;
//         cout << "num = " << num << "   " << "den = " << den << endl;
    if (num > den) { 
        ratio = (long double)(factorial(num, den));
    } else if (num < den) {
        ratio = (long double)(factorial_inv(den, num));
//         cout << "ratio = " << ratio << endl;
    } else {
        ratio = 1.0;
    }
    
    return ratio;
}

int double_factorial(const int &n) { 
    if (n == 0 || n == 1) 
        return 1; 
    return n * double_factorial(n - 2); 
} 

long double SH_Normalization_Constant(const int &l, const int &m) {
    long double K_l_m;
    K_l_m = (2.0*l+1.0)/(4.0*PI);
//     K_l_m *= (long double)(factorial(l - m))/(long double)(factorial(l + m));
    K_l_m *= factorial_ratios(l - m, l + m);
    K_l_m = sqrt(K_l_m);
    
//     cout << "factorial_ratios(l - m, l + m) = " << factorial_ratios(l - m, l + m) << endl;
    
    return K_l_m;
}

long double P_l_m_Polynomials(const long double &x, const int &l, const int &m) {
    long double P_l_m;
    switch (l) {
        case 0:
            switch (m) {
                case 0:
                    P_l_m = 1.0;
                    break;
            };
            break;
        default:
            if (l == m) {
                P_l_m = pow(-1.0, m)*double_factorial(2.0*m-1.0)*pow((1.0 - pow(x, 2)), m/2.0);
            } else if (l == m+1) {
                P_l_m = x*(2.0*m+1.0)*P_l_m_Polynomials(x, l-1, m);
            }
            break;
    };
    return P_l_m;
}

long double Associated_Legendre_Polynomials(const long double &x, const int &l, const int &m) {
    long double P_l_m;
    long double P_lminus1_m, P_lminus2_m;
    long double coeff1, coeff2;
    
    if (l == fabs(m) || l == fabs(m)+1) {
        P_l_m = P_l_m_Polynomials(x, l, fabs(m));
    } else {
        P_lminus2_m = Associated_Legendre_Polynomials(x, l - 2, fabs(m));
        P_lminus1_m = Associated_Legendre_Polynomials(x, l - 1, fabs(m));
        
        coeff1 = x*(long double)(2.0*l - 1.0)/(long double)(l - fabs(m));
        coeff2 = (long double)(l + fabs(m) - 1.0)/(long double)(l - fabs(m));
            
        P_l_m = coeff1*P_lminus1_m - coeff2*P_lminus2_m;
    }
    
    return P_l_m;
}

long double Precompute_SH_Basis_Coeffs(const long double &l, const long double &m, const long double &i, const long double &j, const long double &k) {
    long double Binomial_Coeffs;
    long double factorial_ratio;
    long double SH_Basis_Coeffs;
    
    if (m == 0) {
        SH_Basis_Coeffs = pow(2,l);
        SH_Basis_Coeffs *= sqrt((2.0*l+1.0)/(4.0*PI));
        
        Binomial_Coeffs = binomialCoefficients(l, k);
        SH_Basis_Coeffs *= Binomial_Coeffs;
        
        Binomial_Coeffs = binomialCoefficients((long double)((l+k-1)/2.0), l);
        SH_Basis_Coeffs *= Binomial_Coeffs;
    } else if (m > 0) {
        SH_Basis_Coeffs = m*pow(-1.0,k+j)*pow(2.0, l+m-2*j-1);
        SH_Basis_Coeffs *= sqrt(2.0)*sqrt((2.0*l+1.0)/(4.0*PI));
        
        factorial_ratio = factorial_ratios(l - m, l + m);
        SH_Basis_Coeffs *= sqrt(factorial_ratio);
        
        factorial_ratio = factorial_ratios(i, i - m);
        SH_Basis_Coeffs *= factorial_ratio;
        
        factorial_ratio = factorial_ratios(m-j-1, m-2*j);
        SH_Basis_Coeffs *= factorial_ratio/factorial(j);
        
        Binomial_Coeffs = binomialCoefficients(l, i);
        SH_Basis_Coeffs *= Binomial_Coeffs;
        
        Binomial_Coeffs = binomialCoefficients((long double)((l+i-1)/2.0), l);
        SH_Basis_Coeffs *= Binomial_Coeffs;
        
        Binomial_Coeffs = binomialCoefficients(j, k);
        SH_Basis_Coeffs *= Binomial_Coeffs;
    } else {
        SH_Basis_Coeffs = pow(-1.0,k+j)*pow(2.0, l+fabs(m)-2*j-1);
        SH_Basis_Coeffs *= sqrt(2.0)*sqrt((2.0*l+1.0)/(4.0*PI));
        
        factorial_ratio = factorial_ratios(l - fabs(m), l + fabs(m));
        SH_Basis_Coeffs *= sqrt(factorial_ratio);
        
        factorial_ratio = factorial_ratios(i, i - fabs(m));
        SH_Basis_Coeffs *= factorial_ratio;
        
        Binomial_Coeffs = binomialCoefficients((long double)(l), i);
        SH_Basis_Coeffs *= Binomial_Coeffs;
        
        Binomial_Coeffs = binomialCoefficients((long double)((l+i-1)/2.0), l);
        SH_Basis_Coeffs *= Binomial_Coeffs;
        
        Binomial_Coeffs = binomialCoefficients((long double)(fabs(m)-j-1), j);
        SH_Basis_Coeffs *= Binomial_Coeffs;
        
        Binomial_Coeffs = binomialCoefficients((long double)(j), k);
        SH_Basis_Coeffs *= Binomial_Coeffs;
    }
    
    return SH_Basis_Coeffs;
}

long double Precompute_First_Kind_Chebyshev_Basis_Coeffs(const int &n, const int &k) {
    long double factorial_ratio;
    long double Chebyshev_Basis_Coeffs;
    
    if (n == 0) {
        Chebyshev_Basis_Coeffs = 1.0;
    } else {
        Chebyshev_Basis_Coeffs = (n/2.0)*pow(-1.0,k)*pow(2.0, n-2*k);
        
        factorial_ratio = factorial_ratios(n-k-1, n-2*k);
        Chebyshev_Basis_Coeffs *= factorial_ratio/(long double)(factorial(k));
    }
    
    return Chebyshev_Basis_Coeffs;
}

void Chebyshev_First_Kind_to_Monomial_Basis_ratio_E_No_SH(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_mu) {
    long double Coeffs_Cheby;
    long double *Coefficients_Fit_Monomials_Basis;
    int N_Coeffs_Total = N_Coeffs_E*N_Coeffs_f*N_Coeffs_mu;
    int index_Monomial_Basis, index_Orthog_Basis;
    
    Coefficients_Fit_Monomials_Basis = new long double[N_Coeffs_Total];
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Monomials_Basis[i] = 0.0;
    }
    
    for (int i_fit_f = 0; i_fit_f < N_Coeffs_f; i_fit_f++) {
        for (int i_fit_mu = 0; i_fit_mu < N_Coeffs_mu; i_fit_mu++) {
            
            for (int i_fit_E = 0; i_fit_E < N_Coeffs_E; i_fit_E++) {
                index_Monomial_Basis = (i_fit_E*N_Coeffs_f + i_fit_f)*N_Coeffs_mu + i_fit_mu;
                
                for (int i_fit_E_Orthog = 0; i_fit_E_Orthog < N_Coeffs_E; i_fit_E_Orthog++) {
                    index_Orthog_Basis = (i_fit_E_Orthog*N_Coeffs_f + i_fit_f)*N_Coeffs_mu + i_fit_mu;
                    
                    for (int k = 0; k <= floor(i_fit_E_Orthog/2.0); k++) {
                        if (i_fit_E_Orthog - 2*k == i_fit_E) {
                            Coeffs_Cheby = Precompute_First_Kind_Chebyshev_Basis_Coeffs(i_fit_E_Orthog, k);
                            Coefficients_Fit_Monomials_Basis[index_Monomial_Basis] += Coeffs_Cheby*Coefficients_Fit_Orthog_Basis[index_Orthog_Basis];
                        }
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Orthog_Basis[i] = Coefficients_Fit_Monomials_Basis[i];
    }
    
    delete[] Coefficients_Fit_Monomials_Basis;
}

void Chebyshev_First_Kind_to_Monomial_Basis_Norm_f_No_SH(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_mu) {
    long double Coeffs_Cheby;
    long double *Coefficients_Fit_Monomials_Basis;
    int N_Coeffs_Total = N_Coeffs_E*N_Coeffs_f*N_Coeffs_mu;
    int index_Monomial_Basis, index_Orthog_Basis;
    
    Coefficients_Fit_Monomials_Basis = new long double[N_Coeffs_Total];
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Monomials_Basis[i] = 0.0;
    }
    
    for (int i_fit_E = 0; i_fit_E < N_Coeffs_E; i_fit_E++) {
        for (int i_fit_mu = 0; i_fit_mu < N_Coeffs_mu; i_fit_mu++) {
            
            for (int i_fit_f = 0; i_fit_f < N_Coeffs_f; i_fit_f++) {
                index_Monomial_Basis = (i_fit_E*N_Coeffs_f + i_fit_f)*N_Coeffs_mu + i_fit_mu;
                
                for (int i_fit_f_Orthog = 0; i_fit_f_Orthog < N_Coeffs_f; i_fit_f_Orthog++) {
                    index_Orthog_Basis = (i_fit_E*N_Coeffs_f + i_fit_f_Orthog)*N_Coeffs_mu + i_fit_mu;
                    for (int k = 0; k <= floor(i_fit_f_Orthog/2.0); k++) {
                        if (i_fit_f_Orthog - 2*k == i_fit_f) {
                            Coeffs_Cheby = Precompute_First_Kind_Chebyshev_Basis_Coeffs(i_fit_f_Orthog, k);
                            Coefficients_Fit_Monomials_Basis[index_Monomial_Basis] += Coeffs_Cheby*Coefficients_Fit_Orthog_Basis[index_Orthog_Basis];
                        }
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Orthog_Basis[i] = Coefficients_Fit_Monomials_Basis[i];
    }
    
    delete[] Coefficients_Fit_Monomials_Basis;
}

void Chebyshev_First_Kind_to_Monomial_Basis_mu_No_SH(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_mu) {
    long double Coeffs_Cheby;
    long double *Coefficients_Fit_Monomials_Basis;
    int N_Coeffs_Total = N_Coeffs_E*N_Coeffs_f*N_Coeffs_mu;
    int index_Monomial_Basis, index_Orthog_Basis;
    
    Coefficients_Fit_Monomials_Basis = new long double[N_Coeffs_Total];
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Monomials_Basis[i] = 0.0;
    }
    
    for (int i_fit_E = 0; i_fit_E < N_Coeffs_E; i_fit_E++) {
        for (int i_fit_f = 0; i_fit_f < N_Coeffs_f; i_fit_f++) {
        
            for (int i_fit_mu = 0; i_fit_mu < N_Coeffs_mu; i_fit_mu++) {
                index_Monomial_Basis = (i_fit_E*N_Coeffs_f + i_fit_f)*N_Coeffs_mu + i_fit_mu;
                
                for (int i_fit_mu_Orthog = 0; i_fit_mu_Orthog < N_Coeffs_mu; i_fit_mu_Orthog++) {
                    index_Orthog_Basis = (i_fit_E*N_Coeffs_f + i_fit_f)*N_Coeffs_mu + i_fit_mu_Orthog;
                    for (int k = 0; k <= i_fit_mu_Orthog; k++) {
                        if (2*i_fit_mu_Orthog - 2*k == 2*i_fit_mu) {
                            Coeffs_Cheby = Precompute_First_Kind_Chebyshev_Basis_Coeffs(2*i_fit_mu_Orthog, k);
                            Coefficients_Fit_Monomials_Basis[index_Monomial_Basis] += Coeffs_Cheby*Coefficients_Fit_Orthog_Basis[index_Orthog_Basis];
                        }
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Orthog_Basis[i] = Coefficients_Fit_Monomials_Basis[i];
    }
    
    delete[] Coefficients_Fit_Monomials_Basis;
}

void Chebyshev_First_Kind_to_Monomial_Basis_ratio_E(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_SH) {
    long double Coeffs_Cheby;
    long double *Coefficients_Fit_Monomials_Basis;
    int N_Coeffs_Total = N_Coeffs_E*N_Coeffs_f*N_Coeffs_SH;
    int N_Coeffs_Total_Without_SH = N_Coeffs_E*N_Coeffs_f;
    int index_Monomial_Basis, index_Orthog_Basis;
    
    Coefficients_Fit_Monomials_Basis = new long double[N_Coeffs_Total];
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Monomials_Basis[i] = 0.0;
    }
    
    for (int i_fit_f = 0; i_fit_f < N_Coeffs_f; i_fit_f++) {
        for (int i_fit_SH = 0; i_fit_SH < N_Coeffs_SH; i_fit_SH++) {
            
            for (int i_fit_E = 0; i_fit_E < N_Coeffs_E; i_fit_E++) {
                index_Monomial_Basis = i_fit_E*N_Coeffs_f + i_fit_f;
                
                for (int i_fit_E_Orthog = 0; i_fit_E_Orthog < N_Coeffs_E; i_fit_E_Orthog++) {
                    index_Orthog_Basis = i_fit_E_Orthog*N_Coeffs_f + i_fit_f;
                    
                    for (int k = 0; k <= floor(i_fit_E_Orthog/2.0); k++) {
                        if (i_fit_E_Orthog - 2*k == i_fit_E) {
                            Coeffs_Cheby = Precompute_First_Kind_Chebyshev_Basis_Coeffs(i_fit_E_Orthog, k);
                            Coefficients_Fit_Monomials_Basis[i_fit_SH*N_Coeffs_Total_Without_SH+index_Monomial_Basis] += Coeffs_Cheby*Coefficients_Fit_Orthog_Basis[i_fit_SH*N_Coeffs_Total_Without_SH+index_Orthog_Basis];
                        }
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Orthog_Basis[i] = Coefficients_Fit_Monomials_Basis[i];
    }
    
    delete[] Coefficients_Fit_Monomials_Basis;
}

void Chebyshev_First_Kind_to_Monomial_Basis_Norm_f(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_SH) {
    long double Coeffs_Cheby;
    long double *Coefficients_Fit_Monomials_Basis;
    int N_Coeffs_Total = N_Coeffs_E*N_Coeffs_f*N_Coeffs_SH;
    int N_Coeffs_Total_Without_SH = N_Coeffs_E*N_Coeffs_f;
    int index_Monomial_Basis, index_Orthog_Basis;
    
    Coefficients_Fit_Monomials_Basis = new long double[N_Coeffs_Total];
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Monomials_Basis[i] = 0.0;
    }
    
    for (int i_fit_E = 0; i_fit_E < N_Coeffs_E; i_fit_E++) {
        for (int i_fit_SH = 0; i_fit_SH < N_Coeffs_SH; i_fit_SH++) {
            
            for (int i_fit_f = 0; i_fit_f < N_Coeffs_f; i_fit_f++) {
                index_Monomial_Basis = i_fit_E*N_Coeffs_f + i_fit_f;
                
                for (int i_fit_f_Orthog = 0; i_fit_f_Orthog < N_Coeffs_f; i_fit_f_Orthog++) {
                    index_Orthog_Basis = i_fit_E*N_Coeffs_f + i_fit_f_Orthog;
                    for (int k = 0; k <= i_fit_f_Orthog; k++) {
                        if (2*i_fit_f_Orthog - 2*k == 2*i_fit_f) {
                            Coeffs_Cheby = Precompute_First_Kind_Chebyshev_Basis_Coeffs(2*i_fit_f_Orthog, k);
                            Coefficients_Fit_Monomials_Basis[i_fit_SH*N_Coeffs_Total_Without_SH+index_Monomial_Basis] += Coeffs_Cheby*Coefficients_Fit_Orthog_Basis[i_fit_SH*N_Coeffs_Total_Without_SH+index_Orthog_Basis];
                        }
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Orthog_Basis[i] = Coefficients_Fit_Monomials_Basis[i];
    }
    
    delete[] Coefficients_Fit_Monomials_Basis;
}

void Spherical_Harmonics_to_Monomial_Basis(long double *Coefficients_Fit_Orthog_Basis, const int &N_Coeffs_E, const int &N_Coeffs_f, const int &N_Coeffs_SH, const int &Order_SH, const int *Array_l_SH, const int *Array_m_SH) {
    long double Coeffs_SH;
    long double *Coefficients_Fit_Monomials_Basis;
    int N_Coeffs_Total = N_Coeffs_E*N_Coeffs_f*N_Coeffs_SH;
    int N_Coeffs_Total_Without_SH = N_Coeffs_E*N_Coeffs_f;
    int index_Basis;
    
    int index_SH_Monomial_Basis, index_SH_Orthog_Basis;
    
    Coefficients_Fit_Monomials_Basis = new long double[N_Coeffs_Total];
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Monomials_Basis[i] = 0.0;
    }
    
    for (int i_fit_E = 0; i_fit_E < N_Coeffs_E; i_fit_E++) {
        for (int i_fit_f = 0; i_fit_f < N_Coeffs_f; i_fit_f++) {
            index_Basis = i_fit_E*N_Coeffs_f + i_fit_f;
            
            index_SH_Monomial_Basis = 0;
            for (int p = 0; p <= Order_SH; p+=2) {
                for (int q = 0; q <= Order_SH - p; q+=2) {
                    index_SH_Orthog_Basis = 0;
                    for (int l = 0; l <= Order_SH; l+=2) {
                        for (int m = 0; m <= l; m+=2) {
                            if (m == 0) {
                                for (int k = 0; k <= l; k++) {
                                    if (p == k && q == 0) {
                                        Coeffs_SH = Precompute_SH_Basis_Coeffs(l, m, 0, 0, k);
                                        Coefficients_Fit_Monomials_Basis[index_SH_Monomial_Basis*N_Coeffs_Total_Without_SH+index_Basis] += Coeffs_SH*Coefficients_Fit_Orthog_Basis[index_SH_Orthog_Basis*N_Coeffs_Total_Without_SH+index_Basis];
                                    }
                                }
                            } else {
                                for (int i = m; i <= l; i++) {
                                    for (int j = 0; j <= floor(m/2.0); j++) {
                                        for (int k = 0; k <= j; k++) {
                                            if (p == i-m+2*k && q == m - 2*j) {
                                                Coeffs_SH = Precompute_SH_Basis_Coeffs(l, m, i, j, k);
                                                Coefficients_Fit_Monomials_Basis[index_SH_Monomial_Basis*N_Coeffs_Total_Without_SH+index_Basis] += Coeffs_SH*Coefficients_Fit_Orthog_Basis[index_SH_Orthog_Basis*N_Coeffs_Total_Without_SH+index_Basis];  
                                            } 
//                                             else if ((i-m+2*k) % 2 != 0 || (m-2*j) % 2 != 0) {
//                                                 Coeffs_SH = Precompute_SH_Basis_Coeffs(l, m, i, j, k);
//                                                 cout << "l = " << l << "   " << "m = " << m << "   " << "(i-m+2*k) = " << (i-m+2*k) << "   " << "(m-2*j) = " << (m-2*j) << "   " << "Coeffs_SH = " << Coeffs_SH << endl;
//                                             }
                                        }
                                    }
                                }
                            }
                            index_SH_Orthog_Basis++;
                        }
                    }   
                    index_SH_Monomial_Basis++;
                }
            }
        }
    }
    
    for (int i = 0; i < N_Coeffs_Total; i++) {
        Coefficients_Fit_Orthog_Basis[i] = Coefficients_Fit_Monomials_Basis[i];
    }
    
    delete[] Coefficients_Fit_Monomials_Basis;
}

long double Cartesian_Spherical_harmonics(const long double &x, const long double &y, const long double &z, const int &l, const int &m) {
    long double Coeffs_SH, Cartes_SH;
    long double Cartes_SH_j, Cartes_SH_k;
    long double x2, z2;
    
    if (m == 0) {
        for (int k = l; k >= 0; k--) {
            Coeffs_SH = Precompute_SH_Basis_Coeffs(l, m, 0, 0, k);
//             cout << "k = " << k << "   " << "l = " << l << "   " << "Coeffs_SH = " << Coeffs_SH << endl;
            if (k == l) {
                Cartes_SH = Coeffs_SH;
            } else {
                Cartes_SH = Coeffs_SH + z*Cartes_SH;
            }
        }
    } else if (m > 0) {
        x2 = x*x;
        z2 = z*z;
        for (int i = l; i >= m; i--) {
            for (int j = 0; j <= floor(m/2); j++) {
                for (int k = j; k >= 0; k--) {
                    Coeffs_SH = Precompute_SH_Basis_Coeffs(l, m, i, j, k);
                    if (k == j) {
                        Cartes_SH_k = Coeffs_SH;
                    } else {
                        Cartes_SH_k = Coeffs_SH + z2*Cartes_SH_k;
                    }
                }
                if (j == 0) {
                    Cartes_SH_j = Cartes_SH_k;
                } else {
                    Cartes_SH_j = Cartes_SH_k + x2*Cartes_SH_j;
                }
            }
            
            if ((m % 2) != 0) { //then m is odd
                Cartes_SH_j *= x;
            }
            
            if (i == l) {
                Cartes_SH = Cartes_SH_j;
            } else {
                Cartes_SH = Cartes_SH_j + z*Cartes_SH;
            }
        }
    } else {
        x2 = x*x;
        z2 = z*z;
        for (int i = l; i >= fabs(m); i--) {
            for (int j = 0; j <= floor((fabs(m)-1)/2); j++) {
                // j should take values up to (m-1)/2
                for (int k = j; k >= 0; k--) {
                    Coeffs_SH = Precompute_SH_Basis_Coeffs(l, m, i, j, k);
                    if (k == j) {
                        Cartes_SH_k = Coeffs_SH;
                    } else {
                        Cartes_SH_k = Coeffs_SH + z2*Cartes_SH_k;
                    }
                }
                
                if (j == 0) {
                    Cartes_SH_j = Cartes_SH_k;
                } else {
                    Cartes_SH_j = Cartes_SH_k + x2*Cartes_SH_j;
                }
            }
            if ((-m % 2) == 0) { //then m is even
                Cartes_SH_j *= x;
            }
            if (i == l) {
                Cartes_SH = Cartes_SH_j;
            } else {
                Cartes_SH = Cartes_SH_j + z*Cartes_SH;
            }
        }
        Cartes_SH *= y;
    }
    
//     cout << "l = " << l << "   " << "m = " << m << "   " << "x = " << x << "   " << "y = " << y << "   " << "z = " << z << "   " << "Cartes_SH = " << Cartes_SH << endl;
    
    return Cartes_SH;
}

// long double Cartesian_Spherical_harmonics(const long double &x, const long double &y, const long double &z, const int &l, const int &m) {
//     long double K_l_m, Y_l_m;
//     long double cos_m_phi, sin_minus_m_phi;
//     long double cos_phi, sin_phi, cos_theta;
//     
// //     cos_phi = x/sqrt(1.0 - pow(z,2));
// //     sin_phi = y/sqrt(1.0 - pow(z,2));
// //     cos_theta = z;
// //     
// //     if (fabs(z) > 1.0 - 1.0e-4) {
// //        cos_phi = 1.0;
// //        sin_phi = 0.0;
// //     }
//     
//     cos_phi = y/sqrt(1.0 - pow(x,2));
//     sin_phi = z/sqrt(1.0 - pow(x,2));
//     cos_theta = x;
//     
//     if (fabs(x) > 1.0 - 1.0e-4) {
//        cos_phi = 1.0;
//        sin_phi = 0.0;
//     }
//     
//     K_l_m = SH_Normalization_Constant(l,fabs(m));
//     
//     if (m > 0) {
//         cos_m_phi = Chebyshev_Polynomial_Basis(cos_phi, m);
//         Y_l_m = pow(-1.0, m)*sqrt(2.0)*K_l_m*cos_m_phi*Associated_Legendre_Polynomials(cos_theta,l, m);
//     } else if (m == 0) {
//         Y_l_m = K_l_m*Associated_Legendre_Polynomials(cos_theta,l, m);
//     } else {
//         sin_minus_m_phi = sin_phi*Chebyshev_Second_Kind_Polynomial_Basis(cos_phi, fabs(m)-1);
//         Y_l_m = pow(-1.0, m)*sqrt(2.0)*K_l_m*sin_minus_m_phi*Associated_Legendre_Polynomials(cos_theta,l, fabs(m));
//         
// //         cout << "pow(-1.0, m) = " << pow(-1.0, m) << "   " << "l = " << l << "    " << "m = " << m << "    " << "sin_minus_m_phi = " << sin_minus_m_phi << "    " << "sin_phi = " << sin_phi << "    " << "Chebyshev = " << Chebyshev_Second_Kind_Polynomial_Basis(cos_phi, -m-1) << "    " << "cos_phi = " << cos_phi << "    " << "K_l_m = " << K_l_m << "   " << "Assoc_Leg_Pol = " << Associated_Legendre_Polynomials(cos_theta,l, fabs(m)) << endl;
//     }
//     return Y_l_m;
// }

long double VanderMonde_Matrix_1_var(const long double &var_1, const int &Order_poly_1) {
    long double Van_Entry;
    long double poly_1;
    poly_1 = Chebyshev_Polynomial_Basis(var_1, Order_poly_1);
    
    Van_Entry = poly_1;
    return Van_Entry;
}

long double VanderMonde_Matrix_2_vars(const long double &var_1, const long double &var_2, const int &Order_poly_1, const int &Order_poly_2) {
    long double Van_Entry;
    long double poly_1, poly_2;
    poly_1 = Chebyshev_Polynomial_Basis(var_1, Order_poly_1);
    poly_2 = Chebyshev_Polynomial_Basis(var_2, Order_poly_2);
    
    Van_Entry = poly_1*poly_2;
    return Van_Entry;
}

long double VanderMonde_Vector_N_vars(const int &Index_Entry, const int &Index_Point) {
    long double Van_Vector;
    
    if (Index_Entry == Index_Point) {
        Van_Vector = 1.0;
    } else {
        Van_Vector = 0.0;
    }
    return Van_Vector;
}

void Solve_A_x_b(const long double *VanderMonde_Matrix, long double *Coeff_Vand_Matrix, const long double *VanderMonde_Vector, const int &n, const int &index_Cheby) {
    // First perform Cholesky decomposition of A
    long double *Temp_Mat_L, *Temp_Mat_U, *Temp_Vec;
    Temp_Mat_L = new long double[n*n];
    Temp_Mat_U = new long double[n*n];
    Temp_Vec = new long double[n];
    
    LUdecomposition(VanderMonde_Matrix, Temp_Mat_L, Temp_Mat_U, n);
    
    LU_Solve_Forward(n, Temp_Vec, Temp_Mat_L, VanderMonde_Vector);
    
    LU_Solve_Back(n, Temp_Vec, Temp_Mat_U, Temp_Vec);
    
    for (int i = 0; i < n; i++) {
        Coeff_Vand_Matrix[i*n + index_Cheby] = Temp_Vec[i];
    }
    
    delete[] Temp_Mat_L;
    delete[] Temp_Mat_U;
    delete[] Temp_Vec;
}

void Check_A_x_b(const long double *VanderMonde_Matrix, const long double *Coeff_Vand_Matrix, const int &n) {
    // First perform Cholesky decomposition of A
    long double temp_val;
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            temp_val = 0.0;
            for (int k = 0; k < n; k++) {
                temp_val += VanderMonde_Matrix[i*n + k] * Coeff_Vand_Matrix[k*n + j];
            }
            if (i == j) {
                if (fabs(1.0 - temp_val) > 1.0e-6) {
                    cout << "Issue with Vandermonde System:" << endl;
                    cout << "i = " << i << "  " << "j = " << j << "  " << "temp_val = " << temp_val << endl;
                    exit(0);
                }
            } else {
                if (fabs(temp_val) > 1.0e-6) {
                    cout << "Issue with Vandermonde System:" << endl;
                    cout << "i = " << i << "  " << "j = " << j << "  " << "temp_val = " << temp_val << endl;
                    exit(0);
                }
            }
        }
    }
}

void LUdecomposition(const long double *A, long double *L, long double *U, const int &n) {
   for (int i = 0; i < n; i++) {
       for (int j = 0; j < n; j++) {
           if (j < i) {
               L[j*n+i] = 0;
           } else {
               L[j*n+i] = A[j*n+i];
               for (int k = 0; k < i; k++) {
                   L[j*n+i] = L[j*n+i] - L[j*n+k] * U[k*n+i];
            }
         }
      }
      for (int j = 0; j < n; j++) {
          if (j < i) {
              U[i*n+j] = 0;
          } else if (j == i) {
              U[i*n+j] = 1;
          } else {
              U[i*n+j] = A[i*n+j] / L[i*n+i];
              for (int k = 0; k < i; k++) {
                  U[i*n+j] = U[i*n+j] - ((L[i*n+k] * U[k*n+j]) / L[i*n+i]);
            }
         }
      }
   }
}

void LU_Solve_Back(const int &n, long double *x, const long double *A, const long double *b)
{
    // Back solve A x = b
    for (int i = n-1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i+1; j < n; j++) {
            x[i] -= A[i*n+j] * x[j];
        }
//         if (fabs(A[i*n+i]) < 1.0e-10) {
//             x[i] = 0.0;
//         } else {
            x[i] /= A[i*n+i];
//         }
    }
}

void LU_Solve_Forward(const int &n, long double *x, const long double *A, const long double *b)
{
    // Forward solve A x = b
    for (int i = 0; i < n; i++) {
        x[i] = b[i];
        for (int j = 0; j < i; j++) {
            x[i] -= A[i*n+j] * x[j];
        }
        
//         if (fabs(A[i*n+i]) < 1.0e-10) {
//             x[i] = 0.0;
//         } else {
            x[i] /= A[i*n+i];
//         }
    }
}

void Finite_Difference_Fit(long double *vec_diff, const long double *mat_u_y, const int &finite_difference_type, const int &index_e, const int &index_f, const int &N_points_E, const int &N_points_f, const int &N_points_Theta, const int &N_points_Phi, const long double &h, const int &order, const int &VAR_NUM) {
    int N_pts;
    
    switch (VAR_NUM) {
        case VAR_E :
            N_pts = N_points_E;
            break;
        case VAR_N1 :
            N_pts = N_points_f;
            break;
        case VAR_THETA :
            N_pts = N_points_Theta;
            break;
        case VAR_PHI :
            N_pts = N_points_Phi;
            break;
    }
    
//     switch (finite_difference_type) {
//         case FORWARD:
//             for (int i = 0; i < N_pts; i++) {
//                 vec_diff[i] = forward_finite_difference(mat_u_y, index_e, index_f, N_points_E, N_points_f, N_points_Theta, N_points_Phi, h, order, VAR_NUM);
//             }
//             break;
//         case BACKWARD:
//             for (int i = 0; i < N_pts; i++) {
//                 vec_diff[i] = backward_finite_difference(mat_u_y, index_e, index_f, N_points_E, N_points_f, N_points_Theta, N_points_Phi, h, order, VAR_NUM);
//             }
//             break;
//         case CENTRAL:
//             cout << endl;
//             break;
//     };
}

long double forward_finite_difference(const long double *u_y, const int &index_e, const int &index_f, const int &N_points_E, const int &N_points_f, const int &N_points_Theta, const int &N_points_Phi, const long double &h, const int &order, const int &VAR_NUM) {
    long double diff_1 = 0.0;
    int prec;
    int n;
    long double *c, *x;
    prec = 3;
    n = order + prec;
    c = new long double[n];
    x = new long double[n];
    
    differ_forward ( h, order, prec, c, x );
    // forward difference for points on left boundary
    switch (VAR_NUM) {
        case VAR_N1 :
            for (int i = 0; i < n; i++) {
                diff_1 += c[i] * u_y[index_e * N_points_f + i];
            }
            break;
        case VAR_E :
            for (int i = 0; i < n; i++) {
                diff_1 += c[i] * u_y[i * N_points_f + index_f];
            }
            break;
    }
    delete[] c;
    delete[] x;
    return diff_1;
}

long double backward_finite_difference(const long double *u_y, const int &index_e, const int &index_f, const int &N_points_E, const int &N_points_f, const int &N_points_Theta, const int &N_points_Phi, const long double &h, const int &order, const int &VAR_NUM) {
    int prec;
    int n;
    long double diff_1 = 0.0;
    long double *c, *x;
    prec = 5;
    n = order + prec;
    c = new long double[n];
    x = new long double[n];
    
    differ_backward ( h, order, prec, c, x );
    // forward difference for points on left boundary
    switch (VAR_NUM) {
        case VAR_N1 :
            for (int i = 0; i < n; i++) {
                diff_1 += c[n - i - 1] * u_y[index_e * N_points_f + (N_points_f - 1) - i];
            }
            break;
        case VAR_E :
            for (int i = 0; i < n; i++) {
                diff_1 += c[n - i - 1] * u_y[((N_points_E - 1) - i) * N_points_f + index_f];
            }
            break;
    }
    delete[] c;
    delete[] x;
    
    return diff_1;
}

long double central_finite_difference(const long double *u_y, const int &index_e, const int &index_f, const int &N_points_E, const int &N_points_f, const int &N_points_Theta, const int &N_points_Phi, const long double &h, const int &order) {
    int prec;
    int n;
    long double diff_1 = 0.0;
    long double *c, *x;
    prec = 3;
    n = order + prec;
    c = new long double[n];
    x = new long double[n];
    
    differ_central ( h, order, prec, c, x );
    // forward difference for points on left boundary
    for (int i = 0; i < n; i++) {
        diff_1 += c[n - i - 1] * u_y[index_e * N_points_f + (N_points_f - 1) - i];
    }
    delete[] c;
    delete[] x;
    
    return diff_1;
}

int grid_E(const long double &r_E, const long double *rE_uniform, const int &Npts_E, const int &Npts_f, const int &index_N1) {
    int index_rE;
    int i_rE_L, i_rE_R;
    int iP_rE;
    i_rE_L = 0;
    i_rE_R = Npts_E - 1;
    iP_rE = (i_rE_R + i_rE_L)/2;
    
    while (i_rE_R - i_rE_L > 1) {
        if (r_E < rE_uniform[iP_rE*Npts_f + index_N1]){
            i_rE_R = iP_rE;
        } else {
            i_rE_L = iP_rE;
        }
        iP_rE = (i_rE_R + i_rE_L)/2;
    }
    index_rE = iP_rE;
    return index_rE;
}

// long double Cartesian_Spherical_harmonics(const long double &x, const long double &y, const long double &z, const int &Index) {
//     long double f_SH;
//     long double y_0_0;
//     long double y_1_m1, y_1_0, y_1_1;
//     long double y_2_m2, y_2_m1, y_2_0, y_2_1, y_2_2;
//     long double y_3_m3, y_3_m2, y_3_m1, y_3_0, y_3_1, y_3_2, y_3_3;
//     long double y_4_m4, y_4_m3, y_4_m2, y_4_m1, y_4_0, y_4_1, y_4_2, y_4_3, y_4_4;
//     
//     y_0_0 = sqrt(1.0/(4.0*PI));
//     y_1_m1 = sqrt(3.0/(4.0*PI))*y;
//     y_1_0 = sqrt(3.0/(4.0*PI))*z;
//     y_1_1 = sqrt(3.0/(4.0*PI))*x;
//     
//     y_2_m2 = sqrt(15.0/(4.0*PI))*x*y;
//     y_2_m1 = sqrt(15.0/(4.0*PI))*y*z;
//     y_2_0 = sqrt(5.0/(16.0*PI))*(3.0*pow(z,2) - 1.0);
//     y_2_1 = sqrt(15.0/(4.0*PI))*x*z;
//     y_2_2 = sqrt(15.0/(16.0*PI))*(pow(x,2) - pow(y,2));
//     
//     y_3_m3 = sqrt(35.0/(32.0*PI))*y*(3.0*pow(x,2) - pow(y,2));
//     y_3_m2 = sqrt(105.0/(4.0*PI))*x*y*z;
//     y_3_m1 = sqrt(21.0/(32.0*PI))*y*(5.0*pow(z,2) - 1.0);
//     y_3_0 = sqrt(7.0/(16.0*PI))*z*(5.0*pow(z,2) - 3.0);
//     y_3_1 = sqrt(21.0/(32.0*PI))*x*(5.0*pow(z,2) - 1.0);
//     y_3_2 = sqrt(105.0/(16.0*PI))*z*(pow(x,2) - pow(y,2));
//     y_3_3 = sqrt(35.0/(32.0*PI))*x*(pow(x,2) - 3.0*pow(y,2));
//     
//     y_4_m4 = 3.0*sqrt(35.0/(16.0*PI))*x*y*(pow(x,2) - pow(y,2));
//     y_4_m3 = 3.0*sqrt(35.0/(32.0*PI))*y*z*(3.0*pow(x,2) - pow(y,2));
//     y_4_m2 = 3.0*sqrt(5.0/(16.0*PI))*x*y*(7.0*pow(z,2) - 1.0);
//     y_4_m1 = 3.0*sqrt(5.0/(32.0*PI))*y*z*(7.0*pow(z,2) - 3.0);
//     y_4_0 = (3.0/16.0)*sqrt(1.0/PI)*(35.0*pow(z,4) - 30.0*pow(z,2) + 3.0);
//     y_4_1 = 3.0*sqrt(5.0/(32.0*PI))*x*z*(7.0*pow(z,2) - 3.0);
//     y_4_2 = (3.0/8.0)*sqrt(5.0/PI)*(pow(x,2) - pow(y,2))*(7.0*pow(z,2) - 1.0);
//     y_4_3 = (3.0/4.0)*sqrt(35.0/(2.0*PI))*x*z*(pow(x,2) - 3.0*pow(y,2));
//     y_4_4 = (3.0/16.0)*sqrt(35.0/PI)*(pow(x,2)*(pow(x,2) - 3.0*pow(y,2)) - pow(y,2)*(3.0*pow(x,2) - pow(y,2)));
//     
//     if (Index == 0) { 
//         f_SH = y_0_0;
//     } else if (Index == 1) {
//         f_SH = y_1_m1;
//     } else if (Index == 2) {
//         f_SH = y_1_0;
//     } else if (Index == 3) {
//         f_SH = y_1_1;
//     } else if (Index == 4) {
//         f_SH = y_2_m2;
//     } else if (Index == 5) {
//         f_SH = y_2_m1;
//     } else if (Index == 6) {
//         f_SH = y_2_0;
//     } else if (Index == 7) {
//         f_SH = y_2_1;
//     } else if (Index == 8) {
//         f_SH = y_2_2;
//     } else if (Index == 9) {
//         f_SH = y_3_m3;
//     } else if (Index == 10) {
//         f_SH = y_3_m2;
//     } else if (Index == 11) {
//         f_SH = y_3_m1;
//     } else if (Index == 12) {
//         f_SH = y_3_0;
//     } else if (Index == 13) {
//         f_SH = y_3_1;
//     } else if (Index == 14) {
//         f_SH = y_3_2;
//     } else if (Index == 15) {
//         f_SH = y_3_3;
//     } else if (Index == 16) {
//         f_SH = y_4_m4;
//     } else if (Index == 17) {
//         f_SH = y_4_m3;
//     } else if (Index == 18) {
//         f_SH = y_4_m2;
//     } else if (Index == 19) {
//         f_SH = y_4_m1;
//     } else if (Index == 20) {
//         f_SH = y_4_0;
//     } else if (Index == 21) {
//         f_SH = y_4_1;
//     } else if (Index == 22) {
//         f_SH = y_4_2;
//     } else if (Index == 23) {
//         f_SH = y_4_3;
//     } else if (Index == 24) {
//         f_SH = y_4_4;
//     }
//     
//     return f_SH;
// }

long double Lebedev_Quadrature_Matrix_SH_Temp(const long double *Matrix, const int &index_SH, const int &degree_SH, const int &Quad_Rule) {
    int order;
    long double *x, *y, *z, *w;
    long double integral_approx, Y_l_m;
    long double norm_f;
    order = Quad_Rule;
    
    x = new long double[2*order*order];
    y = new long double[2*order*order];
    z = new long double[2*order*order];
    w = new long double[2*order*order];
    
    setup_spherical_harmonics_data (order, x, y, z, w);
    
    integral_approx = 0.0;
    
    for ( int m = 0; m < 2*order*order; m++ ) {
        Y_l_m = Cartesian_Spherical_harmonics(x[m], y[m], z[m], index_SH, degree_SH);
        integral_approx = integral_approx + w[m] * Matrix[m] * Y_l_m;
    }
    
//     if (index_SH > 0 && degree_SH != 0) {
//         if (fabs(integral_approx) > 1.0e-8) {
//                 cout << "index_SH = " << index_SH << "   " << "degree_SH = " << degree_SH << endl;
//             for ( int m = 0; m < order; m++ ) {
//                 Y_l_m = Cartesian_Spherical_harmonics(x[m], y[m], z[m], index_SH, degree_SH);
//                 cout << "Matrix = " << Matrix[m] << "   " << "w = " << w[m] << "   " << "   " << "x = " << x[m] << "   " << "y[m] = " << y[m] << "   " << "z[m] = " << z[m] << "   " << "Y_l_m = " << Y_l_m << "   " << "integral = " << integral_approx << endl;
//             }
//         }
//         cout << endl;
//         cout << endl;
//     }
    
//     cout << "index_SH = " << index_SH << "   " << "degree_SH = " << degree_SH << "   " << "integral_approx = " << integral_approx << endl;
    
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] w;
    
    return integral_approx;
}

long double Lebedev_Quadrature_Matrix_SH(const long double *Matrix, const int &index_SH, const int &order, const int &degree) {
    long double *x, *y, *z, *w;
    long double integral_approx, Y_l_m;
    
    x = new long double[2*order*order];
    y = new long double[2*order*order];
    z = new long double[2*order*order];
    w = new long double[2*order*order];
    
    setup_spherical_harmonics_data (order, x, y, z, w);
    
    integral_approx = 0.0;
    
    for ( int m = 0; m < 2*order*order; m++ ) {
        Y_l_m = Cartesian_Spherical_harmonics(x[m], y[m], z[m], index_SH, degree);
        integral_approx = integral_approx + w[m] * Matrix[m] * Y_l_m;
    }
    
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] w;
    
    return integral_approx;
}

long double Lebedev_Quadrature_Orthog_SH(const int &order_SH1, const int &order_SH2, const int &m1, const int &m2, const int &order) {
    long double *x, *y, *z, *w;
    long double integral_approx, Y_l_m_1, Y_l_m_2;
    
    x = new long double[2*order*order];
    y = new long double[2*order*order];
    z = new long double[2*order*order];
    w = new long double[2*order*order];
    
    setup_spherical_harmonics_data (order, x, y, z, w);
    
    integral_approx = 0.0;
    
    for ( int m = 0; m < 2*order*order; m++ ) {
        Y_l_m_1 = Cartesian_Spherical_harmonics(x[m], y[m], z[m], order_SH1, m1);
        Y_l_m_2 = Cartesian_Spherical_harmonics(x[m], y[m], z[m], order_SH2, m2);
//         cout << "Y_l_m_1 = " << Y_l_m_1 << "   " << "Y_l_m_2 = " << Y_l_m_2 << endl;
        integral_approx = integral_approx + w[m] * Y_l_m_1 * Y_l_m_2;
    }
    
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] w;
    
    return integral_approx;
}

void Orthogonality_Test(const int &Quad_Rule) {
    long double delta_ij;
    for (int i_SH1 = 0; i_SH1 < 10; i_SH1++) {
        for (int i_SH2 = 0; i_SH2 < 10; i_SH2++) {
            for (int m_SH1 = -i_SH1; m_SH1 <= i_SH1; m_SH1++) {
                for (int m_SH2 = -i_SH2; m_SH2 <= i_SH2; m_SH2++) {
                    delta_ij = Lebedev_Quadrature_Orthog_SH(i_SH1, i_SH2, m_SH1, m_SH2, Quad_Rule);
                    if ((i_SH1 == i_SH2) && (m_SH1 == m_SH2)) {
                        if (fabs(delta_ij - 1.0) > 1.0e-8) {
                            cout << "i_SH1, i_SH2 = " << i_SH1 << "    " << i_SH2 << "      " << "m_SH1, m_SH2 = " << m_SH1 << "    " << m_SH2 << "     " << "delta_ij = " << delta_ij << endl;
                            exit(0);
                        }
                    } else {
                        if (fabs(delta_ij) > 1.0e-8) {
                            cout << "i_SH1, i_SH2 = " << i_SH1 << "    " << i_SH2 << "      " << "m_SH1, m_SH2 = " << m_SH1 << "    " << m_SH2 << "     " << "delta_ij = " << delta_ij << endl;
                            exit(0);
                        }
                    }
                }
            }
        }
    }
}

long double Test_Spherical_Harmonics(const long double *Vals, const int &index_Cheby, const int &N_points_Cheby, const long double &x, const long double &y, const long double &z, const int &Order_SH) {
    long double f_SH;
    f_SH = 0.0;
    
    for (int i_SH = 0; i_SH < Order_SH; i_SH++) {
        for (int l_SH = -i_SH; l_SH <= i_SH; l_SH++) {
            f_SH += Vals[SH_Linear_Index(i_SH,l_SH)*N_points_Cheby + index_Cheby]*Cartesian_Spherical_harmonics(x,y,z,i_SH,l_SH);
        }
    }
    return f_SH;
}
