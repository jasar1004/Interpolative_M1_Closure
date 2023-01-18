#ifndef _M1_STATE_PARAMETERS_H_INCLUDED
#define _M1_STATE_PARAMETERS_H_INCLUDED

//*********************************************
// M1_State_Parameters.h
//*********************************************

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
#include <cmath>
#include <errno.h>
#include <limits>

#include <mpi.h>

#ifndef _ONE_DIMENSION_QUADRATURE_HPP_INCLUDED
#include "../../Packages/One_Dimensional_Quadratures/quadrule.hpp"
#endif // _ONE_DIMENSION_QUADRATURE_HPP_INCLUDED

#ifndef _CIRCLE_RULE_HPP_INCLUDED
#include "../../Packages/Circle_Quadrature/circle_rule.hpp"
#endif // _CIRCLE_RULE_HPP_INCLUDED


#define PI 3.141592653589793238462643383279502884
#define TOLER 1e-12
#define TOLER_CONSTRAINTS 1e-12

#define PATHVAR "NLOPT_MN_OPTIMIZATION_Path"

#define LEAST_SQUARES_L_CHI2_OFF                          30
#define LEAST_SQUARES_L_CHI2_ON                           31

#define GRAY                 900
#define NON_GRAY             901

#define BOSE_EINSTEIN        1000
#define HYPERBOLIC_LIMIT     1001
#define LOGARITHMIC_LIMIT    1002

#define VAR_E                2000
#define VAR_N1                2001

#define PREDEFINED                3000
#define DISCRETE                  3001

#define ONE_DIMENSIONAL           4000
#define THREE_DIMENSIONAL         4001

#define VAR_THETA                 5000
#define VAR_PHI                   5001

#define FORWARD                   6000
#define BACKWARD                  6001
#define CENTRAL                   6001

#define N1_1_EQUAL_ZERO           7000
#define N1_1_EQUAL_N1             7001
#define FULL_3D                   7001

#define DISCRETE_SET_ENERGY       8000

#define UNIFORM_ALGEBRAIC_MAPPING_N1        9000
#define CHEBYSHEV_ALGEBRAIC_MAPPING_N1      9001

#define SYMMETRY_N1_1                       10000
#define SYMMETRY_N1_2                       10001
#define SYMMETRY_N1_3                       10002

#define CLOSING_FLUX                        100000
#define PARTIAL_MOMENTS                     100001

#define TAYLOR_SERIES_EXPANSION                             1000000
#define DERIVATIVE_BOUNDARY                                 1000001
#define DERIVATIVE_BOUNDARY_RADIATIVE_ENERGY                1000002

#define FREE_STREAMING_LIMIT                 1000010

#define LENGTH_SCALE_DIST_FIT               10000000
#define LENGTH_SCALE_DIST_UNIF              10000001

using namespace std;

extern int DISPLAY_ID;
extern int PRIMARY_ID;

extern long double finite_diff_h;
extern long double tol_grad;

struct M1_State_Param {
    bool display = false;
    bool flag_use_mpi = false;
    int id_proc = 0, num_proc = 1;
    long double theta_val_peak, phi_val_peak;
    static const int NFUNS_PM_3D = 6;
    static int Max_Ent_Solution_Type;
    int index_Moment = 0, index_Lag_i = 0, index_Lag_j = 0, index_Lag_k = 0;
    long double I0 = 0, N1 = 0;
    long double N1_1, N1_2, N1_3;
    long double x_val, y_val, z_val;
    long double ratio_I0;
    int NVARS = 2;
    int order_quad = 0;
    long double *MOM_VEC = NULL;
    long double *MOM_VEC_iso = NULL;
    long double *Sk = NULL;
    long double *Sk_cur = NULL;
    long double *ak = NULL;
    long double *x = NULL;
    long double *tol_x = NULL;
    int index_p = 0;
    int index_a = 0;
    int Index_Partial_Moments_Derivatives = 0;
    int flag_Realizability = 0;
    int Boundary_Point;
    bool flag_Taylor_Series_Expansion = false;
    static int Regime; 
    int Problem_Type = 0, Domain_Type = 0;
    int Node_Dist_E, Node_Dist_f, Node_Dist_Phi;
    int Dimension;
    int N1_1_VAL_DOMAIN;
    
    int finite_difference_type_r_I0;
    int finite_difference_type_N1;
    int flag_finite_difference_r_I0 = false;
    int flag_finite_difference_N1 = false;
    
    long double Omega1_peak, Omega2_peak, Omega3_peak;
    
    // Parameters for finite differencing procedure
    int Finite_Diff_Spectrum_Boundary_Type;
    long double ratio_I0_knot;
    long double f_test_knot;
    long double h_ratio_I0, h_f_test;
    long double x_finite_diff_N1, x_finite_diff_r_I0;
    int index_finite_diff;
    
    // Weight and abscissas for One Dimensional quadrature
    long double *x_quad = NULL;
    long double *w_quad_1D = NULL;
    
    // Weight and abscissas for three-dimensional quadrature
    long double *phi_quad = NULL;
    long double *w_phi_quad = NULL;
    long double *mu_quad = NULL;
    long double *w_mu_quad = NULL;
    
    long double *Omega1 = NULL, *Omega2 = NULL, *Omega3 = NULL;
    long double *w_quad_3D = NULL;
    
    // Parameters for quadrature refinement
    static const int max_refinement_level = 20;
    long double bounds_quad_mu_refin[max_refinement_level+1];
    long double bounds_quad_phi_refin[max_refinement_level+1];
    int refinement_level = 0;
    int N_intervals_mu_quad = 1, N_intervals_phi_quad = 1;
    long double mu_peak, phi_peak;
    
    void Find_Peak_Locations();
    
    void Cartesian_to_Spherical_Coordinates(long double &theta, long double &phi, const long double &Omega1, const long double &Omega2, const long double &Omega3);
    
    void set_bounds_quad_refin(const int &id_refinement_level);
    void set_bounds_quad_refin_mu();
    void set_bounds_quad_refin_phi();
    
    int N_dirs(const int &N_Dirs_mu, const int &N_Dirs_Phi);
    int Id_Angle(const int &Id_angle_mu, const int &Id_angle_Phi);
    
    long double linear_interpolation(const long double &a_new, const long double &b_new, const long double &val_orig, const long double &a_orig, const long double &b_orig);
    
    long double diff_linear_interpolation(const long double &a_new, const long double &b_new, const long double &a_orig, const long double &b_orig);
    
    // void Peak_Location();
    // void Peak_Location(long double &theta_p, long double &phi_p);
    
    void compute_x_quad(const long double &mu_start, const long double &mu_end);
    
    void compute_Omegas();
    
    void compute_Omegas(const long double &mu_start, const long double &mu_end, const long double &phi_start, const long double &phi_end);
    void Set_3D_Quad(void);
    
    // Allocation/deallocation for quadrature in 1D
    void Allocate_Quad_1D();
    void Deallocate_Quad_1D();
    
    // Allocation/deallocation for quadrature in 3D
    void Allocate_Quad_3D();
    void Deallocate_Quad_3D();
    
    void Allocate();
    void Deallocate();
    
    M1_State_Param(const int &NUM_VARS); // Constructor
    
    void copy(const M1_State_Param &M1_State);
    M1_State_Param(const M1_State_Param &M1_State); // Constructor
    
    void Set_Num_Vars();
    
    void set_x(const long double *x_iter);
    void set_Sk(const long double *Sk_iter);
    void set_Sk_cur(const long double *Sk_cur_iter);
    void set_ak(const long double *ak_cur);
};

inline M1_State_Param :: M1_State_Param(const M1_State_Param &M1_State) {
    copy(M1_State);
    Allocate();
    
    switch (Dimension) {
        case ONE_DIMENSIONAL:
            Allocate_Quad_1D();
            break;
        case THREE_DIMENSIONAL:
            Allocate_Quad_3D();
            break;
        default:
            cout << "Domain type not specified for quadrature refinement" << endl;
            exit(0);
            break;
    }
    
    set_bounds_quad_refin(refinement_level);
}

inline void M1_State_Param :: copy(const M1_State_Param &M1_State) {
    NVARS = M1_State.NVARS;
    order_quad = M1_State.order_quad;
    
    ratio_I0 = M1_State.ratio_I0;
    I0 = M1_State.I0;
    N1 = M1_State.N1;
    N1_1 = M1_State.N1_1;
    N1_2 = M1_State.N1_2;
    N1_3 = M1_State.N1_3;
    
    display = M1_State.display;
    flag_use_mpi = M1_State.flag_use_mpi;
    id_proc = M1_State.id_proc;
    num_proc = M1_State.num_proc;
    theta_val_peak = M1_State.theta_val_peak;
    phi_val_peak = M1_State.phi_val_peak;
    
    index_Moment = M1_State.index_Moment;
    index_Lag_i = M1_State.index_Lag_i;
    index_Lag_j = M1_State.index_Lag_j;
    index_Lag_k = M1_State.index_Lag_k;
    
    x_val = M1_State.x_val;
    y_val = M1_State.y_val;
    z_val = M1_State.z_val;
    
    index_p = M1_State.index_p;
    index_a = M1_State.index_a;
    Index_Partial_Moments_Derivatives = M1_State.Index_Partial_Moments_Derivatives;
    flag_Realizability = M1_State.flag_Realizability;
    Boundary_Point = M1_State.Boundary_Point;
    flag_Taylor_Series_Expansion = M1_State.flag_Taylor_Series_Expansion;
    Regime = M1_State.Regime;
    Problem_Type = M1_State.Problem_Type;
    Domain_Type = M1_State.Domain_Type;
    
    Node_Dist_E = M1_State.Node_Dist_E;
    Node_Dist_f = M1_State.Node_Dist_f;
    Node_Dist_Phi = M1_State.Node_Dist_Phi;
    Dimension = M1_State.Dimension;
    N1_1_VAL_DOMAIN = M1_State.N1_1_VAL_DOMAIN;
    
    finite_difference_type_r_I0 = M1_State.finite_difference_type_r_I0;
    finite_difference_type_N1 = M1_State.finite_difference_type_N1;
    flag_finite_difference_r_I0 = M1_State.flag_finite_difference_r_I0;
    flag_finite_difference_N1 = M1_State.flag_finite_difference_N1;
    
    Omega1_peak = M1_State.Omega1_peak;
    Omega2_peak = M1_State.Omega2_peak;
    Omega3_peak = M1_State.Omega3_peak;
    
    Finite_Diff_Spectrum_Boundary_Type = M1_State.Finite_Diff_Spectrum_Boundary_Type;
    ratio_I0_knot = M1_State.ratio_I0_knot;
    f_test_knot = M1_State.f_test_knot;
    h_ratio_I0 = M1_State.h_ratio_I0;
    h_f_test = M1_State.h_f_test;
    x_finite_diff_N1 = M1_State.x_finite_diff_N1;
    x_finite_diff_r_I0 = M1_State.x_finite_diff_r_I0;
    index_finite_diff = M1_State.index_finite_diff;
    
    refinement_level = M1_State.refinement_level;
}

//**************************************************************************************
// This routine computes the location of the potential peak for the M1 distribution
//**************************************************************************************
inline void M1_State_Param :: Find_Peak_Locations() {
    long double Omega1_peak, Omega2_peak, Omega3_peak;
    long double theta_peak_val, phi_peak_val;
    
    // For the M1 distribution, peaks only occur near the free-streaming limit
    // in which case the distribution approaches one that is zero everywhere except
    // along the direction:
    // s = [\omega_1, \omega_2, \omega_3] = [N1_1, N1_2, N1_3]
    if (N1 < 1.0e-8) {
        Omega1_peak = 1.0;
        Omega2_peak = 0.0;
        Omega3_peak = 0.0;
    } else {
        // Normalize so as to have a unit vector
        Omega1_peak = N1_1/N1;
        Omega2_peak = N1_2/N1;
        Omega3_peak = N1_3/N1;
    }
    
    Cartesian_to_Spherical_Coordinates(theta_peak_val, phi_peak_val, Omega1_peak, Omega2_peak, Omega3_peak);
    
    phi_peak = phi_peak_val;
    mu_peak = cos(theta_peak_val);
    
    if (phi_peak < 0.0) {
        phi_peak = 2.0*PI + phi_peak;
    } else if (phi_peak > 2.0*PI) {
        phi_peak = phi_peak - 2.0*PI;
    }
}

inline void M1_State_Param :: Cartesian_to_Spherical_Coordinates(long double &theta, long double &phi, const long double &Omega1, const long double &Omega2, const long double &Omega3) {
    long double norm_Omegas, norm_Omegas_partial;
    long double cos_phi, sin_phi;
    
    norm_Omegas = pow(Omega1, 2) + pow(Omega2, 2) + pow(Omega3, 2);
    norm_Omegas = sqrt(norm_Omegas);
    
    norm_Omegas_partial = pow(Omega2, 2) + pow(Omega3, 2);
    norm_Omegas_partial = sqrt(norm_Omegas_partial);
    
    theta = Omega1/norm_Omegas;
    theta = acos(theta);
    
    if (norm_Omegas_partial < 1.0e-8) {
        cos_phi = 1.0;
        sin_phi = 0.0;
    } else {
        cos_phi = Omega2/norm_Omegas_partial;
        sin_phi = Omega3/norm_Omegas_partial;
    }
    
    phi = acos(cos_phi);
    
    if (sin_phi < 0.0) {
       phi = 2.0*PI - phi; 
    }
}

inline void M1_State_Param :: set_bounds_quad_refin_mu() {
    long double mu_min, mu_max;
    long double temp_val;
    bool flag_refine_block;
    int i_start = 0;
    mu_min = -1.0;
    mu_max = 1.0;
    
    long double bounds_quad_mu_refin_coarse[N_intervals_mu_quad+1];
    
    if (refinement_level == 0) {
        // This corresponds to the original block structure
        bounds_quad_mu_refin[0] = mu_min;
        bounds_quad_mu_refin[1] = mu_max;
        N_intervals_mu_quad = 1;
    } else {
        // This corresponds to the subsequent levels of refinement
        for (int i = 0; i < N_intervals_mu_quad + 1; i++) {
            bounds_quad_mu_refin_coarse[i] = bounds_quad_mu_refin[i];
        }
        
        for (int i = 0; i < N_intervals_mu_quad; i++) {
            mu_min = bounds_quad_mu_refin_coarse[i];
            mu_max = bounds_quad_mu_refin_coarse[i+1];
            
            if (mu_peak >= mu_min && mu_peak <= mu_max) {
                // In this case the peak is within the interval of interest
                flag_refine_block = true;
            } else {
                // In this case the peak is not within the interval of interest
                flag_refine_block = false;
            }
            
            if (flag_refine_block) {
                // Then the interval [mu_min, mu_max] must be refined
                temp_val = (mu_min + mu_max)/2.0;
                bounds_quad_mu_refin[i_start] = mu_min; 
                bounds_quad_mu_refin[i_start+1] = temp_val;
                bounds_quad_mu_refin[i_start+2] = mu_max;
                i_start = i_start + 2;
            } else {
                // Then the interval [mu_min, mu_max] must not be refined
                bounds_quad_mu_refin[i_start] = mu_min;
                bounds_quad_mu_refin[i_start+1] = mu_max;
                i_start = i_start + 1;
            }
        }
        N_intervals_mu_quad = i_start;
    }
}

inline void M1_State_Param :: set_bounds_quad_refin_phi() {
    long double phi_min, phi_max;
    long double temp_val;
    bool flag_refine_block;
    int i_start = 0;
    phi_min = 0.0;
    phi_max = 2.0*PI;
    
    long double bounds_quad_phi_refin_coarse[N_intervals_phi_quad+1];
    
    if (refinement_level == 0) {
        // This corresponds to the original block structure
        bounds_quad_phi_refin[0] = phi_min;
        bounds_quad_phi_refin[1] = phi_max;
        N_intervals_phi_quad = 1;
    } else {
        // This corresponds to the subsequent levels of refinement
        for (int i = 0; i < N_intervals_phi_quad + 1; i++) {
            bounds_quad_phi_refin_coarse[i] = bounds_quad_phi_refin[i];
        }
        
        for (int i = 0; i < N_intervals_phi_quad; i++) {
            phi_min = bounds_quad_phi_refin_coarse[i];
            phi_max = bounds_quad_phi_refin_coarse[i+1];
            
            if (phi_peak >= phi_min && phi_peak <= phi_max) {
                // In this case the peak is within the interval of interest
                flag_refine_block = true;
            } else {
                // In this case the peak is not within the interval of interest
                flag_refine_block = false;
            }
            
            if (flag_refine_block) {
                // Then the interval [phi_min, phi_max] must be refined
                temp_val = (phi_min + phi_max)/2.0;
                bounds_quad_phi_refin[i_start] = phi_min; 
                bounds_quad_phi_refin[i_start+1] = temp_val;
                bounds_quad_phi_refin[i_start+2] = phi_max;
                i_start = i_start + 2;
            } else {
                // Then the interval [phi_min, phi_max] must not be refined
                bounds_quad_phi_refin[i_start] = phi_min;
                bounds_quad_phi_refin[i_start+1] = phi_max;
                i_start = i_start + 1;
            }
        }
        N_intervals_phi_quad = i_start;
    }
}

inline void M1_State_Param :: set_bounds_quad_refin(const int &id_refinement_level) {
    Find_Peak_Locations();
    
    switch (Dimension) {
        case ONE_DIMENSIONAL:
            for (int i = 0; i <= id_refinement_level; i++) {
                refinement_level = i;
                set_bounds_quad_refin_mu();
            }
            break;
        case THREE_DIMENSIONAL:
            for (int i = 0; i <= id_refinement_level; i++) {
                refinement_level = i;
                set_bounds_quad_refin_mu();
                set_bounds_quad_refin_phi();
            }
            break;
        default:
            cout << "Domain type not specified for quadrature refinement" << endl;
            exit(0);
            break;
    }
    
//     long double phi_start, phi_end;
//     long double mu_start, mu_end;
//     if (refinement_level > 6) {
//         cout << "refinement_level = " << refinement_level << "  "  << "id_refinement_level = " << id_refinement_level << endl;
//         
//         cout << "phi_peak = " << phi_peak << "  " << "mu_peak = " << mu_peak<< endl;
//         
//         if (Dimension == THREE_DIMENSIONAL) {
//             for (int i_phi_quad = 0; i_phi_quad < N_intervals_phi_quad; i_phi_quad++) {
//                 phi_start = bounds_quad_phi_refin[i_phi_quad];
//                 phi_end = bounds_quad_phi_refin[i_phi_quad+1];
//                 cout << "i_phi_quad = " << i_phi_quad << "  " << "N_intervals_phi_quad = " << N_intervals_phi_quad << "  " << "phi_start = " << phi_start << "  " << "phi_end = " << phi_end << "  " << endl;
//             }
//         }
//         
//         for (int i_mu_quad = 0; i_mu_quad < N_intervals_mu_quad; i_mu_quad++) {
//             mu_start = bounds_quad_mu_refin[i_mu_quad];
//             mu_end = bounds_quad_mu_refin[i_mu_quad+1];
//             cout << "i_mu_quad = " << i_mu_quad << "  " << "N_intervals_mu_quad = " << N_intervals_mu_quad << "  " << "mu_start = " << mu_start << "  " << "mu_end = " << mu_end << "  " << endl;
//         }
//         exit(0);
//     }
}

inline int M1_State_Param :: N_dirs(const int &N_Dirs_mu, const int &N_Dirs_Phi) {
    int Ndirs;
    Ndirs = N_Dirs_mu*N_Dirs_Phi;
    return Ndirs;
}

inline int M1_State_Param :: Id_Angle(const int &Id_angle_mu, const int &Id_angle_Phi) {
    int index;
    index = Id_angle_Phi*order_quad + Id_angle_mu;
    return index;
}

inline long double M1_State_Param :: linear_interpolation(const long double &a_new, const long double &b_new, const long double &val_orig, const long double &a_orig, const long double &b_orig) {
    long double interp;
    
    interp = a_new + (b_new - a_new)*(val_orig - a_orig)/(b_orig - a_orig);
    
    return interp;
}

inline long double M1_State_Param :: diff_linear_interpolation(const long double &a_new, const long double &b_new, const long double &a_orig, const long double &b_orig) {
    long double dinterp;
        
    dinterp = (b_new - a_new)/(b_orig - a_orig);
    
    return dinterp;
}

// inline void M1_State_Param :: Peak_Location() {
//     Peak_Location(theta_val_peak, phi_val_peak);
//     
//     Omega1_peak = cos(theta_val_peak);
//     Omega2_peak = sin(theta_val_peak)*cos(phi_val_peak);
//     Omega3_peak = sin(theta_val_peak)*sin(phi_val_peak);
// }

// inline void M1_State_Param :: Peak_Location(long double &theta_p, long double &phi_p) {
//     long double norm_f, norm_f_N1_1_0;
//     norm_f = pow(N1_1, 2) + pow(N1_2, 2) + pow(N1_3, 2);
//     norm_f = sqrt(norm_f);
//     
//     norm_f_N1_1_0 = pow(N1_2, 2) + pow(N1_3, 2);
//     norm_f_N1_1_0 = sqrt(norm_f_N1_1_0);
//     
//     if (norm_f > 1.0e-1) {
//         mu_p = N1_1/norm_f;
//         if (norm_f_N1_1_0 > 1.0e-6) {
//             phi_p = acos(N1_2/norm_f_N1_1_0);
//         } else {
//             phi_p = 0.0;
//         }
//     } else {
//         mu_p = 1.0;
//         phi_p = 0.0;
//     }
// }

inline void M1_State_Param :: compute_x_quad(const long double &mu_start, const long double &mu_end) {
    long double mu_temp, w_mu_temp;
    lobatto_compute ( order_quad, mu_quad, w_mu_quad );
    w_mu_temp = diff_linear_interpolation(mu_start, mu_end, -1.0, 1.0);
    for (int Id_mu = 0; Id_mu < order_quad; Id_mu++) {
        mu_temp = linear_interpolation(mu_start, mu_end, mu_quad[Id_mu], -1.0, 1.0);
        x_quad[Id_mu] = mu_temp;
        
        w_quad_1D[Id_mu] = w_mu_temp * w_mu_quad[Id_mu];
        // For the integration over phi in [0, 2 PI]
        w_quad_1D[Id_mu] *= 2.0 * PI;
    }
}

inline void M1_State_Param :: compute_Omegas() {
    int index;
    
    for (int Id_Phi = 0; Id_Phi < 2*order_quad; Id_Phi++) {
        for (int Id_mu = 0; Id_mu < order_quad; Id_mu++) {
            index = Id_Phi*order_quad + Id_mu;
            Omega1[index] = mu_quad[Id_mu];
            Omega2[index] = sqrt(1.0 - pow(mu_quad[Id_mu], 2))*cos(phi_quad[Id_Phi]);
            Omega3[index] = sqrt(1.0 - pow(mu_quad[Id_mu], 2))*sin(phi_quad[Id_Phi]); 
                
            w_quad_3D[index] = w_mu_quad[Id_mu];
            w_quad_3D[index] *= 2.0 * PI * w_phi_quad[Id_Phi];
        }
    }
}

inline void M1_State_Param :: compute_Omegas(const long double &mu_start, const long double &mu_end, const long double &phi_start, const long double &phi_end) {
    int index;
    long double mu_temp, phi_temp, w_mu_temp, w_phi_temp;
    
    for (int Id_Phi = 0; Id_Phi < 2*order_quad; Id_Phi++) {
        for (int Id_mu = 0; Id_mu < order_quad; Id_mu++) {
            index = Id_Phi*order_quad + Id_mu;
            
            mu_temp = linear_interpolation(mu_start, mu_end, mu_quad[Id_mu], -1.0, 1.0);
            w_mu_temp = diff_linear_interpolation(mu_start, mu_end, -1.0, 1.0);
            
            phi_temp = linear_interpolation(phi_start, phi_end, phi_quad[Id_Phi], 0.0, 2.0*PI);
            w_phi_temp = diff_linear_interpolation(phi_start, phi_end, 0.0, 2.0*PI);
            
            Omega1[index] = mu_temp;
            Omega2[index] = sqrt(1.0 - pow(mu_temp, 2))*cos(phi_temp);
            Omega3[index] = sqrt(1.0 - pow(mu_temp, 2))*sin(phi_temp); 
            
            w_quad_3D[index] = w_mu_temp * w_mu_quad[Id_mu];
            w_quad_3D[index] *= 2.0 * PI * w_phi_temp * w_phi_quad[Id_Phi];
        }
    }
}

inline void M1_State_Param :: Set_3D_Quad(void) {
    circle_rule ( 2*order_quad, w_phi_quad, phi_quad );
    // legendre_set ( order_quad, mu_quad, w_mu_quad ); 
    lobatto_compute ( order_quad, mu_quad, w_mu_quad );
}

inline void M1_State_Param :: Allocate_Quad_3D() {
    // Deallocate first
    Deallocate_Quad_3D();
    
    // Weight and abscissas for 3D quadrature
    if (phi_quad == NULL) {
        phi_quad = new long double[2*order_quad];
    }
    
    if (w_phi_quad == NULL) {
        w_phi_quad = new long double[2*order_quad];
    }
    
    if (mu_quad == NULL) {
        mu_quad = new long double[order_quad];
    }
    
    if (w_mu_quad == NULL) {
        w_mu_quad = new long double[order_quad];
    }
    
    if (Omega1 == NULL) {
        Omega1 = new long double[2*order_quad*order_quad];
    }
    
    if (Omega2 == NULL) {
        Omega2 = new long double[2*order_quad*order_quad];
    }
    
    if (Omega3 == NULL) {
        Omega3 = new long double[2*order_quad*order_quad];
    }
    
    if (w_quad_3D == NULL) {
        w_quad_3D = new long double[2*order_quad*order_quad];
    }
    
    circle_rule ( 2*order_quad, w_phi_quad, phi_quad );
    lobatto_compute ( order_quad, mu_quad, w_mu_quad );
    // legendre_set ( order_quad, mu_quad, w_mu_quad ); 
}

inline void M1_State_Param :: Deallocate_Quad_3D() {
    if (phi_quad != NULL) {
        delete[] phi_quad; phi_quad = NULL;
    }
    if (w_phi_quad != NULL) {
        delete[] w_phi_quad; w_phi_quad = NULL;
    }
    if (mu_quad != NULL) {
        delete[] mu_quad; mu_quad = NULL;
    }
    if (w_mu_quad != NULL) {
        delete[] w_mu_quad; w_mu_quad = NULL;
    }
    if (Omega1 != NULL) {
        delete[] Omega1; Omega1 = NULL;
    }
    if (Omega2 != NULL) {
        delete[] Omega2; Omega2 = NULL;
    }
    if (Omega3 != NULL) {
        delete[] Omega3; Omega3 = NULL;
    }
    if (w_quad_3D != NULL) {
        delete[] w_quad_3D; w_quad_3D = NULL;
    }
}

inline void M1_State_Param :: Allocate() {
    // Deallocate first
    Deallocate();
    
    // Allocate now
    if (MOM_VEC == NULL) {
        MOM_VEC = new long double [NVARS];
    }
    if (MOM_VEC_iso == NULL) {
        MOM_VEC_iso = new long double [NVARS];
    }
    if (Sk == NULL) {
        Sk = new long double [NVARS*NVARS];
    }
    if (Sk_cur == NULL) {
        Sk_cur = new long double [NVARS*NVARS];
    }
    if (ak == NULL) {
        ak = new long double [NVARS*NVARS];
    }
    if (x == NULL) {
        x = new long double [NVARS];
    }
    if (tol_x == NULL) {
        tol_x = new long double [NVARS];
    }
}

inline void M1_State_Param :: Allocate_Quad_1D() {
    // Deallocate first
    Deallocate_Quad_1D();
    
    // Weight and abscissas for quadrature
    if (w_mu_quad == NULL) {
        w_mu_quad = new long double[order_quad];
    }
    if (mu_quad == NULL) {
        mu_quad = new long double[order_quad];
    }
    
    if (x_quad == NULL) {
        x_quad = new long double[order_quad];
    }
    
    if (w_quad_1D == NULL) {
        w_quad_1D = new long double[order_quad];
    }
    
    compute_x_quad(-1.0, 1.0);
    
    // In situtations where the domain of integration in split into several subdomains, 
    // gaussian quadratures can be used to compute the integral within each subdomain (even
    // the bounds of integration would not be [-1, 1]). In fact Gaussian quadrature can be 
    // applied to any type of interval domain because they can be thought of as an interpolation
    // of the integrand by a polynomial which is susequently integrated.
    
    // legendre_set ( order_quad, mu_quad, w_mu_quad );
    // clenshaw_curtis_compute ( order_quad, mu_quad, w_mu_quad );
    // lobatto_compute ( order_quad, mu_quad, w_mu_quad );
    // lobatto_set ( order_quad, mu_quad, w_mu_quad );
}

inline void M1_State_Param :: Deallocate() {
    if (MOM_VEC != NULL) {
        delete[] MOM_VEC; MOM_VEC = NULL;
    }
    if (MOM_VEC_iso != NULL) {
        delete[] MOM_VEC_iso; MOM_VEC_iso = NULL;
    }
    if (Sk != NULL) {
        delete[] Sk; Sk = NULL;
    }
    if (Sk_cur != NULL) {
        delete[] Sk_cur; Sk_cur = NULL;
    }
    if (ak != NULL) {
        delete[] ak; ak = NULL;
    }
    if (x != NULL) {
        delete[] x; x = NULL;
    }
    if (tol_x != NULL) {
        delete[] tol_x; tol_x = NULL;
    }
}

inline void M1_State_Param :: Deallocate_Quad_1D() {
    if (w_mu_quad != NULL) {
        delete[] w_mu_quad; w_mu_quad = NULL;
    }
    if (mu_quad != NULL) {
        delete[] mu_quad; mu_quad = NULL;
    }
    if (x_quad != NULL) {
        delete[] x_quad; x_quad = NULL;
    }
    if (w_quad_1D != NULL) {
        delete[] w_quad_1D; w_quad_1D = NULL;
    }
}

inline M1_State_Param :: M1_State_Param(const int &NUM_VARS) { // Constructor
    NVARS = NUM_VARS;
    Allocate();
}

inline void M1_State_Param :: Set_Num_Vars() {
    switch (Dimension) {
        case ONE_DIMENSIONAL:
            NVARS = 2;
            break;
        case THREE_DIMENSIONAL:
            NVARS = 4;
            break;
        default:
            cout << "Dimension type for M1_State not specified !!!!!!" << endl;
            exit(0);
            break;
    }
    
    Allocate();
}

inline void M1_State_Param :: set_x(const long double *x_iter) {
    for (int i = 0; i < NVARS; i++) {
        x[i] = x_iter[i];
    }
}

inline void M1_State_Param :: set_Sk(const long double *Sk_iter) {
    for (int i = 0; i < NVARS*NVARS; i++) {
        Sk[i] = Sk_iter[i];
    }
}

inline void M1_State_Param :: set_Sk_cur(const long double *Sk_cur_iter) {
    for (int i = 0; i < NVARS*NVARS; i++) {
        Sk_cur[i] = Sk_cur_iter[i];
    }
} 

inline void M1_State_Param :: set_ak(const long double *ak_cur) {
    for (int i = 0; i < NVARS*NVARS; i++) {
        ak[i] = ak_cur[i];
    }
}

#endif // _M1_STATE_PARAMETERS_H_INCLUDED
