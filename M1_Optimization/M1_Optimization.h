#ifndef _M1_Model_1D_H_INCLUDED
#define _M1_Model_1D_H_INCLUDED

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
// #include "../../../../../CFFC/nlopt/include/nlopt.hpp"
#include "/home/g/groth/jasar10/CFFC/nlopt/include/nlopt.hpp"
#include <mpi.h>
#include "../../Packages/cubature_1_0_3/cubature.h"
#include "../../Packages/Chebyshev/chebyshev.hpp"

#ifndef _QUADRULE_HPP_INCLUDED
#include "../../Packages/One_Dimensional_Quadratures/quadrule.hpp"
#endif // _QUADRULE_HPP_INCLUDED

#ifndef _SPHERE_LEBEDEV_RULE_HPP_INCLUDED
#include "../../Packages/Lebedev_Quadrature/sphere_lebedev_rule.hpp"
#endif // _SPHERE_LEBEDEV_RULE_HPP_INCLUDED

#ifndef _CIRCLE_RULE_HPP_INCLUDED
#include "../../Packages/Circle_Quadrature/circle_rule.hpp"
#endif // _CIRCLE_RULE_HPP_INCLUDED

#ifndef _M1_STATE_PARAMETERS_H_INCLUDED
#include "./M1_State_Parameters.h"
#endif // _M1_STATE_PARAMETERS_H_INCLUDED

#include "../../Packages/Spline/spline.h"

#ifndef _DIFFER_HPP_INCLUDED
#include "../../Packages/Finite_Difference/differ.hpp"
#endif //_DIFFER_HPP_INCLUDED

#include "../../Packages/Permutations/Permutations_With_Order.h"

using namespace std;
using namespace nlopt;

extern long double I0_val_max;

extern long double f_L_Chi2_vals[10];
extern long double f_L_In_Plus_vals[10];

extern long double L_vals[11], E_vals[11];
extern long double r_l[4];
extern int Lebed_Rule_Set[6];
extern int N_quad_points_Circle_Set[6];

typedef int (*max_ent_obj_grad_hess)(long double *F_obj_grad_hess, const int &NFUN, const int &Index_Angle, void *fdata);

typedef long double (*Mapping_N1) (const long double &, const long double &);

struct record_Ncoeffs {                                                                                                   
    int N_Points_E, N_Points_f, N_pts_Leg, N_pts_quad;
};

struct record_Chi2 {                                                                                                   
    long double I0, N1, Chi2;
    long double dChi2_dN1, dChi2_dI0, dChi2_drI0;
    long double d2Chi2_dN12, d3Chi2_dN13, d2_Chi2_drI0_dN1, d2_Chi2_dI0_dN1;
    long double d3Chi2_drI0_dN1_2, d3Chi2_dI0_dN1_2;
    long double ratio_I0;
};

struct Algebraic_Mapping_N1 {                                                                                                   
    long double f_L_Chi2;                                                                                               
    long double f_L_PM;
    bool flag_Algebraic_Map_N1 = false;
    Mapping_N1 Mapping_L_Chi2;
    Mapping_N1 Inverse_mapping_L_Chi2;
    Mapping_N1 Mapping_L_In_Plus;
    Mapping_N1 Inverse_mapping_L_In_Plus;
};

struct record_Partial_Moments {      
    long double ratio_I0;                                                                                             
    long double I0, N1, N1_1, N1_2, N1_3;
    long double x_val, y_val, z_val;
    long double I0_Plus, N1_1_Plus, N1_2_Plus, N2_11_Plus, N2_12_Plus, N2_22_Plus;
    
    long double dI0_Plus_dI0, dI0_Plus_dnorm_f, dI0_Plus_dN1_1, dI0_Plus_dN1_2, dI0_Plus_dN1_3;
    long double dN1_1_Plus_dI0, dN1_1_Plus_dnorm_f, dN1_1_Plus_dN1_1, dN1_1_Plus_dN1_2, dN1_1_Plus_dN1_3;
    long double dN1_2_Plus_dI0, dN1_2_Plus_dnorm_f, dN1_2_Plus_dN1_1, dN1_2_Plus_dN1_2, dN1_2_Plus_dN1_3;
    long double dN2_11_Plus_dI0, dN2_11_Plus_dnorm_f, dN2_11_Plus_dN1_1, dN2_11_Plus_dN1_2, dN2_11_Plus_dN1_3;
    long double dN2_12_Plus_dI0, dN2_12_Plus_dnorm_f, dN2_12_Plus_dN1_1, dN2_12_Plus_dN1_2, dN2_12_Plus_dN1_3;
    long double dN2_22_Plus_dI0, dN2_22_Plus_dnorm_f, dN2_22_Plus_dN1_1, dN2_22_Plus_dN1_2, dN2_22_Plus_dN1_3;
    
    long double d2I0_Plus_dI0, d2I0_Plus_dnorm_f_2, d2I0_Plus_dN1_1, d2I0_Plus_dN1_2, d2I0_Plus_dN1_3;
    long double d2N1_1_Plus_dI0, d2N1_1_Plus_dnorm_f_2, d2N1_1_Plus_dN1_1, d2N1_1_Plus_dN1_2, d2N1_1_Plus_dN1_3;
    long double d2N1_2_Plus_dI0, d2N1_2_Plus_dnorm_f_2, d2N1_2_Plus_dN1_1, d2N1_2_Plus_dN1_2, d2N1_2_Plus_dN1_3;
    long double d2N2_11_Plus_dI0, d2N2_11_Plus_dnorm_f_2, d2N2_11_Plus_dN1_1, d2N2_11_Plus_dN1_2, d2N2_11_Plus_dN1_3;
    long double d2N2_12_Plus_dI0, d2N2_12_Plus_dnorm_f_2, d2N2_12_Plus_dN1_1, d2N2_12_Plus_dN1_2, d2N2_12_Plus_dN1_3;
    long double d2N2_22_Plus_dI0, d2N2_22_Plus_dnorm_f_2, d2N2_22_Plus_dN1_1, d2N2_22_Plus_dN1_2, d2N2_22_Plus_dN1_3;
};

struct my_constraint_data {
    long double a = 0.0, b = 0.0, c = 0.0, f = 0.0;
    
    long double operator[](int &index) {
        long double temp_val;
        
        switch(index) {
            case 0:
                temp_val = a;
                break;
            case 1:
                temp_val = b;
                break;
            case 2:
                temp_val = c;
                break;
            case 3:
                temp_val = f;
                break;
        }
        return temp_val;
    }
};

struct Finite_Diff_Parameters {
    long double x0;
    long double Chi2_x0;
    long double dChi2_dN1, d2Chi2_dN1, d3Chi2_dN1, d4Chi2_dN1;
    int flag_finite_diff_setup = 0;
    
    long double I0_plus_x0, N1_1_plus_x0, N1_2_plus_x0, N2_11_plus_x0, N2_12_plus_x0, N2_22_plus_x0;
    long double dI0_Plus_dnorm_f, d2_I0_Plus_dnorm_f_2, d3_I0_Plus_dnorm_f_3, d4_I0_Plus_dnorm_f_4;
    long double dN1_1_Plus_dnorm_f, d2_N1_1_Plus_dnorm_f_2, d3_N1_1_Plus_dnorm_f_3, d4_N1_1_Plus_dnorm_f_4;
    long double dN1_2_Plus_dnorm_f, d2_N1_2_Plus_dnorm_f_2, d3_N1_2_Plus_dnorm_f_3, d4_N1_2_Plus_dnorm_f_4;
    long double dN2_11_Plus_dnorm_f, d2_N2_11_Plus_dnorm_f_2, d3_N2_11_Plus_dnorm_f_3, d4_N2_11_Plus_dnorm_f_4;
    long double dN2_12_Plus_dnorm_f, d2_N2_12_Plus_dnorm_f_2, d3_N2_12_Plus_dnorm_f_3, d4_N2_12_Plus_dnorm_f_4;
    long double dN2_22_Plus_dnorm_f, d2_N2_22_Plus_dnorm_f_2, d3_N2_22_Plus_dnorm_f_3, d4_N2_22_Plus_dnorm_f_4;
    
    long double Taylor_Series(const long double &norm_f) {
        long double delta_x;
        long double Chi2_Approx;
        delta_x = norm_f - x0;
        
        Chi2_Approx = Chi2_x0 + delta_x*dChi2_dN1 + pow(delta_x,2)*d2Chi2_dN1/2.0 + pow(delta_x,3)*d3Chi2_dN1/6.0 + pow(delta_x,4)*d4Chi2_dN1/24.0;
        
        return Chi2_Approx;
    }
    
    void Taylor_Series_PM(long double &I0_plus, long double &N1_1_plus, long double &N1_2_plus, long double &N2_11_plus, long double &N2_12_plus, long double &N2_22_plus, const long double &norm_f) {
        long double delta_x;
        delta_x = norm_f - x0;
        
        // Series expansion for I0_plus
        I0_plus = I0_plus_x0 + delta_x*dI0_Plus_dnorm_f + pow(delta_x,2)*d2_I0_Plus_dnorm_f_2/2.0 + pow(delta_x,3)*d3_I0_Plus_dnorm_f_3/6.0 + pow(delta_x,4)*d4_I0_Plus_dnorm_f_4/24.0;
        
        // Series expansion for N1_1_plus
        N1_1_plus = N1_1_plus_x0 + delta_x*dN1_1_Plus_dnorm_f + pow(delta_x,2)*d2_N1_1_Plus_dnorm_f_2/2.0 + pow(delta_x,3)*d3_N1_1_Plus_dnorm_f_3/6.0 + pow(delta_x,4)*d4_N1_1_Plus_dnorm_f_4/24.0;
        
        // Series expansion for N1_2_plus
        N1_2_plus = N1_2_plus_x0 + delta_x*dN1_2_Plus_dnorm_f + pow(delta_x,2)*d2_N1_2_Plus_dnorm_f_2/2.0 + pow(delta_x,3)*d3_N1_2_Plus_dnorm_f_3/6.0 + pow(delta_x,4)*d4_N1_2_Plus_dnorm_f_4/24.0;
        
        // Series expansion for N2_11_plus
        N2_11_plus = N2_11_plus_x0 + delta_x*dN2_11_Plus_dnorm_f + pow(delta_x,2)*d2_N2_11_Plus_dnorm_f_2/2.0 + pow(delta_x,3)*d3_N2_11_Plus_dnorm_f_3/6.0 + pow(delta_x,4)*d4_N2_11_Plus_dnorm_f_4/24.0;
        
        // Series expansion for N2_12_plus
        N2_12_plus = N2_12_plus_x0 + delta_x*dN2_12_Plus_dnorm_f + pow(delta_x,2)*d2_N2_12_Plus_dnorm_f_2/2.0 + pow(delta_x,3)*d3_N2_12_Plus_dnorm_f_3/6.0 + pow(delta_x,4)*d4_N2_12_Plus_dnorm_f_4/24.0;
        
        // Series expansion for N2_22_plus
        N2_22_plus = N2_22_plus_x0 + delta_x*dN2_22_Plus_dnorm_f + pow(delta_x,2)*d2_N2_22_Plus_dnorm_f_2/2.0 + pow(delta_x,3)*d3_N2_22_Plus_dnorm_f_3/6.0 + pow(delta_x,4)*d4_N2_22_Plus_dnorm_f_4/24.0;
    }
    
    long double Taylor_Series_First_Derivative(const long double &norm_f) {
        long double delta_x;
        long double Chi2_Approx;
        delta_x = norm_f - x0;
        
        Chi2_Approx = dChi2_dN1 + delta_x*d2Chi2_dN1 + pow(delta_x,2)*d3Chi2_dN1/2.0 + pow(delta_x,3)*d4Chi2_dN1/6.0;
        
        return Chi2_Approx;
    }
};

struct Mobius_Scale_Parameters {
    int N_pts_Mob_Scale;
    int Length_Scale_Dist_Type;
    int Least_Squares_L_Chi2_Mode;
    long double *Coefficients_Mobius_Scale_Fit = NULL;
    
    // Constructor
    Mobius_Scale_Parameters(const int &Num_points_Mob_Scale, const long double *Coeffs_Mobius_Scale) {
        N_pts_Mob_Scale = Num_points_Mob_Scale;
        allocate();
        Set_Coefficients_Mob_Scale_Fit(Coeffs_Mobius_Scale);
    }
    
    // Destructor
    ~Mobius_Scale_Parameters() {
        deallocate();
    }
    
    void Set_Coefficients_Mob_Scale_Fit(const long double *Coeffs_Mobius_Scale) {
        for (int i = 0; i < N_pts_Mob_Scale; i++) {
            Coefficients_Mobius_Scale_Fit[i] = Coeffs_Mobius_Scale[i];
        }
    }
    
    void allocate() {
        deallocate();
        Coefficients_Mobius_Scale_Fit = new long double[N_pts_Mob_Scale];
    }
    
    void deallocate() {
        if (Coefficients_Mobius_Scale_Fit != NULL) {
            delete Coefficients_Mobius_Scale_Fit; Coefficients_Mobius_Scale_Fit = NULL;
        }
    }
    
    long double Evaluate_Length_Scale(const long double &norm_f) {
        long double Length_Scale;
        long double norm_f_2 = norm_f*norm_f;
        int index = N_pts_Mob_Scale - 1;
        
        long double poly_val;
        
        switch (Least_Squares_L_Chi2_Mode) {
            case LEAST_SQUARES_L_CHI2_ON:
                Length_Scale = 0.0;
                for (int i_fit_f = 0; i_fit_f < N_pts_Mob_Scale; i_fit_f++) {
                    poly_val = Chebyshev_Polynomial_Basis(norm_f, 2*i_fit_f);
                    Length_Scale += Coefficients_Mobius_Scale_Fit[i_fit_f] * poly_val;
                }
                break;
            case LEAST_SQUARES_L_CHI2_OFF:
                for (int i_fit_f = N_pts_Mob_Scale - 1; i_fit_f >= 0; i_fit_f--) {
                    if (i_fit_f == N_pts_Mob_Scale - 1) {
                        Length_Scale = Coefficients_Mobius_Scale_Fit[i_fit_f];
                    } else {
                        Length_Scale = Coefficients_Mobius_Scale_Fit[i_fit_f] + Length_Scale*norm_f_2;
                    }
                    index--;
                }
                break;
            default:
                cout << "Invalid value for Least_Squares_L_Chi2_Mode !!!!!!!!!!!" << endl;
                exit(0);
                break;
        };
        
        if (Length_Scale_Dist_Type == LENGTH_SCALE_DIST_FIT) {
            Length_Scale = exp(Length_Scale);
        } else if (Length_Scale_Dist_Type != LENGTH_SCALE_DIST_UNIF) {
            cout << "Length scale distribution type not specified!!!!!!!!!!!!!!!" << endl;
            exit(0);
        }
        
        // cout << "N_pts_Mob_Scale = " << N_pts_Mob_Scale << "  " << "norm_f = " << norm_f << "  " << "Length_Scale = " << Length_Scale << endl;
        
        return Length_Scale;
    }
    
//     long double Evaluate_Length_Scale_r_N1(const long double &r_N1) {
//         long double Length_Scale;
//         int index = N_pts_Mob_Scale - 1;
//         
//         for (int i_fit_f = N_pts_Mob_Scale - 1; i_fit_f >= 0; i_fit_f--) {
//             if (i_fit_f == N_pts_Mob_Scale - 1) {
//                 Length_Scale = Coefficients_Mobius_Scale_Fit[i_fit_f];
//             } else {
//                 Length_Scale = Coefficients_Mobius_Scale_Fit[i_fit_f] + Length_Scale*r_N1;
//             }
//             index--;
//         }
//         
//         if (N_pts_Mob_Scale > 1) {
//             Length_Scale = exp(Length_Scale);
//         }
//         
//         return Length_Scale;
//     }
};

struct M1_Var_Num_Points{
    int E = 0, f = 0;
    int order_quad = 0;
    int N_pts_Leg = 0;
    int N_pts_quad = 0;
    int theta = 0, phi = 0;
    int Symmetry_Type = 0;
    int Maximum_Entropy_Solution_Regime = 0;
    int Length_Scale_Dist_Type = 0;
    int Least_Squares_L_Chi2_Mode = LEAST_SQUARES_L_CHI2_OFF;
};

struct record_Lag_Mult {                                            
    long double ratio_I0, I0, N1;
    long double x0, x1;
    int Regime;
    
    long double &operator[](int &index) {
        static long double temp_val;
        
        switch(index) {
            case 0:
                temp_val = x0;
                break;
            case 1:
                temp_val = x1;
                break;
        }
        return temp_val;
    }
    
    void reset() {
        x0 = 0.0; x1 = 0.0;
    }
};

void xyz_to_tp_mine ( long double x, long double y, long double z, long double *t, long double *p );
void Rotate_Frame_Vector(long double *x, const int &NVARS, const long double &x_Lebed, const long double &y_Lebed, const long double &z_Lebed);
void Rotate_Frame_Vector_back(long double *x, const int &NVARS, const long double &x_Lebed, const long double &y_Lebed, const long double &z_Lebed);

// Function to write a record at a specific position
template<class record_Vals>
ostream& write(ostream& ios, const record_Vals &rec) {
    ios.clear(); // clear any errors
    ios.write(reinterpret_cast<const char*>(&rec), sizeof(rec)); //write the record to file
    return ios; // return the stream (for easy error detection/chaining)
}

// Function to write a record at a specific position
template<class record_Vals>
iostream& read(iostream& ios, record_Vals &rec, int &id) {
 ios.clear(); // clear any errors
 ios.seekg(sizeof(record_Vals) * id); // move to record's position
 ios.read(reinterpret_cast<char *>(&rec), sizeof(rec)); //read record from file
 return ios; // return the stream (for easy error detection/chaining)
}

// Function to write a record at a specific position
template<class record_Vals_1, class record_Vals_2>
iostream& read(iostream& ios, record_Vals_1 &rec, int &id_1, int &id_2) {
 ios.clear(); // clear any errors
 ios.seekg(sizeof(record_Vals_1) * id_1 + sizeof(record_Vals_2) * id_2); // move to record's position
 ios.read(reinterpret_cast<char *>(&rec), sizeof(rec)); //read record from file
 return ios; // return the stream (for easy error detection/chaining)
}

long double roundval( const long double &val );

long double roundn( const long double &val, const int &n );

long double Cartesian_to_spherical_Coordinates_Conversion_Jacobian(const int &i, const int &j, const long double &r, const long double &cos_theta, const long double &sin_theta, const long double &cos_phi, const long double &sin_phi);

long double Cartesian_to_spherical_Coordinates_Conversion_Jacobian_Second_Derivatives(const int &i, const int &j, const long double &r, const long double &cos_theta, const long double &sin_theta, const long double &cos_phi, const long double &sin_phi);

long double Cartesian_to_spherical_Coordinates_Jacobian(const long double &norm_f, const long double &x_val, const long double &y_val, const long double &z_val, const long double &dIn_Plus_dN1_1, const long double &dIn_Plus_dN1_2, const long double &dIn_Plus_dN1_3, const int &VAR_NUM);
    
long double Cartesian_to_spherical_Coordinates_Jacobian_Second_Derivatives(const long double &norm_f, const long double &x_val, const long double &y_val, const long double &z_val, const long double &dIn_Plus_dN1_1, const long double &dIn_Plus_dN1_2, const long double &dIn_Plus_dN1_3, const int &VAR_NUM);

long double Heaviside(const long double &x);
 
void Modified_Gram_Schmidt_Factorization(const int &n, const long double *A, long double *Q, long double *R);

long double Rotation_Matrix_X_axis(const int &i, const int &j, const long double &cos_angle, const long double &sin_angle);

long double Rotation_Matrix_Y_axis(const int &i, const int &j, const long double &cos_angle, const long double &sin_angle);

long double Rotation_Matrix_Z_axis(const int &i, const int &j, const long double &cos_angle, const long double &sin_angle);

void Rotate_Frame(long double *x, const int &NVARS, const long double &x_Lebed, const long double &y_Lebed, const long double &z_Lebed);

long double Inverse_Mobius_Transformation(const long double &ratio, const long double &Length_Scale);
 
long double Mobius_Transformation(const long double &Moment, const long double &Length_Scale);

long double d_ratio_E_d_I0(long double &Moment, const long double &Length_Scale);

long double d_ratio_E_d_L_Chi2(long double &Moment, const long double &Length_Scale);
 
long double e_skewed(int &i, int &Np);

long double f_skewed(int &i, int &Np);

void setup_chebyshev_nodes_data(const int &N_points_cheby, long double *x_cheby);

void setup_spherical_harmonics_data (const int &N_points_Leg, long double *x_SH, long double *y_SH, long double *z_SH, long double *w_quad);

void generate_polynomials_1D_M1_2D(long double* poly,const long double &mu, const long double &phi);
 
void generate_polynomials_3D(long double* poly, const int &Id_angle, const M1_State_Param &M1_State);

void generate_polynomials_3D(long double* poly, const long double &Omega1, const long double &Omega2, const long double &Omega3);

long double generate_partial_moments_monomials_basis(const long double &Omega1, const long double &Omega2, const long double &Omega3, const int &index);

long double generate_partial_moments_monomials_basis(const int &Id_angle, const M1_State_Param &M1_State);

long double sum_Lagrange_M1_2D(const long double *x, const long double &mu, const long double &phi);
 
void Distribution_M1_BE_2D(int &NDIM, const long double *Int_Vars, int &NFUN, long double *DISTRIBUTION, const long double *x, int &NVARS_M1);

void Distribution_M1_HL_2D(int &NDIM, const long double *Int_Vars, int &NFUN, long double *DISTRIBUTION, const long double *x, int &NVARS_M1);
 
void Distribution_M1_LL_2D(int &NDIM, const long double *Int_Vars, int &NFUN, long double *DISTRIBUTION, const long double *x, int &NVARS_M1);

long double Map_Jacob(const M1_State_Param &M1_State);

void add_constraints(nlopt_opt &opt, my_constraint_data *data, M1_State_Param &M1_State);

long double myconstraint(unsigned n, const long double *x, const long double *Sk, long double *grad, void *data);

long double sum_Lag_Mom(const long double *x, const int &NVARS, const long double* MOM_VEC);

long double sum_Lagrange(const long double *x, const long double *poly_Sk, const int &NVARS);

void generate_polynomials_1D(long double* poly, const int &Id_angle, const M1_State_Param &M1_State);

void One_Dimensional_Quadrature(long double *Vals, const int &NFUN, max_ent_obj_grad_hess func, void *fdata);

void One_Dimensional_Quadrature_Partial_Moments(long double *Vals, const int &NFUN, max_ent_obj_grad_hess func, void *fdata);

void Quadrature_3D_Func(long double *Vals, const int &NFUN, max_ent_obj_grad_hess func, void *fdata);

void Partial_Quadrature_3D_Func(long double *Vals, const int &NFUN, max_ent_obj_grad_hess func, void *fdata);

int Check_Realizability(M1_State_Param *M1_State);

int obj_function(unsigned NDIM, const long double *Int_Vars, void *fdata, unsigned NFUN, long double *F_obj_grad);
 
int gradient(unsigned NDIM, const long double *Int_Vars, void *fdata, unsigned NFUN, long double *grad);
 
int Chi2_M1_1D(long double *MOMENTS, const int &NFUN, const int &Id_angle, void *fdata);

int Angular_Moments_First_Derivatives(long double *MOMENTS, const int &NFUN, const int &Id_angle, void *fdata);

int Angular_Moments_Second_Derivatives(long double *MOMENTS, const int &NFUN, const int &Id_angle, void *fdata);

int Angular_Moments_Third_Derivatives(long double *MOMENTS, const int &NFUN, const int &Id_angle, void *fdata);

void transpose_matrix(const int &nl, const int &nc, long double *A_transpose, const long double *A);

void Backward_Substitution(const int &nl, const int &nc, long double *A, long double *Q, const long double *R);

int H_ONE_Matrix_Entries(long double *H_ONE, const int &NFUN, const int &Id_angle, void *fdata);

int H_TWO_Matrix_Entries(long double *F_Hess, const int &NFUN, const int &Id_angle, void *fdata);

int H_THREE_Matrix_Entries(long double *H_THREE, const int &NFUN, const int &Id_angle, void *fdata);

int Partial_Moments_3D(long double *Partial_Moms, const int &NFUN, const int &Id_angle, void *fdata);

void Partial_Moments_M1_NON_GRAY_3D(long double *PARTIAL_MOMS, const int &NFUN, const long double *x, const long double *Sk, M1_State_Param &M1_State);
 
long double myfunc(unsigned n, const long double *x, const long double *Sk, long double *grad, void *my_func_data);

void Compute_H_ONE_Matrix(long double *H_ONE, const long double *x, const long double *Sk, M1_State_Param &M1_State);

void Compute_H_TWO_Matrix(long double *H_TWO, const long double *x, const long double *Sk, M1_State_Param &M1_State);

void Compute_H_THREE_Matrix(long double *H_THREE, const long double *x, const long double *Sk, M1_State_Param &M1_State);

void Compute_A_ONE_Matrix(long double *H_ONE, const long double *x, const long double *Sk, M1_State_Param &M1_State);

void Compute_A_TWO_Matrix(long double *H_TWO, const long double *x, const long double *Sk, M1_State_Param &M1_State);

void Compute_A_THREE_Matrix(long double *H_THREE, const long double *x, const long double *Sk, M1_State_Param &M1_State);

void Calculate_Inverse_Hessian_Matrix(long double *Q_data, const long double *x, const long double *Sk, M1_State_Param &M1_State);

void Q_func_Gram_Schmidt(unsigned n, const long double *x, const long double *Sk, long double *Q_data, void *my_func_data, void *my_func_MN_State);

long double dIn_dInm1(const int &i, const int &j, const long double &I0, const long double &N1_1, const long double &N1_2, const long double &N1_3);

long double Calculate_Angular_Moments_Second_Derivatives(long double *dChi2dlambda, long double *A_ONE, const int &j, const int &k, M1_State_Param &M1_State);

long double Calculate_Angular_Moments_Third_Derivatives(long double *dChi2dlambda, long double *A_ONE, const int &i, const int &j, const int &k, M1_State_Param &M1_State);

void Calculate_Eddington_Factor_Derivatives(record_Chi2 &rec_Chi2, const long double *x, const long double *Sk, M1_State_Param &M1_State);

void Calculate_Partial_Moments_Derivatives(record_Partial_Moments &rec_PM, const long double *x, const long double *Sk, M1_State_Param &M1_State);

void Setup_Taylor_Series_Coefficients_PM(Finite_Diff_Parameters &Finite_Diff_PM, M1_State_Param &M1_State, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, long double *Sk_final);

void Calculate_Partial_Moments_Derivatives_Boundary(record_Partial_Moments &rec_PM, M1_State_Param &M1_State, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, const long double *x_SH, const long double *y_SH, const long double *z_SH, long double *Sk_final);

void Regularize_Moment(const long double r_l, M1_State_Param *M1_State);

void Set_Regime(M1_State_Param *M1_State);

void Set_Closure_Free_Streaming(record_Chi2 *rec_Chi2_local, const int &id_count, M1_State_Param &M1_State);

void Set_Partial_Moments_Free_Streaming(record_Partial_Moments *rec_PM_local, const int &id_count, M1_State_Param &M1_State);

void Compute_Eddington_Factor_Derivatives(record_Chi2 *rec_Chi2_local, const int &index_count, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, M1_State_Param &M1_State, const long double *x, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, long double *Sk_final);

void Compute_Eddington_Factor(record_Chi2 *rec_Chi2_local, const int &index_count, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, M1_State_Param &M1_State, const long double *x, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, long double *Sk_final);

void Compute_Partial_Moments_Derivatives(record_Partial_Moments *rec_PM_local, const int &index_count, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, M1_State_Param &M1_State, const long double *x, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, const long double *x_SH, const long double *y_SH, const long double *z_SH, long double *Sk_final);

void Compute_Partial_Moments(record_Partial_Moments *rec_PM_local, const int &index_count, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, M1_State_Param &M1_State, const long double *x, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, const long double *x_SH, const long double *y_SH, const long double *z_SH, long double *Sk_final);

void Set_Moments(M1_State_Param *M1_State, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, Mobius_Scale_Parameters *Mobius_Scale_Params, const Algebraic_Mapping_N1 *struct_map_N1);

void Set_Moments_3D(M1_State_Param *M1_State, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed, Mobius_Scale_Parameters *Mobius_Scale_Params, const Algebraic_Mapping_N1 *struct_map_N1);

void Set_Initial_Guess_And_Tol(M1_State_Param &M1_State);

void Eddington_Factor(long double &Chi_2, const long double *x, const long double *Sk, M1_State_Param &M1_State);

void Setup_Taylor_Series_Coefficients(Finite_Diff_Parameters &Finite_Diff_Chi2, M1_State_Param &M1_State, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, long double *Sk_final);

void Calculate_Eddington_Factor_Derivatives_Boundary(record_Chi2 &rec_Chi2, M1_State_Param &M1_State, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, Mobius_Scale_Parameters *Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, long double *Sk_final);

void Calculate_Eddington_Factor_Derivatives_Boundary_Spectrum_Energy(record_Chi2 &rec_Chi2, M1_State_Param &M1_State, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, Mobius_Scale_Parameters *Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, long double *Sk_final);

void NLOPT_Optimization_Algorithm_NON_GRAY_M1(record_Chi2 *rec_Chi2_local, record_Partial_Moments *rec_PM_array, const int &index_count, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, const int &id_proc, M1_State_Param &M1_State, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, long double *Sk_final);

void NLOPT_Optimization_Algorithm_NON_GRAY_M1_3D(record_Chi2 *rec_Chi2_local, record_Partial_Moments *rec_Partial_Moments, const int &index_count, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, const int &id_proc, M1_State_Param &M1_State, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, long double *Sk_final, const long double *x_Lebed, const long double *y_Lebed, const long double *z_Lebed);

void OPTIM_NON_GRAY_M1_Array(record_Chi2 *rec_Chi2_global, record_Partial_Moments *rec_PM_global, const int &Max_Ent_Solution_Type, const M1_Var_Num_Points *num_points, const int &Problem_Type, const int &Dimension, const int &Node_Distribution_E, const int &Node_Distribution_f, const int &id_proc, const int &num_proc, const int &N_pts_Mob_Scale, const long double *Coefficients_Fit_Mob_Scale, const Algebraic_Mapping_N1 *struc_map_N1, const bool &display, const bool flag_use_mpi);

void OPTIM_NON_GRAY_M1(const M1_Var_Num_Points *num_points, const int &Max_Ent_Solution_Type, const int &Problem_Type, const int &Dimension, const int &Node_Distribution_E, const int &Node_Distribution_f, const int &id_proc, const int &num_proc, const int &N_pts_Mob_Scale, const long double *Coefficients_Mobius_Scale_Fit, fstream &in_Chi2_out, const Algebraic_Mapping_N1 *struct_map_N1, const bool &display, const bool flag_use_mpi);

void OPTIM_NON_GRAY_M1_Fixed_N1(record_Chi2 *rec_Chi2_array, record_Partial_Moments *rec_PM_array, const int &Max_Ent_Solution_Type, const M1_Var_Num_Points *num_points, const int &Problem_Type, const int &Dimension, const int &Node_Distribution_E, const int &Node_Distribution_f, const int &id_proc, const int &num_proc, const int &N_pts_Mob_Scale, const long double *Coefficients_Fit_Mob_Scale, const int &index_N1, const Algebraic_Mapping_N1 *struc_map_N1, const bool &display);

void OPTIM_NON_GRAY_M1_1D_BCs(const M1_Var_Num_Points *num_points, const int &Problem_Type, const int &Dimension, const int &Node_Distribution_E, const int &Node_Distribution_f, const int &id_proc, const int &N_pts_Mob_Scale, const long double *Coefficients_Mobius_Scale_Fit, fstream &in_Chi2_out);

void OPTIM_NON_GRAY_M1_3D_Array(record_Chi2 *rec_Chi2_global, record_Partial_Moments *rec_PM_global, const M1_Var_Num_Points *num_points, M1_State_Param &M1_State, const int &N_pts_Mob_Scale, const long double *Coefficients_Mobius_Scale_Fit, const Algebraic_Mapping_N1 *struc_map_N1);

void OPTIM_NON_GRAY_M1_3D(const M1_Var_Num_Points *num_points, M1_State_Param &M1_State, const int &N_pts_Mob_Scale, const long double *Coefficients_Fit_Mob_Scale, fstream &in_Partial_Moments_out, const Algebraic_Mapping_N1 *struc_map_N1);

#endif // _M1_Model_1D_H_INCLUDED
