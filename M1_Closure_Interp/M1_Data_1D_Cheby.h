#ifndef _M1_DATA_1D_CHEBY_H_INCLUDED
#define _M1_DATA_1D_CHEBY_H_INCLUDED

/*******************************************************************
  File: M1_Data_1D_Cheby.h

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
#include <mpi.h>
#include <unistd.h>

#ifndef _M1_Model_1D_H_INCLUDED
#include "../M1_Optimization/M1_Optimization.h"
#endif // _M1_Model_1D_H_INCLUDED

#include "../M1_Optimization/M1_Model_1D_Utilities.h"

#ifndef _EIGENSTRUCTURE_3D_H_INCLUDED
#include "./Eigenstructure_3D.h"
#endif // _EIGENSTRUCTURE_3D_H_INCLUDED

#define INTERPOLANT_UNIFORM_DISTRIBUTION                  10
#define INTERPOLANT_CHEBYSHEV_DISTRIBUTION                11

#define IMPLEMENTATION_DERIVATIVES_N1                     20
#define IMPLEMENTATION_DERIVATIVES_N1_DERIVATIVES_R_I0    21

#define POLYNOMIAL_TYPE                                   40
#define RATIONAL_TYPE                                     41

using namespace nlopt;
using namespace std;

class M1_Data_1D_Cheby;
struct M1_1D_Data_Pointer;

long double dmapping_L_Chi2_dnorm_f_2(const long double &L, const long double &norm_f);
long double dmapping_L_Chi2_d_f_L_Chi_2(const long double &L, const long double &norm_f);
long double mapping_L_Chi2(const long double &L, const long double &norm_f);
long double Inverse_mapping_L_Chi2(const long double &L, const long double &r_N1);

// Routines for use with optimization algorithm
long double myfunc_Least_Squares(unsigned n, const long double *x, long double *grad, void *my_func_data);
void add_constraints(nlopt_opt &opt, my_constraint_data *data, M1_1D_Data_Pointer *data_realiz);
long double myconstraint(unsigned n, const long double *x, long double *grad, void *data);
long double myconstraint_realizability_Lower_Bound(unsigned n, const long double *x, long double *grad, void *data);
long double myconstraint_realizability_Upper_Bound(unsigned n, const long double *x, long double *grad, void *data);
long double myconstraint_Hyperbolicity_Lambda_s(unsigned n, const long double *x, long double *grad, void *data);

long double myfunc_Least_Squares_f_L_Chi2(unsigned n, const long double *x, long double *grad, void *my_func_data);
void add_constraints_f_L_Chi2(nlopt_opt &opt, my_constraint_data *data, M1_1D_Data_Pointer *data_realiz);
long double myconstraint_f_L_Chi2(unsigned n, const long double *x, long double *grad, void *data);
long double myconstraint_f_L_Chi2_realizability_Lower_Bound(unsigned n, const long double *x, long double *grad, void *data);
long double myconstraint_f_L_Chi2_realizability_Upper_Bound(unsigned n, const long double *x, long double *grad, void *data);
long double myconstraint_f_L_Chi2_Hyperbolicity_Lambda_s(unsigned n, const long double *x, long double *grad, void *data);

long double myfunc_Nested_Least_Squares_L_Chi2(unsigned n, const long double *x, long double *grad, void *my_func_data);
void add_constraints_Nested_Least_Squares_L_Chi2(nlopt_opt &opt, my_constraint_data *data, M1_1D_Data_Pointer *data_realiz);
long double myconstraint_Nested_Least_Squares_L_Chi2_realizability_Lower_Bound(unsigned n, const long double *x, long double *grad, void *data);
long double myconstraint_Nested_Least_Squares_L_Chi2_realizability_Upper_Bound(unsigned n, const long double *x, long double *grad, void *data);
long double myconstraint_Nested_Least_Squares_L_Chi2_Hyperbolicity_Lambda_s(unsigned n, const long double *x, long double *grad, void *data);
long double myconstraint_Nested_Least_Squares_L_Chi2_Speed_Limit(unsigned n, const long double *x, long double *grad, void *data);
long double myconstraint_Nested_Least_Squares_L_Chi2_Hyperbolicity_Lambda_s_2D(unsigned n, const long double *x, long double *grad, void *data);

struct M1_1D_Data_Pointer {                                                                                                   
    M1_Data_1D_Cheby *M1_Data_Uniform_BE;                                                                                   
    M1_Data_1D_Cheby *M1_Data_Uniform_HL;                                                                                   
    M1_Data_1D_Cheby *M1_Data_Uniform_LL;
    M1_Data_1D_Cheby *M1_Data_BE;
    M1_Data_1D_Cheby *M1_Data_HL;
    M1_Data_1D_Cheby *M1_Data_LL;
};

class M1_Data_1D_Cheby {
public:
    char path_out[256], prefix[256], extension[256];
    fstream in_out_Optim;
    
    record_Ncoeffs rec_Ncoeffs;
    record_Chi2 rec_Chi2;
    
    static int Problem_Type;
    static bool flag_test_Implementation;
    static int Implementation_type;
    static int Least_Squares_L_Chi2_Mode;
    static int Weighting_Function_fchi2_Interpolation_Type;
    
    int id_proc;
    int num_proc;
    
    int Node_Distribution_E;
    int Node_Distribution_N1;
    
    long double L_inf_Norm;
    long double L2_Norm_Chi2;
    
    long double L_Chi2_0 = 1.0;
    
    int index_N1;
    
    int N_Points_E;
    int N_Points_f;
    int N_pts_Mob_Scale;
    
    int N_pts_L_Chi2_LS = 1;
    
    int index_Ncoeffs;
    int index_rec_Chi2;
    
    long double *E_NON_GRAY = NULL;
    long double *ratio_I0_NON_GRAY = NULL;
    long double *N1_1_NON_GRAY = NULL;
    
    long double *Coefficients_Vander_Matrix_2vars = NULL;
    
    long double *Chi2_NON_GRAY = NULL;
    
    long double *dChi2_drI0 = NULL;
    long double *dChi2_dI0 = NULL;
    long double *dChi2_dN1 = NULL;
    
    long double *d2Chi2_dN12 = NULL;
    long double *d3Chi2_dN13 = NULL;
    long double *d2Chi2_drI0_dN1 = NULL;
    long double *d3Chi2_drI0_dN1_2 = NULL;
    long double *d2Chi2_dI0_dN1 = NULL;
    long double *d3Chi2_dI0_dN1_2 = NULL;
    
    long double *f_Chi2_NON_GRAY = NULL;
    long double *error_fit = NULL;
    
    long double *Coefficient_Matrix_Fit_Chi2 = NULL;
    long double *Coefficient_Matrix_Fit_Mob_Scale = NULL;
    
    long double *Coefficient_Matrix_Fit_Chi2_HL = NULL;
    long double *Coefficient_Matrix_Fit_Chi2_LL = NULL;
    
    // Constructor
    M1_Data_1D_Cheby(const int &Num_Coeffs_E, const int &Num_Coeffs_f, int Num_pts_Mob_Scale = 1, int Num_pts_L_Chi2_LS = 1) {
        N_Points_E = Num_Coeffs_E;
        N_Points_f = Num_Coeffs_f;
        N_pts_Mob_Scale = Num_pts_Mob_Scale;
        N_pts_L_Chi2_LS = Num_pts_L_Chi2_LS;
        allocate();
    }
    
    // Destructor
    ~M1_Data_1D_Cheby() {
        deallocate();
    }
    
    void Copy_to(M1_Data_1D_Cheby *New_M1_Data_1D_Cheby);
    
    void allocate();
    
    void deallocate();
    
    void Allocate_Coefficients_L_Chi2();
    
    void Deallocate_Coefficients_L_Chi2();
    
    void OpenInputFile(char *filename);
    
    void ReadInputData();
    
    void CloseInputFile();
    
    void Check_Hyperbolicity_Max_Ent_Solutions(const int &Type);
    
    void SetupInterpolant_Values_HL_LL();
    
    void SetupInterpolant_Values_BE(M1_Data_1D_Cheby &M1_Data_1D_Cheby_HL, M1_Data_1D_Cheby &M1_Data_1D_Cheby_LL);
    
    void SetupInterpolant_Values_BE_Impl_Derivs_N1(M1_Data_1D_Cheby &M1_Data_1D_Cheby_HL, M1_Data_1D_Cheby &M1_Data_1D_Cheby_LL);
    
    void SetupInterpolant_Values_BE_Impl_Derivs_N1_Derivs_r_I0(M1_Data_1D_Cheby &M1_Data_1D_Cheby_HL, M1_Data_1D_Cheby &M1_Data_1D_Cheby_LL);
    
    void SetupInterpolant_Values_Fixed_Mob_Scale_Dist(const long double *index_Mobius_Perms);
    
    void Setup_Vandermonde_Matrix_2vars();
    
    void Vandermonde_Interpolation_2vars();
    
    void Setup_Coefficients_HL_LL(M1_Data_1D_Cheby &Chi2_M1_1D_HL, M1_Data_1D_Cheby &Chi2_M1_1D_LL);
    
    void Polynomial_Interpolation_HL_LL(const int &Maximum_Entropy_Solution_Regime);
    
    void Polynomial_Interpolation_BE(M1_Data_1D_Cheby &M1_Data_1D_Cheby_HL, M1_Data_1D_Cheby &M1_Data_1D_Cheby_LL);
    
    long double Evaluate_h_Chi2(const long double &ratio_E, const long double &norm_f);
    
    long double Evaluate_dh_Chi2_dratio_E(const long double &ratio_E, const long double &norm_f);
    
    long double Evaluate_dh_Chi2_d_sqr_N1(const long double &ratio_E, const long double &norm_f);
    
    long double Evaluate_theta_Chi2(const long double &ratio_E, const long double &norm_f);
    
    long double Evaluate_dtheta_Chi2_dratio_E(const long double &ratio_E, const long double &norm_f);
    
    long double Evaluate_dtheta_Chi2_d_sqr_N1(const long double &ratio_E, const long double &norm_f);
    
    long double Evaluate_g_Chi2_HL_LL(const long double &ratio_E, const long double &norm_f, const int &Regime);
    
    long double Evaluate_dg_Chi2_HL_LL_d_sqr_N1(const long double &ratio_E, const long double &norm_f, const int &Regime);
    
    long double Evaluate_g_Chi2(const long double &ratio_E, const long double &norm_f);
    
    long double Evaluate_dg_Chi2_dratio_E(const long double &ratio_E, const long double &norm_f);
    
    long double Evaluate_dg_Chi2_d_sqr_N1(const long double &ratio_E, const long double &norm_f);
    
    long double Evaluate_f_Chi2(const long double &ratio_E, const long double &norm_f);
    
    long double Evaluate_df_Chi2_dratio_E(const long double &ratio_E, const long double &norm_f);
    
    long double Evaluate_df_Chi2_d_sqr_N1(const long double &ratio_E, const long double &norm_f);
    
    long double Evaluate_Chi2(const long double &ratio_E, const long double &norm_f);
    
    long double Evaluate_dChi2_dratio_E(const long double &ratio_E, const long double &norm_f);
    
    long double Evaluate_dChi2_dN1(const long double &ratio_E, const long double &norm_f);
    
    void Evaluate_Chi2_derivatives(long double &E_dChi2_dE, long double &dChi2_df, const long double &r, const long double &N1);
    
    void Evaluate_Chi2_derivatives_Finite_Difference(long double &E_dChi2_dE, long double &dChi2_df, const long double &r, const long double &N1);
    
    long double Evaluate_Length_Scale(const long double &norm_f);
    
    long double Evaluate_diff_Length_Scale_dN1(const long double &norm_f);
    
    long double Evaluate_diff_Length_Scale_d_sqr_N1(const long double &norm_f);
    
    void Polynomial_Interpolation(M1_Data_1D_Cheby &Chi2_M1_1D_Uniform, ofstream &output_Opt_Coefficients);
    
    void Fit_Convergence_Check(const int &VAR_NUM);
    
    void Fit_Convergence_Check(M1_Data_1D_Cheby *Chi2_M1_1D_Uniform, const int &flag_Write_Output, ofstream &out_L_inf_Norm, const int &VAR_NUM);
    
    long double Recompute_I0_Mapping(const long double &ratio_E, const long double &L_N3_ijk_Actual, const long double &L_N3_ijk_Fit);
    
    long double Compute_drI01_drI02(const long double &ratio_E, const long double &L_N3_ijk_Actual, const long double &L_N3_ijk_Fit);
    
    //********************************************************************
    // Routines Needed for Least Squares Module for Optimization of L_Chi2
    //********************************************************************
    long double Evaluate_diff_Length_Scale_dN1_Least_Squares(const long double &norm_f);
    long double Evaluate_Length_Scale_Least_Squares(const long double &norm_f);
    
    void Vandermonde_Interpolation_1var_Least_Squares();
    void Vandermonde_Interpolation_2vars_Least_Squares();
    // long double Evaluate_dratio_E_dLength_Scale(const long double &ratio_E, const long double &Length_Scale);
    void Precompute_Max_Ent_Solution_Least_Squares(M1_Data_1D_Cheby *Chi2_M1_1D_Uniform);
    void Precompute_Max_Ent_Solution_Least_Squares_Uniform_E_Chebyshev_N1();
    void Precompute_Final_Max_Ent_Solution(const int &Maximum_Entropy_Solution_Regime = 0);
    long double Least_Squares_Optimization_L_Chi2(M1_1D_Data_Pointer &M1_1D_Data);
    
    void Nested_Least_Squares_Optimization_L_Chi2(M1_1D_Data_Pointer &M1_1D_Data);
    
    void Polynomial_Interpolation_Non_Gray_M1(M1_1D_Data_Pointer &M1_1D_Data, ofstream &output_Opt_Coefficients);
    
    void Write_Coefficients_M1_Closure_Interp(ofstream &output_Opt_Coefficients);
    
    void Compute_dI2ij_dIn_Finite_Difference(const long double &r_I0, const long double &N1, const long double &phi);
    
    void Compute_L_ONE_L_TWO_Errors_Least_Squares_L_Chi2(M1_Data_1D_Cheby *Chi2_M1_1D_Uniform, const int &flag_Write_Output, ofstream &out_L_inf_Norm, const int &VAR_NUM);
    
    void Compute_L_ONE_L_TWO_Errors_Nested_Least_Squares_L_Chi2(M1_Data_1D_Cheby *Chi2_M1_1D_Uniform, const int &flag_Write_Output, ofstream &out_L_inf_Norm, const int &VAR_NUM);
    
    void Write_Max_Ent_Data_Matlab();
    
    void Write_Output_Data_Fit_Chi2_Matlab();
    
//     long double Chebyshev_First_Kind_Basis(const long double &x, const int &Index);
    
//     long double Chebyshev_First_Kind_Derivatives_Basis(const long double &x, const int &Index);

//     long double Chebyshev_Second_Kind_Basis(const long double &x, const int &Index);
};

inline void M1_Data_1D_Cheby :: Setup_Coefficients_HL_LL(M1_Data_1D_Cheby &Chi2_M1_1D_HL, M1_Data_1D_Cheby &Chi2_M1_1D_LL) {
    for (int i = 0; i < N_Points_f; i++) {
        Coefficient_Matrix_Fit_Chi2_HL[i] = Chi2_M1_1D_HL.Coefficient_Matrix_Fit_Chi2[i];
        Coefficient_Matrix_Fit_Chi2_LL[i] = Chi2_M1_1D_LL.Coefficient_Matrix_Fit_Chi2[i];
    }   
}

inline void M1_Data_1D_Cheby :: Copy_to(M1_Data_1D_Cheby *New_M1_Data_1D_Cheby) {
//     New_M1_Data_1D_Cheby->path_out = &path_out;
//     New_M1_Data_1D_Cheby->prefix = &prefix;
//     New_M1_Data_1D_Cheby->extension = &extension;
    
//     New_M1_Data_1D_Cheby->in_out_Optim = in_out_Optim;
    
    int N_pts_total;
    New_M1_Data_1D_Cheby->rec_Ncoeffs = rec_Ncoeffs;
    New_M1_Data_1D_Cheby->rec_Chi2 = rec_Chi2;
    
    New_M1_Data_1D_Cheby->id_proc = id_proc;
    New_M1_Data_1D_Cheby->num_proc = num_proc;
    
    New_M1_Data_1D_Cheby->Node_Distribution_E = Node_Distribution_E;
    New_M1_Data_1D_Cheby->Node_Distribution_N1 = Node_Distribution_N1;
    
    New_M1_Data_1D_Cheby->L_inf_Norm = L_inf_Norm;
    New_M1_Data_1D_Cheby->L2_Norm_Chi2 = L2_Norm_Chi2;
    
    New_M1_Data_1D_Cheby->index_N1 = index_N1;
    
    New_M1_Data_1D_Cheby->N_Points_E = N_Points_E;
    New_M1_Data_1D_Cheby->N_Points_f = N_Points_f;
    New_M1_Data_1D_Cheby->N_pts_Mob_Scale = N_pts_Mob_Scale;
    
    New_M1_Data_1D_Cheby->N_pts_L_Chi2_LS = N_pts_L_Chi2_LS;
    
    New_M1_Data_1D_Cheby->L_Chi2_0 = L_Chi2_0;
    
    New_M1_Data_1D_Cheby->index_Ncoeffs = index_Ncoeffs;
    New_M1_Data_1D_Cheby->index_rec_Chi2 = index_rec_Chi2;
    
    for (int i = 0; i < N_Points_f; i++) {
        New_M1_Data_1D_Cheby->Coefficient_Matrix_Fit_Chi2_HL[i] = Coefficient_Matrix_Fit_Chi2_HL[i];
        New_M1_Data_1D_Cheby->Coefficient_Matrix_Fit_Chi2_LL[i] = Coefficient_Matrix_Fit_Chi2_LL[i];
    }
    
    for (int i = 0; i < N_pts_L_Chi2_LS; i++) {
        New_M1_Data_1D_Cheby->Coefficient_Matrix_Fit_Mob_Scale[i] = Coefficient_Matrix_Fit_Mob_Scale[i];
    }
    
    N_pts_total = N_Points_E*N_Points_f*N_pts_Mob_Scale;
    
    for (int i = 0; i < N_pts_total; i++) {
        New_M1_Data_1D_Cheby->flag_test_Implementation = flag_test_Implementation;
        New_M1_Data_1D_Cheby->ratio_I0_NON_GRAY[i] = ratio_I0_NON_GRAY[i];
        New_M1_Data_1D_Cheby->E_NON_GRAY[i] = E_NON_GRAY[i];
        New_M1_Data_1D_Cheby->N1_1_NON_GRAY[i] = N1_1_NON_GRAY[i];
        New_M1_Data_1D_Cheby->Chi2_NON_GRAY[i] = Chi2_NON_GRAY[i];
        New_M1_Data_1D_Cheby->dChi2_drI0[i] = dChi2_drI0[i];
        New_M1_Data_1D_Cheby->dChi2_dI0[i] = dChi2_dI0[i];
        New_M1_Data_1D_Cheby->dChi2_dN1[i] = dChi2_dN1[i];
        New_M1_Data_1D_Cheby->d2Chi2_dN12[i] = d2Chi2_dN12[i];
        New_M1_Data_1D_Cheby->d3Chi2_dN13[i] = d3Chi2_dN13[i];
        New_M1_Data_1D_Cheby->d2Chi2_drI0_dN1[i] = d2Chi2_drI0_dN1[i];
        New_M1_Data_1D_Cheby->d3Chi2_drI0_dN1_2[i] = d3Chi2_drI0_dN1_2[i];
        New_M1_Data_1D_Cheby->d2Chi2_dI0_dN1[i] = d2Chi2_dI0_dN1[i];
        New_M1_Data_1D_Cheby->d3Chi2_dI0_dN1_2[i] = d3Chi2_dI0_dN1_2[i];
        
        New_M1_Data_1D_Cheby->f_Chi2_NON_GRAY[i] = f_Chi2_NON_GRAY[i];
        New_M1_Data_1D_Cheby->error_fit[i] = error_fit[i];
    }
    
    N_pts_total = N_Points_E*N_Points_f;
    for (int i = 0; i < N_pts_total; i++) {
        New_M1_Data_1D_Cheby->Coefficient_Matrix_Fit_Chi2[i] = Coefficient_Matrix_Fit_Chi2[i];
    }
    for (int i = 0; i < N_pts_total*N_pts_total; i++) {
        New_M1_Data_1D_Cheby->Coefficients_Vander_Matrix_2vars[i] = Coefficients_Vander_Matrix_2vars[i];
    }
}

inline void M1_Data_1D_Cheby :: Write_Max_Ent_Data_Matlab() {
    char path_out[256];
    int index;
    strcpy(path_out, getenv(PATHVAR));
    
    strcat(path_out, "/M1_Model/Non_Gray_M1_Closure/Max_Ent_Data_for_Matlab.dat");
        
    ofstream in_out;
    in_out.open(path_out);
        
    if (!in_out) {
        cout << "Max_Ent_Data_for_Matlab.dat could not be accessed!" << endl;    
    }
    
    if (in_out.good()) {
        for (int id_Mob_Scale = 0 ; id_Mob_Scale < N_pts_Mob_Scale; id_Mob_Scale++) {
            for (int index_e = 0 ; index_e < N_Points_E; index_e++) {
                for (int i_f = 0; i_f < N_Points_f; i_f++) {
                    index = (id_Mob_Scale*N_Points_E + index_e)*N_Points_f + i_f;
                    
                    in_out << ratio_I0_NON_GRAY[index] << setw(18) << N1_1_NON_GRAY[index] << setw(18) << Chi2_NON_GRAY[index] << setw(18) << f_Chi2_NON_GRAY[index] << endl;  
                    // << setw(18) << error_fit[index] << endl;
                }
            }
        }
        in_out.close();
    }
}

inline void M1_Data_1D_Cheby :: Write_Output_Data_Fit_Chi2_Matlab() {
    char path_out[256];
    int index;
    int num_points;
    double r_I0, N1;
    double Chi2, f_Chi2;
    num_points = 100;
    
    strcpy(path_out, getenv(PATHVAR));
    
    strcat(path_out, "/M1_Model/Non_Gray_M1_Closure/Output_Data_Fit_Chi2_Matlab.dat");
        
    ofstream in_out;
    in_out.open(path_out);
        
    if (!in_out) {
        cout << "Output_Data_Fit_Chi2_Matlab.dat could not be accessed!" << endl;    
    }
    
    if (in_out.good()) {
        cout << "***************************************************** \n";
        cout << "******* Writing Output Data Fit Chi2 Matlab ********* \n";
        cout << "***************************************************** \n";
        for (int i = 0; i < num_points; i++) {
            for (int j = 0; j < num_points; j++) {
                r_I0 = zeros_shifted(i, num_points, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
                N1 = zeros_shifted(j + num_points - 1, 2*(num_points - 1) + 1, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
                
                Chi2 = Evaluate_Chi2(r_I0, N1);
                f_Chi2 = Evaluate_g_Chi2(r_I0, N1);
                
                in_out << r_I0 << setw(18) << N1 << setw(18) << Chi2 << setw(18) << f_Chi2 << endl;
            }
        }
        in_out.close();
    }
}

#endif // _M1_DATA_1D_CHEBY_H_INCLUDED
