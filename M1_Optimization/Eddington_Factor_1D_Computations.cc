#ifndef _M1_Model_1D_H_INCLUDED
#include "../M1_Optimization/M1_Optimization.h"
#endif // _M1_Model_1D_H_INCLUDED

void Set_Closure_Free_Streaming(record_Chi2 *rec_Chi2_local, const int &id_count, M1_State_Param &M1_State) {
    if (M1_State.Regime == HYPERBOLIC_LIMIT) {
        M1_State.MOM_VEC[0] = 0.0;
    } else if (M1_State.Regime == LOGARITHMIC_LIMIT) {
        M1_State.MOM_VEC[0] = 1.0e32;
    }
    
    rec_Chi2_local[id_count].I0 = M1_State.MOM_VEC[0];
    rec_Chi2_local[id_count].ratio_I0 = M1_State.ratio_I0;
    rec_Chi2_local[id_count].N1 = M1_State.N1;
    rec_Chi2_local[id_count].Chi2 = pow(M1_State.N1, 2);            
}

void Compute_Eddington_Factor(record_Chi2 *rec_Chi2_local, const int &index_count, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, M1_State_Param &M1_State, const long double *x, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, long double *Sk_final) {
    long double Chi2;
    Finite_Diff_Parameters Finite_Diff_Chi2;
    if (M1_State.Boundary_Point == FREE_STREAMING_LIMIT) {
        Set_Closure_Free_Streaming(rec_Chi2_local, index_count, M1_State);
    } else if (M1_State.flag_Taylor_Series_Expansion && M1_State.flag_finite_difference_N1 == false) {
        rec_Chi2_local[index_count].I0 = M1_State.MOM_VEC[0];
        rec_Chi2_local[index_count].ratio_I0 = M1_State.ratio_I0;
        rec_Chi2_local[index_count].N1 = M1_State.N1;
        // In this case we perform Taylor series expansion in the vicinity of the free-streaming 
        // limit to compute the Eddington factor
        Setup_Taylor_Series_Coefficients(Finite_Diff_Chi2, M1_State, Var_index, num_points, Mobius_Scale_Params, struc_map_N1, Sk_final);
    
        rec_Chi2_local[index_count].Chi2 = Finite_Diff_Chi2.Taylor_Series(M1_State.N1);
        rec_Chi2_local[index_count].dChi2_dN1 = Finite_Diff_Chi2.Taylor_Series_First_Derivative(M1_State.N1);
    } else {
        Eddington_Factor(Chi2, x, Sk_final, M1_State);
        
        if (M1_State.Regime == HYPERBOLIC_LIMIT) {
            M1_State.MOM_VEC[0] = 0.0;
        } else if (M1_State.Regime == LOGARITHMIC_LIMIT) {
            M1_State.MOM_VEC[0] = 1.0e32;
        }
        
        if (M1_State.id_proc == DISPLAY_ID) {
            if (M1_State.display) {
                cout << "Chi_2 = " << Chi2 << endl;
                cout << "\n" << endl;
            }
        }
        
        rec_Chi2_local[index_count].I0 = M1_State.MOM_VEC[0];
        rec_Chi2_local[index_count].ratio_I0 = M1_State.ratio_I0;
        rec_Chi2_local[index_count].N1 = M1_State.N1;
        rec_Chi2_local[index_count].Chi2 = Chi2;
    }
    
    if (rec_Chi2_local[index_count].Chi2 > 1.0 || rec_Chi2_local[index_count].Chi2 < pow(rec_Chi2_local[index_count].N1, 2)) {
        cout << "Non Realizable closure for M1" << endl;
        cout << "E = " << rec_Chi2_local[index_count].I0 << "   "  << "N1_1 = " << rec_Chi2_local[index_count].N1 << "   " << "Chi2 = " << rec_Chi2_local[index_count].Chi2 << "  " << "diff_lower = " << rec_Chi2_local[index_count].Chi2 - pow(rec_Chi2_local[index_count].N1, 2) << endl;
        exit(0);
    }
}

void Compute_Eddington_Factor_Derivatives(record_Chi2 *rec_Chi2_local, const int &index_count, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, M1_State_Param &M1_State, const long double *x, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, long double *Sk_final) {
    if (M1_State.Boundary_Point == FREE_STREAMING_LIMIT) {
        if ((num_points->Maximum_Entropy_Solution_Regime == HYPERBOLIC_LIMIT)  ||
            (num_points->Maximum_Entropy_Solution_Regime == LOGARITHMIC_LIMIT) ||
            (M1_State.Regime == BOSE_EINSTEIN) ||
            (M1_State.Regime == GRAY)) {
            
            if (!M1_State.flag_finite_difference_N1) {
                // We call the following routine only if there is no finite differencing
                // procedure already initiated
                Calculate_Eddington_Factor_Derivatives_Boundary(rec_Chi2_local[index_count], M1_State, Var_index, num_points, &Mobius_Scale_Params, struc_map_N1, Sk_final);
            }
        } else if ((M1_State.Regime == HYPERBOLIC_LIMIT) || 
                   (M1_State.Regime == LOGARITHMIC_LIMIT)   ) {
//             if (!M1_State.flag_finite_difference_r_I0) {
                // We call the following routine only if there is no finite differencing
                // procedure already initiated
                // Calculate_Eddington_Factor_Derivatives_Boundary_Spectrum_Energy(rec_Chi2_local[index_count], M1_State, Var_index, num_points, &Mobius_Scale_Params, struc_map_N1, Sk_final);
//             } else 
            if (!M1_State.flag_finite_difference_N1) {
                // We call the following routine only if finite differencing for r_I0 is already initiated for the free-
                // streaming limit, which means that finite differencing for N1 must now be initiated if not already
                Calculate_Eddington_Factor_Derivatives_Boundary(rec_Chi2_local[index_count], M1_State, Var_index, num_points, &Mobius_Scale_Params, struc_map_N1, Sk_final);
            }
        } else {
            cout << "Invalid Regime type !!!!!!!!!!!!!!!!!!!" << endl;
            exit(0);
        }
    } else {
        if ((num_points->Maximum_Entropy_Solution_Regime == HYPERBOLIC_LIMIT) || 
            (num_points->Maximum_Entropy_Solution_Regime == LOGARITHMIC_LIMIT) ||
            (M1_State.Regime == BOSE_EINSTEIN) ||
            (M1_State.Regime == GRAY)) {
            
            if (!M1_State.flag_Taylor_Series_Expansion) {
                // Derivatives with respect to N1 already computed in the case of taylor series
                // expansion
                Calculate_Eddington_Factor_Derivatives(rec_Chi2_local[index_count], x, Sk_final, M1_State);
            }
        } else if ((M1_State.Regime == HYPERBOLIC_LIMIT) || 
                   (M1_State.Regime == LOGARITHMIC_LIMIT)   ) {
            M1_State.Finite_Diff_Spectrum_Boundary_Type = M1_State.Regime;
//             if (!M1_State.flag_finite_difference_r_I0) {
                // We call the following routine only if there is no finite differencing
                // procedure already initiated for r_I0
                // Calculate_Eddington_Factor_Derivatives_Boundary_Spectrum_Energy(rec_Chi2_local[index_count], M1_State, Var_index, num_points, &Mobius_Scale_Params, struc_map_N1, Sk_final);
//             } else {
                // We call the following routine only if finite differencing for r_I0 is already initiated
                Calculate_Eddington_Factor_Derivatives(rec_Chi2_local[index_count], x, Sk_final, M1_State);
//             }
        } else {
            cout << "Invalid Regime type !!!!!!!!!!!!!!!!!!!" << endl;
            exit(0);
        }
    }
}

void Setup_Taylor_Series_Coefficients(Finite_Diff_Parameters &Finite_Diff_Chi2, M1_State_Param &M1_State, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, long double *Sk_final) {
    // Perform finite difference approximation to obtain derivative
    M1_State_Param M1_State_temp(M1_State);
    record_Chi2 rec_Chi2_array;
    record_Partial_Moments rec_PM_local;
    long double *Chi2_vals;
    long double h;
    int order, prec, n, nmax;
    long double *c, *x;
    long double diff_1;
    order = 4;
    prec = 4;
    n = order + prec;
    nmax = n;
    c = new long double[n];
    x = new long double[n];
    Chi2_vals = new long double[n];
    h = finite_diff_h;
    
    Finite_Diff_Chi2.x0 = 1.0;
    Finite_Diff_Chi2.Chi2_x0 = 1.0;
    
    M1_State_temp.f_test_knot = 1.0;
    M1_State_temp.h_f_test = h;
    
    M1_State_temp.flag_finite_difference_N1 = true;
    
    // backward finite difference approximation to fourth derivative
    differ_backward ( h, order, prec, c, x );
    
    for (int id_diff = 0; id_diff < n; id_diff++) {
        M1_State_temp.index_finite_diff = id_diff;
        M1_State_temp.x_finite_diff_N1 = x[id_diff];
        
        M1_State_temp.finite_difference_type_N1 = TAYLOR_SERIES_EXPANSION;
        
        NLOPT_Optimization_Algorithm_NON_GRAY_M1(&rec_Chi2_array, &rec_PM_local, 0, Var_index, num_points, M1_State.id_proc, M1_State_temp, Mobius_Scale_Params, struc_map_N1, Sk_final);
        Chi2_vals[id_diff] = rec_Chi2_array.Chi2;
    }
    
    diff_1 = 0.0;
    for (int i = 0; i < n; i++) {
        diff_1 += c[n-i-1] * Chi2_vals[nmax-i-1];
    }
    Finite_Diff_Chi2.d4Chi2_dN1 = diff_1;
    
    // backward finite difference approximation to third derivative
    order = 3;
    n = order + prec;
    differ_backward ( h, order, prec, c, x );
    diff_1 = 0.0;
    for (int i = 0; i < n; i++) {
        diff_1 += c[n-i-1] * Chi2_vals[nmax-i-1];
    }
    Finite_Diff_Chi2.d3Chi2_dN1 = diff_1;
    
    // backward finite difference approximation to second derivative
    order = 2;
    n = order + prec;
    differ_backward ( h, order, prec, c, x );
    diff_1 = 0.0;
    for (int i = 0; i < n; i++) {
        diff_1 += c[n-i-1] * Chi2_vals[nmax-i-1];
    }
    Finite_Diff_Chi2.d2Chi2_dN1 = diff_1;
    
    // backward finite difference approximation to first derivative
    order = 1;
    n = order + prec;
    differ_backward ( h, order, prec, c, x );
    diff_1 = 0.0;
    for (int i = 0; i < n; i++) {
        diff_1 += c[n-i-1] * Chi2_vals[nmax-i-1];
    }
    Finite_Diff_Chi2.dChi2_dN1 = diff_1;
    
//     Finite_Diff_Chi2.flag_finite_diff_setup = 1;
    
    if (M1_State_temp.display) {
        if (M1_State_temp.id_proc == DISPLAY_ID) {
            for (int id_diff = 0; id_diff < nmax; id_diff++) {
                cout << "id_diff = " << id_diff << "   " << "x = " << x[id_diff] << "   " << "Chi2 = " << Chi2_vals[id_diff] << endl;
            }
            cout << "dChi2_dN1 = " << Finite_Diff_Chi2.dChi2_dN1 << "  " << "d2Chi2_dN1 = " << Finite_Diff_Chi2.d2Chi2_dN1 << "  " << "d3Chi2_dN1 = " << Finite_Diff_Chi2.d3Chi2_dN1 << "  " << "d4Chi2_dN1 = " << Finite_Diff_Chi2.d4Chi2_dN1 << endl;
        }
    }
    
    delete[] c;
    delete[] x;
    delete[] Chi2_vals;
 }

void Calculate_Eddington_Factor_Derivatives_Boundary_Spectrum_Energy(record_Chi2 &rec_Chi2, M1_State_Param &M1_State, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, Mobius_Scale_Parameters *Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, long double *Sk_final) {
    // Perform finite difference approximation to obtain derivative
    M1_State_Param M1_State_temp(M1_State);
    record_Chi2 rec_Chi2_array;
    record_Partial_Moments rec_PM_local;
    long double *Chi2_vals, *d_Chi2_vals_dN1, *d2_Chi2_vals_dN1_2;
    long double h;
    int order, prec, n;
    long double *c, *x;
    long double diff_Chi2_val_dr, diff_2_Chi2_val_drI0_dN1, diff_3_Chi2_val_drI0_dN1_2;
    order = 1;
    prec = 4;
    n = order + prec;
    c = new long double[n];
    x = new long double[n];
    Chi2_vals = new long double[n];
    d_Chi2_vals_dN1 = new long double[n];
    d2_Chi2_vals_dN1_2 = new long double[n];
    h = 1.0*1.0e-3;
    
    // Set finite difference flag to true
    M1_State_temp.flag_finite_difference_r_I0 = true;
    
    switch (M1_State_temp.Finite_Diff_Spectrum_Boundary_Type) {
        case HYPERBOLIC_LIMIT:
            M1_State_temp.ratio_I0_knot = -1.0;
            // forward finite difference approximation to first derivative
            differ_forward ( h, order, prec, c, x );
            break;
        case LOGARITHMIC_LIMIT:
            M1_State_temp.ratio_I0_knot = 1.0;
            // backward finite difference approximation to first derivative
            differ_backward ( h, order, prec, c, x );
            break;
        default:
            cout << "Finite_Diff_Spectrum_Boundary_Type not specified !!!!" << endl;
            exit(0);
            break;
    }
    
    M1_State_temp.h_ratio_I0 = h;
    
    for (int id_diff = 0; id_diff < n; id_diff++) {
        M1_State_temp.index_finite_diff = id_diff;
        M1_State_temp.x_finite_diff_r_I0 = x[id_diff];
        
        M1_State_temp.finite_difference_type_r_I0 = DERIVATIVE_BOUNDARY_RADIATIVE_ENERGY;
        
        NLOPT_Optimization_Algorithm_NON_GRAY_M1(&rec_Chi2_array, &rec_PM_local, 0, Var_index, num_points, 0, M1_State_temp, *Mobius_Scale_Params, struc_map_N1, Sk_final);
        
        Chi2_vals[id_diff] = rec_Chi2_array.Chi2;
        d_Chi2_vals_dN1[id_diff] = rec_Chi2_array.dChi2_dN1;
        d2_Chi2_vals_dN1_2[id_diff] = rec_Chi2_array.d2Chi2_dN12;
    }
    
    diff_Chi2_val_dr = 0.0;
    diff_2_Chi2_val_drI0_dN1 = 0.0;
    diff_3_Chi2_val_drI0_dN1_2 = 0.0;
    
    for (int i = 0; i < n; i++) {
        diff_Chi2_val_dr += c[i] * Chi2_vals[i];
        diff_2_Chi2_val_drI0_dN1 += c[i] * d_Chi2_vals_dN1[i];
        diff_3_Chi2_val_drI0_dN1_2 += c[i] * d2_Chi2_vals_dN1_2[i];
    }
    
    rec_Chi2.dChi2_dN1 = d_Chi2_vals_dN1[0];
    rec_Chi2.d2Chi2_dN12 = d2_Chi2_vals_dN1_2[0];
    
    rec_Chi2.dChi2_drI0 = diff_Chi2_val_dr;
    rec_Chi2.d2_Chi2_drI0_dN1 = diff_2_Chi2_val_drI0_dN1;
    rec_Chi2.d3Chi2_drI0_dN1_2 = diff_3_Chi2_val_drI0_dN1_2;
    
//     long double r_N1, drI0_dI0, Length_Scale;
//     
//     if (M1_State_temp.Node_Dist_f == UNIFORM_ALGEBRAIC_MAPPING_N1 || M1_State_temp.Node_Dist_f == CHEBYSHEV_ALGEBRAIC_MAPPING_N1) {
//         r_N1 = Mapping_L_Chi2(struc_map_N1->f_L_Chi2, M1_State_temp.N1);
//         Length_Scale = Mobius_Scale_Params->Evaluate_Length_Scale_r_N1(r_N1);
//     } else if (struc_map_N1->flag_Algebraic_Map_N1) {
//         r_N1 = Mapping_L_Chi2(struc_map_N1->f_L_Chi2, M1_State_temp.N1);
//         Length_Scale = Mobius_Scale_Params->Evaluate_Length_Scale_r_N1(r_N1);
//     } else {
//         Length_Scale = Mobius_Scale_Params->Evaluate_Length_Scale(M1_State_temp.N1);
//     }
//     
//     drI0_dI0 = d_ratio_E_d_I0(M1_State_temp.MOM_VEC[0], Length_Scale);
//     
//     rec_Chi2.dChi2_dI0 = M1_State_temp.MOM_VEC[0]*rec_Chi2.dChi2_drI0*drI0_dI0 - M1_State_temp.N1*rec_Chi2.dChi2_dN1;
//     rec_Chi2.d2_Chi2_dI0_dN1 = M1_State_temp.MOM_VEC[0]*rec_Chi2.d2_Chi2_dI0_dN1*drI0_dI0 - M1_State_temp.N1*rec_Chi2.d2Chi2_dN12;
//     rec_Chi2.d3Chi2_dI0_dN1_2 = M1_State_temp.MOM_VEC[0]*rec_Chi2.d3Chi2_dI0_dN1_2*drI0_dI0 - M1_State_temp.N1*rec_Chi2.d3Chi2_dN13;
    
    // cout << "I0 = " << rec_Chi2.I0 << "   " << "N1 = " << rec_Chi2.N1 << "   " << "dChi2_drI0 = " << rec_Chi2.dChi2_drI0 << "   " << "d2_Chi2_drI0_dN1 = " << rec_Chi2.d2_Chi2_drI0_dN1 << "   " << "d3Chi2_drI0_dN1_2 = " << rec_Chi2.d3Chi2_drI0_dN1_2 << endl;
    
    if (M1_State_temp.display) {
        if (M1_State_temp.id_proc == DISPLAY_ID) {
            cout << "dChi2_drI0 = " << rec_Chi2.dChi2_drI0 << "   " << "I0 = " << M1_State_temp.I0 << "   " << "ratio_I0 = " << M1_State_temp.ratio_I0 << endl;
        }
    }
    
    delete[] c;
    delete[] x;
    delete[] Chi2_vals;
    delete[] d_Chi2_vals_dN1;
    delete[] d2_Chi2_vals_dN1_2;
}
 
void Calculate_Eddington_Factor_Derivatives_Boundary(record_Chi2 &rec_Chi2, M1_State_Param &M1_State, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, Mobius_Scale_Parameters *Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, long double *Sk_final) {
    // Perform finite difference approximation to obtain derivative
    M1_State_Param M1_State_temp(M1_State);
    record_Chi2 rec_Chi2_array;
    record_Partial_Moments rec_PM_Finite_Diff;
    long double *Chi2_vals;
    long double h;
    int order, prec, n, nmax;
    long double *c, *x;
    long double diff_chi2_dN1, diff_2_chi2_dN1_2, diff_3_chi2_dN1_3;
    order = 3;
    prec = 4;
    n = order + prec;
    nmax = n;
    c = new long double[n];
    x = new long double[n];
    Chi2_vals = new long double[n];
    h = finite_diff_h;
    
    M1_State_temp.f_test_knot = 1.0;
    M1_State_temp.h_f_test = h;
    
    // Set finite difference flag to true
    M1_State_temp.flag_finite_difference_N1 = true;
    
    // backward finite difference approximation to third derivative
    differ_backward ( h, order, prec, c, x );
    
    for (int id_diff = 0; id_diff < n; id_diff++) {
        M1_State_temp.index_finite_diff = id_diff;
        M1_State_temp.x_finite_diff_N1 = x[id_diff];
        
        M1_State_temp.finite_difference_type_N1 = DERIVATIVE_BOUNDARY;
        
        NLOPT_Optimization_Algorithm_NON_GRAY_M1(&rec_Chi2_array, &rec_PM_Finite_Diff, 0, Var_index, num_points, 0, M1_State_temp, *Mobius_Scale_Params, struc_map_N1, Sk_final);
        Chi2_vals[id_diff] = rec_Chi2_array.Chi2;
    }
    
    diff_3_chi2_dN1_3 = 0.0;
    for (int i = 0; i < n; i++) {
        diff_3_chi2_dN1_3 += c[n-i-1] * Chi2_vals[nmax-i-1];
    }
    rec_Chi2.d3Chi2_dN13 = diff_3_chi2_dN1_3;
    
    // backward finite difference approximation to second derivative
    order = 2;
    n = order + prec;
    differ_backward ( h, order, prec, c, x );
    diff_2_chi2_dN1_2 = 0.0;
    for (int i = 0; i < n; i++) {
        diff_2_chi2_dN1_2 += c[n-i-1] * Chi2_vals[nmax-i-1];
    }
    rec_Chi2.d2Chi2_dN12 = diff_2_chi2_dN1_2;
    
    // backward finite difference approximation to first derivative
    order = 1;
    n = order + prec;
    differ_backward ( h, order, prec, c, x );
    diff_chi2_dN1 = 0.0;
    for (int i = 0; i < n; i++) {
        diff_chi2_dN1 += c[n-i-1] * Chi2_vals[nmax-i-1];
    }
    rec_Chi2.dChi2_dN1 = diff_chi2_dN1;
    
    if (M1_State_temp.display) {
        if (M1_State_temp.id_proc == DISPLAY_ID) {
            for (int id_diff = 0; id_diff < n; id_diff++) {
                cout << "id_diff = " << id_diff << "   " << "x = " << x[n-id_diff-1] << "   " << "c = " << c[n-id_diff-1] << "   " << "Chi2 = " << Chi2_vals[nmax-id_diff-1] << endl;
            }
            cout << "dChi2_dN1 = " << rec_Chi2.dChi2_dN1 << "   " << "I0 = " << M1_State_temp.I0 << "   " << "ratio_I0 = " << M1_State_temp.ratio_I0 << endl;
        }
    }
    
//     switch (M1_State.Regime) {
//         case GRAY:
//             break;
//          case BOSE_EINSTEIN:
//              break;
//          case LOGARITHMIC_LIMIT:
//              rec_Chi2.dChi2_dN1 = 1.0;
//              break;
//          case HYPERBOLIC_LIMIT:
//              rec_Chi2.dChi2_dN1 = 2.0;
//              break;
//      }
    
    delete[] c;
    delete[] x;
    delete[] Chi2_vals;
 }
 
 void Calculate_Eddington_Factor_Derivatives(record_Chi2 &rec_Chi2, const long double *x, const long double *Sk, M1_State_Param &M1_State) {
     int NF_Chi = 2;
     long double *dlam_dIn, *dIn_dlam, *dChi2_dlam;
     dlam_dIn = new long double[M1_State.NVARS*M1_State.NVARS];
     dIn_dlam = new long double[M1_State.NVARS*M1_State.NVARS];
     dChi2_dlam = new long double[M1_State.NVARS];
     
     M1_State.set_x(x);
     M1_State.set_Sk(Sk);
     M1_State.set_Sk_cur(Sk);
     
     Compute_A_ONE_Matrix(dlam_dIn, x, Sk, M1_State);
     Compute_H_ONE_Matrix(dIn_dlam, x, Sk, M1_State);
     
     long double temp_val_E, temp_val_N1_1;
     long double diff_temp_val_E, diff_temp_val_N1_1;
     for (int i = 0; i < M1_State.NVARS; i++) {
         temp_val_E = 0.0;
         temp_val_N1_1 = 0.0;
         for (int j = 0; j < M1_State.NVARS; j++) {
             temp_val_E += dIn_dlam[0*M1_State.NVARS+j]*dlam_dIn[j*M1_State.NVARS+i];
             temp_val_N1_1 += dIn_dlam[1*M1_State.NVARS+j]*dlam_dIn[j*M1_State.NVARS+i];
        }
        diff_temp_val_E = fabs(temp_val_E - dIn_dInm1(0, i, M1_State.I0, M1_State.N1_1, M1_State.N1_2, M1_State.N1_3));
        diff_temp_val_N1_1 = fabs(temp_val_N1_1 - dIn_dInm1(1, i, M1_State.I0, M1_State.N1_1, M1_State.N1_2, M1_State.N1_3));
        
//         if (diff_temp_val_E > 1.0e-6 || diff_temp_val_N1_1 > 1.0e-6) {
//             cout << "temp_val_E = " << temp_val_E << "   " << "temp_val_N1_1 = " << temp_val_N1_1 << endl;
//             
//             cout << "diff_temp_val_E = " << diff_temp_val_E << "   " << "diff_temp_val_N1_1 = " << diff_temp_val_N1_1 << endl;
//             exit(0);
//         }
     }
     
     One_Dimensional_Quadrature(dChi2_dlam, NF_Chi, Angular_Moments_First_Derivatives, &M1_State);
     
     rec_Chi2.dChi2_dI0 = 0.0;
     rec_Chi2.dChi2_dN1 = 0.0;
     for (int i = 0; i < M1_State.NVARS; i++) {
         rec_Chi2.dChi2_dI0 += dChi2_dlam[i]*dlam_dIn[i*M1_State.NVARS + 0];
         rec_Chi2.dChi2_dN1 += dChi2_dlam[i]*dlam_dIn[i*M1_State.NVARS + 1];
     }
     rec_Chi2.dChi2_dI0 = rec_Chi2.dChi2_dI0; // - rec_Chi2.Chi2;
     
     rec_Chi2.d2Chi2_dN12 = Calculate_Angular_Moments_Second_Derivatives(dChi2_dlam, dlam_dIn, 1, 1, M1_State);
     
     // rec_Chi2.d3Chi2_dN13 = Calculate_Angular_Moments_Third_Derivatives(dChi2_dlam, dlam_dIn, 1, 1, 1, M1_State);
     
    if (M1_State.display) {
        if (M1_State.id_proc == DISPLAY_ID) {
            cout << "dChi2_dN1 = " << rec_Chi2.dChi2_dN1 << "   " << "d2Chi2_dN12 = " << rec_Chi2.d2Chi2_dN12 << endl;
            cout << "********************************************************************" << endl;
            cout << endl;
        }
     }
     
     delete[] dlam_dIn;
     delete[] dIn_dlam;
     delete[] dChi2_dlam;
 }
 
void Eddington_Factor(long double &Chi_2, const long double *x, const long double *Sk, M1_State_Param &M1_State) {
    int NF_Chi;
    long double *Moments_Vals;
    NF_Chi = 2;
    Moments_Vals = new long double[NF_Chi];
    
    M1_State.set_x(x);
    M1_State.set_Sk(Sk);
    
    One_Dimensional_Quadrature(Moments_Vals, NF_Chi, Chi2_M1_1D, &M1_State);
    
//     Moments_Vals[0] /= Sk[0];
//     Moments_Vals[1] = (Moments_Vals[1] - pow(Sk[2],2)*Moments_Vals[0] - 2.0*Sk[2]*Sk[3]*M1_State.N1*Moments_Vals[0])/pow(Sk[3], 2);
    
//     if (id == DISPLAY_ID) {
    if (M1_State.display) {
        if (M1_State.id_proc == DISPLAY_ID) {
            cout << "E = " << Moments_Vals[0] << "    " << "P = " << Moments_Vals[1] << endl;
        }
    }
//     }
    Chi_2 = Moments_Vals[1]/Moments_Vals[0];
    
    if (fabs(M1_State.MOM_VEC[0] - Moments_Vals[0])/M1_State.MOM_VEC[0] > tol_grad) {
        cout << "********************* Tolerance not satisfied in reality !!!!!! ************************" << fabs(M1_State.MOM_VEC[0] - Moments_Vals[0]) << endl;
        exit(0);
    }
    
    delete[] Moments_Vals;
 }

 int Chi2_M1_1D(long double *MOMENTS, const int &NFUN, const int &Id_angle, void *fdata) {
     long double coeff_1;
     long double sum_Lag_temp;
     M1_State_Param *M1_State = (M1_State_Param *) fdata;
     
     coeff_1 = 1.0;
     
     long double *poly, *poly_Sk;
     poly = new long double[M1_State->NVARS];
     
     switch (M1_State->Dimension) {
         case ONE_DIMENSIONAL:
             generate_polynomials_1D(poly, Id_angle, *M1_State);
             break;
         case THREE_DIMENSIONAL:
             generate_polynomials_3D(poly, Id_angle, *M1_State);
             break;
         default:
             cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
             exit(0);
             break;
     }
     
     poly_Sk = new long double[M1_State->NVARS];
     for (int i = 0; i < M1_State->NVARS; i++) {
         poly_Sk[i] = 0.0;
         for (int j = 0; j < M1_State->NVARS; j++) {
             poly_Sk[i] += M1_State->Sk[i*M1_State->NVARS+j]*poly[j];
        }
     }
     
     switch (M1_State->Regime) {
        case GRAY:
            MOMENTS[0] = poly_Sk[0]/(pow(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 4));
            MOMENTS[1] = pow(poly_Sk[1],2)/(pow(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 4));
            break;
         case BOSE_EINSTEIN:
             MOMENTS[0] = poly_Sk[0]*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS))/(1.0 - exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS)));
             MOMENTS[1] = pow(poly_Sk[1],2)*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS))/(1.0 - exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS)));
             break;
         case LOGARITHMIC_LIMIT:
             MOMENTS[0] = poly_Sk[0]/(-sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             MOMENTS[1] = pow(poly_Sk[1],2)/(-sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             break;
         case HYPERBOLIC_LIMIT:
             MOMENTS[0] = poly_Sk[0]*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             MOMENTS[1] = pow(poly_Sk[1],2)*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             break;
     }
     
     delete[] poly;
     delete[] poly_Sk;
     return 0;
 }

// ******************************************************************************
// This routine computes entries of the first-derivatives of the Eddington factor 
// with respect to the Lagrange multipliers
// ===> {\partial I^{(N)}} {\partial \lam_j}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
 int Angular_Moments_First_Derivatives(long double *MOMENTS, const int &NFUN, const int &Id_angle, void *fdata) {
     long double coeff_1;
     M1_State_Param *M1_State = (M1_State_Param *) fdata;
     coeff_1 = 1.0;
    
     long double *poly, *poly_Sk;
     long double *poly_Sk_cur, poly_ak;
     poly = new long double[M1_State->NVARS];
     poly_Sk = new long double[M1_State->NVARS];
     poly_Sk_cur = new long double[M1_State->NVARS];
     
     switch (M1_State->Dimension) {
         case ONE_DIMENSIONAL:
             generate_polynomials_1D(poly, Id_angle, *M1_State);
             
             switch (M1_State->Max_Ent_Solution_Type) {
                 case CLOSING_FLUX:
                     // In this case we compute the Eddington factor in the 1D case
                     poly_ak = pow(M1_State->x_quad[Id_angle], 2);
                     break;
                 case PARTIAL_MOMENTS:
                     cout << " ********************* Partial Moments must be computed in the 3D case !!!!! ********************* " << endl;
                     exit(0);
                     break;
                 default:
                     cout << " ********************* Max_Ent_Solution_Type Not Specified !!!!! ********************* " << endl;
                     exit(0);
                     break;
            }
            break;
         case THREE_DIMENSIONAL:
             generate_polynomials_3D(poly, Id_angle, *M1_State);
             
             switch (M1_State->Max_Ent_Solution_Type) {
                 case CLOSING_FLUX:
                     cout << " ********************* Eddington factor must be computed in the 1D case !!!!! ********************* " << endl;
                     break;
                 case PARTIAL_MOMENTS:
                     // In this case we compute the partial moments in the 3D case
                     poly_ak = generate_partial_moments_monomials_basis(Id_angle, *M1_State);
                     break;
                 default:
                     cout << " ********************* Max_Ent_Solution_Type Not Specified !!!!! ********************* " << endl;
                     exit(0);
                     break;
            }
            break;
         default:
             cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
             exit(0);
             break;
     }
     
     for (int i = 0; i < M1_State->NVARS; i++) {
         poly_Sk[i] = 0.0;
         for (int j = 0; j < M1_State->NVARS; j++) {
             poly_Sk[i] += M1_State->Sk[i*M1_State->NVARS+j]*poly[j];
        }
     }
     
     for (int i = 0; i < M1_State->NVARS; i++) {
         poly_Sk_cur[i] = 0.0;
         for (int j = 0; j < M1_State->NVARS; j++) {
             poly_Sk_cur[i] += M1_State->Sk_cur[i*M1_State->NVARS+j]*poly[j];
         }
     }
     
     switch (M1_State->Regime) {
         case GRAY:
             for (int i = 0; i < M1_State->NVARS; i++) {
                 MOMENTS[i] = -4.0*poly_Sk_cur[i]*poly_ak/(pow(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 5));
             }
             break;
         case BOSE_EINSTEIN:
             for (int i = 0; i < M1_State->NVARS; i++) {
                 MOMENTS[i] = coeff_1*poly_Sk_cur[i]*poly_ak*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS))/pow((1.0 - exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS))), 2.0);
             }
             break;
         case HYPERBOLIC_LIMIT:
             for (int i = 0; i < M1_State->NVARS; i++) {
                 MOMENTS[i] = coeff_1*poly_Sk_cur[i]*poly_ak*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             }
             break;
         case LOGARITHMIC_LIMIT:
             for (int i = 0; i < M1_State->NVARS; i++) {
                 MOMENTS[i] = coeff_1*poly_Sk_cur[i]*poly_ak/pow(-sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 2.0);
            }
            break;
     }
     
     delete[] poly;
     delete[] poly_Sk_cur;
     delete[] poly_Sk;
     return 0;
 }
 
// ******************************************************************************
// This routine computes entries of the second-derivatives of the Eddington 
// factor with respect to the Lagrange multipliers
// ===> {\partial^{2} I^{(N)}} {\partial \lam_j \partial \lam_k}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
 int Angular_Moments_Second_Derivatives(long double *MOMENTS, const int &NFUN, const int &Id_angle, void *fdata) {
     long double coeff_1;
     coeff_1 = 1.0;
     long double exp_x;
     int index_Lag1, index_Lag2, index_Moment;
     long double *poly, *poly_Sk;
     long double poly_ak, poly_Lag_1, poly_Lag_2;
    
     M1_State_Param *M1_State = (M1_State_Param *) fdata;
     
     poly = new long double[M1_State->NVARS];
     poly_Sk = new long double[M1_State->NVARS];
     
     index_Moment = M1_State->index_Moment;
     index_Lag1 = M1_State->index_Lag_i;
     index_Lag2 = M1_State->index_Lag_j;
     
     switch (M1_State->Dimension) {
         case ONE_DIMENSIONAL:
             generate_polynomials_1D(poly, Id_angle, *M1_State);
             
             switch (M1_State->Max_Ent_Solution_Type) {
                 case CLOSING_FLUX:
                     // In this case we compute the Eddington factor in the 1D case
                     poly_ak = pow(M1_State->x_quad[Id_angle], 2);
                     break;
                 case PARTIAL_MOMENTS:
                     cout << " ********************* Partial Moments must be computed in the 3D case !!!!! ********************* " << endl;
                     exit(0);
                     break;
                 default:
                     cout << " ********************* Max_Ent_Solution_Type Not Specified !!!!! ********************* " << endl;
                     exit(0);
                     break;
            }
            break;
         case THREE_DIMENSIONAL:
             generate_polynomials_3D(poly, Id_angle, *M1_State);
             
             switch (M1_State->Max_Ent_Solution_Type) {
                 case CLOSING_FLUX:
                     cout << " ********************* Eddington factor must be computed in the 1D case !!!!! ********************* " << endl;
                     break;
                 case PARTIAL_MOMENTS:
                     // In this case we compute the partial moments in the 3D case
                     poly_ak = generate_partial_moments_monomials_basis(Id_angle, *M1_State);
                     break;
                 default:
                     cout << " ********************* Max_Ent_Solution_Type Not Specified !!!!! ********************* " << endl;
                     exit(0);
                     break;
            }
            break;
         default:
             cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
             exit(0);
             break;
     }
     
     poly_Sk = new long double[M1_State->NVARS];
     for (int i = 0; i < M1_State->NVARS; i++) {
         poly_Sk[i] = 0.0;
         for (int j = 0; j < M1_State->NVARS; j++) {
             poly_Sk[i] += M1_State->Sk[i*M1_State->NVARS+j]*poly[j];
        }
     }
     
     poly_Lag_1 = 0.0;
     poly_Lag_2 = 0.0;
     for (int i = 0; i < M1_State->NVARS; i++) {
         poly_Lag_1 += M1_State->ak[index_Lag1*M1_State->NVARS+i]*poly[i];
         poly_Lag_2 += M1_State->ak[index_Lag2*M1_State->NVARS+i]*poly[i];
     }
     
     switch (M1_State->Regime) {
         case GRAY:
             MOMENTS[0] = 20.0*poly_Lag_1*poly_Lag_2*poly_ak/(pow(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 6));
             break;
         case BOSE_EINSTEIN:
             exp_x = exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             MOMENTS[0] = coeff_1*poly_Lag_1*poly_Lag_2*poly_ak*exp_x*(1.0 + exp_x)/pow((1.0 - exp_x), 3.0);
             break;
         case HYPERBOLIC_LIMIT:
             MOMENTS[0] = coeff_1*poly_Lag_1*poly_Lag_2*poly_ak*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             break;
         case LOGARITHMIC_LIMIT:
             MOMENTS[0] = 2.0*coeff_1*poly_Lag_1*poly_Lag_2*poly_ak/pow(-sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 3.0);
            break;
     }
     
     delete[] poly;
     delete[] poly_Sk;
     return 0;
 }
 
// ******************************************************************************
// This routine computes entries of the third-derivatives of the Eddington factor
// with respect to the Lagrange multipliers
// ===> {\partial^{3} I^{(N)}} {\partial \lam_j \partial \lam_k \partial \lam_l}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
int Angular_Moments_Third_Derivatives(long double *MOMENTS, const int &NFUN, const int &Id_angle, void *fdata) {
     long double coeff_1, exp_x;
     int index_Lag1, index_Lag2, index_Lag3;
     M1_State_Param *M1_State = (M1_State_Param *) fdata;
     index_Lag1 = M1_State->index_Lag_i;
     index_Lag2 = M1_State->index_Lag_j;
     index_Lag3 = M1_State->index_Lag_k;
     coeff_1 = 1.0;
    
     long double *poly, *poly_Sk;
     long double poly_ak, poly_Lag_1, poly_Lag_2, poly_Lag_3;
     poly = new long double[M1_State->NVARS];
     
     switch (M1_State->Dimension) {
         case ONE_DIMENSIONAL:
             generate_polynomials_1D(poly, Id_angle, *M1_State);
             
             switch (M1_State->Max_Ent_Solution_Type) {
                 case CLOSING_FLUX:
                     // In this case we compute the Eddington factor in the 1D case
                     poly_ak = pow(M1_State->x_quad[Id_angle], 2);
                     break;
                 case PARTIAL_MOMENTS:
                     cout << " ********************* Partial Moments must be computed in the 3D case !!!!! ********************* " << endl;
                     exit(0);
                     break;
                 default:
                     cout << " ********************* Max_Ent_Solution_Type Not Specified !!!!! ********************* " << endl;
                     exit(0);
                     break;
            }
            break;
         case THREE_DIMENSIONAL:
             generate_polynomials_3D(poly, Id_angle, *M1_State);
             
             switch (M1_State->Max_Ent_Solution_Type) {
                 case CLOSING_FLUX:
                     cout << " ********************* Eddington factor must be computed in the 1D case !!!!! ********************* " << endl;
                     break;
                 case PARTIAL_MOMENTS:
                     // In this case we compute the partial moments in the 3D case
                     poly_ak = generate_partial_moments_monomials_basis(Id_angle, *M1_State);
                     break;
                 default:
                     cout << " ********************* Max_Ent_Solution_Type Not Specified !!!!! ********************* " << endl;
                     exit(0);
                     break;
            }
            break;
         default:
             cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
             exit(0);
             break;
     }
     
     poly_Sk = new long double[M1_State->NVARS];
     for (int i = 0; i < M1_State->NVARS; i++) {
         poly_Sk[i] = 0.0;
         for (int j = 0; j < M1_State->NVARS; j++) {
             poly_Sk[i] += M1_State->Sk[i*M1_State->NVARS+j]*poly[j];
        }
     }
     
     poly_Lag_1 = 0.0;
     poly_Lag_2 = 0.0;
     poly_Lag_3 = 0.0;
     for (int i = 0; i < M1_State->NVARS; i++) {
         poly_Lag_1 += M1_State->ak[index_Lag1*M1_State->NVARS+i]*poly[i];
         poly_Lag_2 += M1_State->ak[index_Lag2*M1_State->NVARS+i]*poly[i];
         poly_Lag_3 += M1_State->ak[index_Lag3*M1_State->NVARS+i]*poly[i];
     }
     
     switch (M1_State->Regime) {
         case GRAY:
             MOMENTS[0] = -120.0*poly_Lag_1*poly_Lag_2*poly_Lag_3*poly_ak/(pow(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 7));
             break;
         case BOSE_EINSTEIN:
             exp_x = exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             MOMENTS[0] = coeff_1*poly_Lag_1*poly_Lag_2*poly_Lag_3*poly_ak*exp_x*(1.0 + 4.0*exp_x + exp_x*exp_x)/pow((1.0 - exp_x), 4.0);
             break;
         case HYPERBOLIC_LIMIT:
             MOMENTS[0] = coeff_1*poly_Lag_1*poly_Lag_2*poly_Lag_3*poly_ak*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             break;
         case LOGARITHMIC_LIMIT:
             MOMENTS[0] = 6.0*coeff_1*poly_Lag_1*poly_Lag_2*poly_Lag_3*poly_ak/pow(-sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 4.0);
            break;
     }
     
     delete[] poly;
     delete[] poly_Sk;
     return 0;
 }
 
