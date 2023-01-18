#ifndef _M1_Model_1D_H_INCLUDED
#include "../M1_Optimization/M1_Optimization.h"
#endif // _M1_Model_1D_H_INCLUDED

void Set_Partial_Moments_Free_Streaming(record_Partial_Moments *rec_PM_local, const int &id_count, M1_State_Param &M1_State) {
    if (M1_State.Regime == HYPERBOLIC_LIMIT) {
        M1_State.MOM_VEC[0] = 0.0;
    } else if (M1_State.Regime == LOGARITHMIC_LIMIT) {
        M1_State.MOM_VEC[0] = 1.0e32;
    }
    
    rec_PM_local[id_count].I0 = M1_State.MOM_VEC[0];
    rec_PM_local[id_count].ratio_I0 = M1_State.ratio_I0;
    rec_PM_local[id_count].N1 = M1_State.N1;
    rec_PM_local[id_count].N1_1 = M1_State.N1_1;
    rec_PM_local[id_count].N1_2 = M1_State.N1_2;
    rec_PM_local[id_count].N1_3 = M1_State.N1_3;
    rec_PM_local[id_count].x_val = M1_State.x_val;
    rec_PM_local[id_count].y_val = M1_State.y_val;
    rec_PM_local[id_count].z_val = M1_State.z_val;
    
    rec_PM_local[id_count].I0_Plus = Heaviside(M1_State.N1_1);
    rec_PM_local[id_count].N1_1_Plus = M1_State.N1_1*Heaviside(M1_State.N1_1);
    rec_PM_local[id_count].N1_2_Plus = M1_State.N1_2*Heaviside(M1_State.N1_1);
    rec_PM_local[id_count].N2_11_Plus = pow(M1_State.N1_1,2)*Heaviside(M1_State.N1_1);
    rec_PM_local[id_count].N2_12_Plus = M1_State.N1_1*M1_State.N1_2*Heaviside(M1_State.N1_1);
    rec_PM_local[id_count].N2_22_Plus = pow(M1_State.N1_2,2)*Heaviside(M1_State.N1_1);
}  
    
void Compute_Partial_Moments(record_Partial_Moments *rec_PM_local, const int &index_count, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, M1_State_Param &M1_State, const long double *x, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, const long double *x_SH, const long double *y_SH, const long double *z_SH, long double *Sk_final) {
    long double *Partial_Moms;
    Partial_Moms = new long double[M1_State.NFUNS_PM_3D];
    
    Finite_Diff_Parameters Finite_Diff_PM;
    
    if (M1_State.Boundary_Point == FREE_STREAMING_LIMIT) {
        Set_Partial_Moments_Free_Streaming(rec_PM_local, index_count, M1_State);
    } else if (M1_State.flag_Taylor_Series_Expansion && M1_State.flag_finite_difference_N1 == false) {
        // In this case we perform Taylor series expansion in the vicinity of the free-streaming 
        // limit to compute the partial moments
        rec_PM_local[index_count].I0 = M1_State.MOM_VEC[0];
        rec_PM_local[index_count].ratio_I0 = M1_State.ratio_I0;
        rec_PM_local[index_count].N1 = M1_State.N1;
        rec_PM_local[index_count].N1_1 = M1_State.N1_1;
        rec_PM_local[index_count].N1_2 = M1_State.N1_2;
        rec_PM_local[index_count].N1_3 = M1_State.N1_3;
        
        rec_PM_local[index_count].x_val = M1_State.x_val;
        rec_PM_local[index_count].y_val = M1_State.y_val;
        rec_PM_local[index_count].z_val = M1_State.z_val;
        
        // In this case we perform Taylor series expansion in the vicinity of the free-streaming 
        // limit to compute the Eddington factor
        Setup_Taylor_Series_Coefficients_PM(Finite_Diff_PM, M1_State, Var_index, num_points, Mobius_Scale_Params, struc_map_N1, x_SH, y_SH, z_SH, Sk_final);
        
        Finite_Diff_PM.Taylor_Series_PM(rec_PM_local[index_count].I0_Plus, rec_PM_local[index_count].N1_1_Plus, rec_PM_local[index_count].N1_2_Plus, rec_PM_local[index_count].N2_11_Plus, rec_PM_local[index_count].N2_12_Plus, rec_PM_local[index_count].N2_22_Plus, rec_PM_local[index_count].N1);
        
        // rec_PM_local[index_count].dChi2_dN1 = Finite_Diff_PM.Taylor_Series_First_Derivative(M1_State.N1);
    } else {
        rec_PM_local[index_count].I0 = M1_State.MOM_VEC[0];
        rec_PM_local[index_count].ratio_I0 = M1_State.ratio_I0;
        rec_PM_local[index_count].N1 = M1_State.N1;
        rec_PM_local[index_count].N1_1 = M1_State.N1_1;
        rec_PM_local[index_count].N1_2 = M1_State.N1_2;
        rec_PM_local[index_count].N1_3 = M1_State.N1_3;
        
        rec_PM_local[index_count].x_val = M1_State.x_val;
        rec_PM_local[index_count].y_val = M1_State.y_val;
        rec_PM_local[index_count].z_val = M1_State.z_val;
        
        Partial_Moments_M1_NON_GRAY_3D(Partial_Moms, M1_State.NFUNS_PM_3D, x, Sk_final, M1_State);
        
        if (M1_State.Regime == HYPERBOLIC_LIMIT) {
            M1_State.MOM_VEC[0] = 0.0;
        } else if (M1_State.Regime == LOGARITHMIC_LIMIT) {
            M1_State.MOM_VEC[0] = 1.0e32;
        }
        
        rec_PM_local[index_count].I0_Plus = Partial_Moms[0];
        rec_PM_local[index_count].N1_1_Plus = Partial_Moms[1];
        rec_PM_local[index_count].N1_2_Plus = Partial_Moms[2];
        rec_PM_local[index_count].N2_11_Plus = Partial_Moms[3];
        rec_PM_local[index_count].N2_12_Plus = Partial_Moms[4];
        rec_PM_local[index_count].N2_22_Plus = Partial_Moms[5];
        
        if (M1_State.id_proc == DISPLAY_ID) {
            if (M1_State.display) {
                cout << "E_Plus = " << Partial_Moms[0] << "    " << "N1_1_Plus = " << Partial_Moms[1] << "    " << "N1_2_Plus = " << Partial_Moms[2] << "    " << "N2_11_Plus = " << Partial_Moms[3] << "    " << "N2_12_Plus = " << Partial_Moms[4] << "    " << "N2_22_Plus = " << Partial_Moms[5] << endl;
                cout << "\n" << endl;
            }
        }
    }
    
    delete[] Partial_Moms;
}

void Compute_Partial_Moments_Derivatives(record_Partial_Moments *rec_PM_local, const int &index_count, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, M1_State_Param &M1_State, const long double *x, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, const long double *x_SH, const long double *y_SH, const long double *z_SH, long double *Sk_final) {
    if (M1_State.Boundary_Point == FREE_STREAMING_LIMIT) {
        if ((num_points->Maximum_Entropy_Solution_Regime == HYPERBOLIC_LIMIT)  ||
            (num_points->Maximum_Entropy_Solution_Regime == LOGARITHMIC_LIMIT) ||
            (M1_State.Regime == BOSE_EINSTEIN) ||
            (M1_State.Regime == GRAY)) {
            
            if (!M1_State.flag_finite_difference_N1) {
                // We call the following routine only if there is no finite differencing
                // procedure already initiated
                Calculate_Partial_Moments_Derivatives_Boundary(rec_PM_local[index_count], M1_State, Var_index, num_points, Mobius_Scale_Params, struc_map_N1, x_SH, y_SH, z_SH, Sk_final);
            }
        } else if ((M1_State.Regime == HYPERBOLIC_LIMIT) || 
                   (M1_State.Regime == LOGARITHMIC_LIMIT)   ) {
//             if (!M1_State.flag_finite_difference_r_I0) {
                // We call the following routine only if there is no finite differencing
                // procedure already initiated
//                 Calculate_Partial_Moments_Derivatives_Boundary_Spectrum_Energy(rec_Chi2_local[index_count], M1_State, Var_index, num_points, &Mobius_Scale_Params, struc_map_N1, Sk_final);
//             } else 
                if (!M1_State.flag_finite_difference_N1) {
                // We call the following routine only if finite differencing for r_I0 is already initiated for the free-
                // streaming limit, which means that finite differencing for N1 must now be initiated if not already
                Calculate_Partial_Moments_Derivatives_Boundary(rec_PM_local[index_count], M1_State, Var_index, num_points, Mobius_Scale_Params, struc_map_N1, x_SH, y_SH, z_SH, Sk_final);
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
                Calculate_Partial_Moments_Derivatives(rec_PM_local[index_count], x, Sk_final, M1_State);
            }
        } else if ((M1_State.Regime == HYPERBOLIC_LIMIT) || 
                   (M1_State.Regime == LOGARITHMIC_LIMIT)   ) {
            M1_State.Finite_Diff_Spectrum_Boundary_Type = M1_State.Regime;
//             if (!M1_State.flag_finite_difference_r_I0) {
                // We call the following routine only if there is no finite differencing
                // procedure already initiated for r_I0
//                 Calculate_Partial_Moments_Derivatives_Boundary_Spectrum_Energy(rec_Chi2_local[index_count], M1_State, Var_index, num_points, &Mobius_Scale_Params, struc_map_N1, Sk_final);
//             } else {
                // We call the following routine only if finite differencing for r_I0 is already initiated
                Calculate_Partial_Moments_Derivatives(rec_PM_local[index_count], x, Sk_final, M1_State);
//             }
        } else {
            cout << "Invalid Regime type !!!!!!!!!!!!!!!!!!!" << endl;
            exit(0);
        }
    }
}

void Setup_Taylor_Series_Coefficients_PM(Finite_Diff_Parameters &Finite_Diff_PM, M1_State_Param &M1_State, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, const long double *x_SH, const long double *y_SH, const long double *z_SH, long double *Sk_final) {
    // Perform finite difference approximation to obtain derivative
    M1_State_Param M1_State_temp(M1_State);
    record_Partial_Moments rec_Chi2_Finite_Diff;
    record_Partial_Moments rec_PM_Finite_Diff;
    long double *E_Plus_vals, *N1_1_Plus_vals, *N1_2_Plus_vals, *N2_11_Plus_vals, *N2_12_Plus_vals, *N2_22_Plus_vals;
    long double h;
    int order, prec, n, nmax;
    long double *c, *x;
    long double diff_I0_Plus, diff_N1_1_Plus, diff_N1_2_Plus, diff_N2_11_Plus, diff_N2_12_Plus, diff_N2_22_Plus;
    order = 4;
    prec = 6;
    n = order + prec;
    nmax = n;
    c = new long double[n];
    x = new long double[n];
    E_Plus_vals = new long double[n];
    N1_1_Plus_vals = new long double[n];
    N1_2_Plus_vals = new long double[n];
    N2_11_Plus_vals = new long double[n];
    N2_12_Plus_vals = new long double[n];
    N2_22_Plus_vals = new long double[n];
    h = finite_diff_h;
    
    long double x_proj, y_proj, z_proj;
    
    x_proj = M1_State.N1_1/M1_State.N1;
    y_proj = M1_State.N1_2/M1_State.N1;
    z_proj = M1_State.N1_3/M1_State.N1;
    
    Finite_Diff_PM.x0 = 1.0;
    Finite_Diff_PM.I0_plus_x0 = Heaviside(M1_State.N1_1);
    Finite_Diff_PM.N1_1_plus_x0 = x_proj*Heaviside(M1_State.N1_1);
    Finite_Diff_PM.N1_2_plus_x0 = y_proj*Heaviside(M1_State.N1_1);
    Finite_Diff_PM.N2_11_plus_x0 = pow(x_proj,2)*Heaviside(M1_State.N1_1);
    Finite_Diff_PM.N2_12_plus_x0 = x_proj*y_proj*Heaviside(M1_State.N1_1);
    Finite_Diff_PM.N2_22_plus_x0 = pow(y_proj,2)*Heaviside(M1_State.N1_1);
    
    M1_State_temp.f_test_knot = 1.0;
    M1_State_temp.h_f_test = h;
    
    M1_State_temp.flag_finite_difference_N1 = true;
    
    // backward finite difference approximation to fourth derivative
    differ_backward ( h, order, prec, c, x );
    
    for (int id_diff = 0; id_diff < n; id_diff++) {
        M1_State_temp.index_finite_diff = id_diff;
        M1_State_temp.x_finite_diff_N1 = x[id_diff];
        
        M1_State_temp.finite_difference_type_N1 = TAYLOR_SERIES_EXPANSION;
        
        NLOPT_Optimization_Algorithm_NON_GRAY_M1_3D(NULL, &rec_PM_Finite_Diff, 0, Var_index, num_points, M1_State_temp.id_proc, M1_State_temp, Mobius_Scale_Params, struc_map_N1, Sk_final, x_SH, y_SH, z_SH);
        
        E_Plus_vals[id_diff] = rec_PM_Finite_Diff.I0_Plus;
        N1_1_Plus_vals[id_diff] = rec_PM_Finite_Diff.N1_1_Plus;
        N1_2_Plus_vals[id_diff] = rec_PM_Finite_Diff.N1_2_Plus;
        N2_11_Plus_vals[id_diff] = rec_PM_Finite_Diff.N2_11_Plus;
        N2_12_Plus_vals[id_diff] = rec_PM_Finite_Diff.N2_12_Plus;
        N2_22_Plus_vals[id_diff] = rec_PM_Finite_Diff.N2_22_Plus;
    }
    
    diff_I0_Plus = 0.0;
    diff_N1_1_Plus = 0.0;
    diff_N1_2_Plus = 0.0;
    diff_N2_11_Plus = 0.0;
    diff_N2_12_Plus = 0.0;
    diff_N2_22_Plus = 0.0;
    for (int i = 0; i < n; i++) {
        diff_I0_Plus += c[n - 1 - i] * E_Plus_vals[nmax - i - 1];
        diff_N1_1_Plus += c[n - 1 - i] * N1_1_Plus_vals[nmax - i - 1];
        diff_N1_2_Plus += c[n - 1 - i] * N1_2_Plus_vals[nmax - i - 1];
        diff_N2_11_Plus += c[n - 1 - i] * N2_11_Plus_vals[nmax - i - 1];
        diff_N2_12_Plus += c[n - 1 - i] * N2_12_Plus_vals[nmax - i - 1];
        diff_N2_22_Plus += c[n - 1 - i] * N2_22_Plus_vals[nmax - i - 1];
    }
    Finite_Diff_PM.d4_I0_Plus_dnorm_f_4 = diff_I0_Plus;
    Finite_Diff_PM.d4_N1_1_Plus_dnorm_f_4 = diff_N1_1_Plus;
    Finite_Diff_PM.d4_N1_2_Plus_dnorm_f_4 = diff_N1_2_Plus;
    Finite_Diff_PM.d4_N2_11_Plus_dnorm_f_4 = diff_N2_11_Plus;
    Finite_Diff_PM.d4_N2_12_Plus_dnorm_f_4 = diff_N2_12_Plus;
    Finite_Diff_PM.d4_N2_22_Plus_dnorm_f_4 = diff_N2_22_Plus;
    
    // backward finite difference approximation to third derivative
    order = 3;
    n = order + prec;
    differ_backward ( h, order, prec, c, x );
    diff_I0_Plus = 0.0;
    diff_N1_1_Plus = 0.0;
    diff_N1_2_Plus = 0.0;
    diff_N2_11_Plus = 0.0;
    diff_N2_12_Plus = 0.0;
    diff_N2_22_Plus = 0.0;
    for (int i = 0; i < n; i++) {
        diff_I0_Plus += c[n - 1 - i] * E_Plus_vals[nmax - i - 1];
        diff_N1_1_Plus += c[n - 1 - i] * N1_1_Plus_vals[nmax - i - 1];
        diff_N1_2_Plus += c[n - 1 - i] * N1_2_Plus_vals[nmax - i - 1];
        diff_N2_11_Plus += c[n - 1 - i] * N2_11_Plus_vals[nmax - i - 1];
        diff_N2_12_Plus += c[n - 1 - i] * N2_12_Plus_vals[nmax - i - 1];
        diff_N2_22_Plus += c[n - 1 - i] * N2_22_Plus_vals[nmax - i - 1];
    }
    Finite_Diff_PM.d3_I0_Plus_dnorm_f_3 = diff_I0_Plus;
    Finite_Diff_PM.d3_N1_1_Plus_dnorm_f_3 = diff_N1_1_Plus;
    Finite_Diff_PM.d3_N1_2_Plus_dnorm_f_3 = diff_N1_2_Plus;
    Finite_Diff_PM.d3_N2_11_Plus_dnorm_f_3 = diff_N2_11_Plus;
    Finite_Diff_PM.d3_N2_12_Plus_dnorm_f_3 = diff_N2_12_Plus;
    Finite_Diff_PM.d3_N2_22_Plus_dnorm_f_3 = diff_N2_22_Plus;
    
    // backward finite difference approximation to third derivative
    order = 2;
    n = order + prec;
    differ_backward ( h, order, prec, c, x );
    diff_I0_Plus = 0.0;
    diff_N1_1_Plus = 0.0;
    diff_N1_2_Plus = 0.0;
    diff_N2_11_Plus = 0.0;
    diff_N2_12_Plus = 0.0;
    diff_N2_22_Plus = 0.0;
    for (int i = 0; i < n; i++) {
        diff_I0_Plus += c[n - 1 - i] * E_Plus_vals[nmax - i - 1];
        diff_N1_1_Plus += c[n - 1 - i] * N1_1_Plus_vals[nmax - i - 1];
        diff_N1_2_Plus += c[n - 1 - i] * N1_2_Plus_vals[nmax - i - 1];
        diff_N2_11_Plus += c[n - 1 - i] * N2_11_Plus_vals[nmax - i - 1];
        diff_N2_12_Plus += c[n - 1 - i] * N2_12_Plus_vals[nmax - i - 1];
        diff_N2_22_Plus += c[n - 1 - i] * N2_22_Plus_vals[nmax - i - 1];
    }
    Finite_Diff_PM.d2_I0_Plus_dnorm_f_2 = diff_I0_Plus;
    Finite_Diff_PM.d2_N1_1_Plus_dnorm_f_2 = diff_N1_1_Plus;
    Finite_Diff_PM.d2_N1_2_Plus_dnorm_f_2 = diff_N1_2_Plus;
    Finite_Diff_PM.d2_N2_11_Plus_dnorm_f_2 = diff_N2_11_Plus;
    Finite_Diff_PM.d2_N2_12_Plus_dnorm_f_2 = diff_N2_12_Plus;
    Finite_Diff_PM.d2_N2_22_Plus_dnorm_f_2 = diff_N2_22_Plus;
    
    // backward finite difference approximation to third derivative
    order = 1;
    n = order + prec;
    differ_backward ( h, order, prec, c, x );
    diff_I0_Plus = 0.0;
    diff_N1_1_Plus = 0.0;
    diff_N1_2_Plus = 0.0;
    diff_N2_11_Plus = 0.0;
    diff_N2_12_Plus = 0.0;
    diff_N2_22_Plus = 0.0;
    for (int i = 0; i < n; i++) {
        diff_I0_Plus += c[n - 1 - i] * E_Plus_vals[nmax - i - 1];
        diff_N1_1_Plus += c[n - 1 - i] * N1_1_Plus_vals[nmax - i - 1];
        diff_N1_2_Plus += c[n - 1 - i] * N1_2_Plus_vals[nmax - i - 1];
        diff_N2_11_Plus += c[n - 1 - i] * N2_11_Plus_vals[nmax - i - 1];
        diff_N2_12_Plus += c[n - 1 - i] * N2_12_Plus_vals[nmax - i - 1];
        diff_N2_22_Plus += c[n - 1 - i] * N2_22_Plus_vals[nmax - i - 1];
    }
    Finite_Diff_PM.dI0_Plus_dnorm_f = diff_I0_Plus;
    Finite_Diff_PM.dN1_1_Plus_dnorm_f = diff_N1_1_Plus;
    Finite_Diff_PM.dN1_2_Plus_dnorm_f = diff_N1_2_Plus;
    Finite_Diff_PM.dN2_11_Plus_dnorm_f = diff_N2_11_Plus;
    Finite_Diff_PM.dN2_12_Plus_dnorm_f = diff_N2_12_Plus;
    Finite_Diff_PM.dN2_22_Plus_dnorm_f = diff_N2_22_Plus;
    
    if (M1_State_temp.display) {
        cout << "dI0_Plus_dnorm_f = " << Finite_Diff_PM.dI0_Plus_dnorm_f << endl;
    }
    
    delete[] c;
    delete[] x;
    delete[] E_Plus_vals;
    delete[] N1_1_Plus_vals;
    delete[] N1_2_Plus_vals;
    delete[] N2_11_Plus_vals;
    delete[] N2_12_Plus_vals;
    delete[] N2_22_Plus_vals;
 }

void Calculate_Partial_Moments_Derivatives_Boundary(record_Partial_Moments &rec_PM, M1_State_Param &M1_State, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, const long double *x_SH, const long double *y_SH, const long double *z_SH, long double *Sk_final) {
    // Perform finite difference approximation to obtain derivative
    M1_State_Param M1_State_temp(M1_State);
    record_Partial_Moments rec_PM_Finite_Diff;
    long double *E_Plus_vals, *N1_1_Plus_vals, *N1_2_Plus_vals, *N2_11_Plus_vals, *N2_12_Plus_vals, *N2_22_Plus_vals;
    long double h;
    int order, prec, n;
    long double *c, *x;
    long double diff_I0_Plus, diff_N1_1_Plus, diff_N1_2_Plus, diff_N2_11_Plus, diff_N2_12_Plus, diff_N2_22_Plus;
    order = 1;
    prec = 6;
    n = order + prec;
    c = new long double[n];
    x = new long double[n];
    E_Plus_vals = new long double[n];
    N1_1_Plus_vals = new long double[n];
    N1_2_Plus_vals = new long double[n];
    N2_11_Plus_vals = new long double[n];
    N2_12_Plus_vals = new long double[n];
    N2_22_Plus_vals = new long double[n];
    h = finite_diff_h;
    
    M1_State_temp.f_test_knot = 1.0;
    M1_State_temp.h_f_test = h;
    
    M1_State_temp.flag_finite_difference_N1 = true;
    
    // backward finite difference approximation to first derivative
    differ_backward ( h, order, prec, c, x );
    
    for (int id_diff = 0; id_diff < n; id_diff++) {
        M1_State_temp.index_finite_diff = id_diff;
        M1_State_temp.x_finite_diff_N1 = x[id_diff];
        
        M1_State_temp.finite_difference_type_N1 = DERIVATIVE_BOUNDARY;
        
        NLOPT_Optimization_Algorithm_NON_GRAY_M1_3D(NULL, &rec_PM_Finite_Diff, 0, Var_index, num_points, M1_State_temp.id_proc, M1_State_temp, Mobius_Scale_Params, struc_map_N1, Sk_final, x_SH, y_SH, z_SH);
        
        E_Plus_vals[id_diff] = rec_PM_Finite_Diff.I0_Plus;
        N1_1_Plus_vals[id_diff] = rec_PM_Finite_Diff.N1_1_Plus;
        N1_2_Plus_vals[id_diff] = rec_PM_Finite_Diff.N1_2_Plus;
        N2_11_Plus_vals[id_diff] = rec_PM_Finite_Diff.N2_11_Plus;
        N2_12_Plus_vals[id_diff] = rec_PM_Finite_Diff.N2_12_Plus;
        N2_22_Plus_vals[id_diff] = rec_PM_Finite_Diff.N2_22_Plus;
    }
    
    diff_I0_Plus = 0.0;
    diff_N1_1_Plus = 0.0;
    diff_N1_2_Plus = 0.0;
    diff_N2_11_Plus = 0.0;
    diff_N2_12_Plus = 0.0;
    diff_N2_22_Plus = 0.0;
    for (int i = 0; i < n; i++) {
        diff_I0_Plus += c[n - 1 - i] * E_Plus_vals[n - 1 - i];
        diff_N1_1_Plus += c[n - 1 - i] * N1_1_Plus_vals[n - 1 - i];
        diff_N1_2_Plus += c[n - 1 - i] * N1_2_Plus_vals[n - 1 - i];
        diff_N2_11_Plus += c[n - 1 - i] * N2_11_Plus_vals[n - 1 - i];
        diff_N2_12_Plus += c[n - 1 - i] * N2_12_Plus_vals[n - 1 - i];
        diff_N2_22_Plus += c[n - 1 - i] * N2_22_Plus_vals[n - 1 - i];
    }
    rec_PM.dI0_Plus_dnorm_f = diff_I0_Plus;
    rec_PM.dN1_1_Plus_dnorm_f = diff_N1_1_Plus;
    rec_PM.dN1_2_Plus_dnorm_f = diff_N1_2_Plus;
    rec_PM.dN2_11_Plus_dnorm_f = diff_N2_11_Plus;
    rec_PM.dN2_12_Plus_dnorm_f = diff_N2_12_Plus;
    rec_PM.dN2_22_Plus_dnorm_f = diff_N2_22_Plus;
    
    if (M1_State_temp.display) {
        cout << "dI0_Plus_dnorm_f = " << rec_PM.dI0_Plus_dnorm_f << endl;
    }
    
    delete[] c;
    delete[] x;
    delete[] E_Plus_vals;
    delete[] N1_1_Plus_vals;
    delete[] N1_2_Plus_vals;
    delete[] N2_11_Plus_vals;
    delete[] N2_12_Plus_vals;
    delete[] N2_22_Plus_vals;
 }

void Calculate_Partial_Moments_Derivatives(record_Partial_Moments &rec_PM, const long double *x, const long double *Sk, M1_State_Param &M1_State) {
     long double *A_ONE, *H_ONE, *dIn_plus_d_lam;
     A_ONE = new long double[M1_State.NVARS*M1_State.NVARS];
     H_ONE = new long double[M1_State.NVARS*M1_State.NVARS];
     dIn_plus_d_lam = new long double[M1_State.NVARS];
     
     M1_State.set_x(x);
     M1_State.set_Sk(Sk);
     M1_State.set_Sk_cur(Sk);
     
     Compute_A_ONE_Matrix(A_ONE, x, Sk, M1_State);
     Compute_H_ONE_Matrix(H_ONE, x, Sk, M1_State);
     
     long double temp_val_E, temp_val_N1_1, temp_val_N1_2, temp_val_N1_3;
     long double diff_temp_val_E, diff_temp_val_N1_1, diff_temp_val_N1_2, diff_temp_val_N1_3;
     for (int i = 0; i < M1_State.NVARS; i++) {
         temp_val_E = 0.0;
         temp_val_N1_1 = 0.0;
         temp_val_N1_2 = 0.0;
         temp_val_N1_3 = 0.0;
         for (int j = 0; j < M1_State.NVARS; j++) {
             temp_val_E += H_ONE[0*M1_State.NVARS+j]*A_ONE[j*M1_State.NVARS+i];
             temp_val_N1_1 += H_ONE[1*M1_State.NVARS+j]*A_ONE[j*M1_State.NVARS+i];
             temp_val_N1_2 += H_ONE[2*M1_State.NVARS+j]*A_ONE[j*M1_State.NVARS+i];
             temp_val_N1_3 += H_ONE[3*M1_State.NVARS+j]*A_ONE[j*M1_State.NVARS+i];
        }
        diff_temp_val_E = fabs(temp_val_E - dIn_dInm1(0, i, M1_State.I0, M1_State.N1_1, M1_State.N1_2, M1_State.N1_3));
        diff_temp_val_N1_1 = fabs(temp_val_N1_1 - dIn_dInm1(1, i, M1_State.I0, M1_State.N1_1, M1_State.N1_2, M1_State.N1_3));
        diff_temp_val_N1_2 = fabs(temp_val_N1_2 - dIn_dInm1(2, i, M1_State.I0, M1_State.N1_1, M1_State.N1_2, M1_State.N1_3));
        diff_temp_val_N1_3 = fabs(temp_val_N1_3 - dIn_dInm1(3, i, M1_State.I0, M1_State.N1_1, M1_State.N1_2, M1_State.N1_3));
        
        if (diff_temp_val_E > 1.0e-6 || diff_temp_val_N1_1 > 1.0e-6 || diff_temp_val_N1_2 > 1.0e-6 || diff_temp_val_N1_3 > 1.0e-6) {
            cout << "temp_val_E = " << temp_val_E << "   " << "temp_val_N1_1 = " << temp_val_N1_1 << "   " << "temp_val_N1_2 = " << temp_val_N1_2 << "   " << "temp_val_N1_3 = " << temp_val_N1_3 << endl;
            
            cout << "diff_temp_val_E = " << diff_temp_val_E << "   " << "diff_temp_val_N1_1 = " << diff_temp_val_N1_1 << "   " << "diff_temp_val_N1_2 = " << diff_temp_val_N1_2 << "   " << "diff_temp_val_N1_3 = " << diff_temp_val_N1_3 << endl;
            exit(0);
        }
     }
     
     M1_State.Index_Partial_Moments_Derivatives = 0;
     Partial_Quadrature_3D_Func(dIn_plus_d_lam, M1_State.NVARS, Angular_Moments_First_Derivatives, &M1_State);
     
     rec_PM.dI0_Plus_dN1_1 = 0.0;
     rec_PM.dI0_Plus_dN1_2 = 0.0;
     rec_PM.dI0_Plus_dN1_3 = 0.0;
     
     for (int i = 0; i < M1_State.NVARS; i++) {
         rec_PM.dI0_Plus_dN1_1 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 1];
         rec_PM.dI0_Plus_dN1_2 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 2];
         rec_PM.dI0_Plus_dN1_3 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 3];
     }
     
     rec_PM.d2I0_Plus_dN1_1 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 1, 1, M1_State);
     rec_PM.d2I0_Plus_dN1_2 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 2, 2, M1_State);
     rec_PM.d2I0_Plus_dN1_3 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 3, 3, M1_State);
     
     if (M1_State.display) {
         cout << "dI0_Plus_dN1_1 = " << rec_PM.dI0_Plus_dN1_1 << endl;
         cout << "********************************************************************" << endl;
         cout << endl;
     }
     
//      switch (M1_State.N1_1_VAL_DOMAIN) {
//          case FULL_3D:
             rec_PM.dI0_Plus_dnorm_f = Cartesian_to_spherical_Coordinates_Jacobian(rec_PM.N1, 
                                                                                   rec_PM.x_val, rec_PM.y_val, rec_PM.z_val, 
                                                                                   rec_PM.dI0_Plus_dN1_1, 
                                                                                   rec_PM.dI0_Plus_dN1_2, 
                                                                                   rec_PM.dI0_Plus_dN1_3, 
                                                                                   0);
             
             rec_PM.d2I0_Plus_dnorm_f_2 = Cartesian_to_spherical_Coordinates_Jacobian_Second_Derivatives(rec_PM.N1, 
                                                                                                         rec_PM.x_val, 
                                                                                                         rec_PM.y_val, 
                                                                                                         rec_PM.z_val, 
                                                                                                         rec_PM.d2I0_Plus_dN1_1, 
                                                                                                         rec_PM.d2I0_Plus_dN1_2, 
                                                                                                         rec_PM.d2I0_Plus_dN1_3, 
                                                                                                         0);
//              break;
//          case N1_1_EQUAL_ZERO:
//              // Then N1_1 not included
//              rec_PM.dI0_Plus_dnorm_f = Cartesian_to_spherical_Coordinates_Jacobian(rec_PM.N1, 
//                                                                                    rec_PM.x_val, rec_PM.y_val, rec_PM.z_val, 
//                                                                                    0.0, 
//                                                                                    rec_PM.dI0_Plus_dN1_2, 
//                                                                                    rec_PM.dI0_Plus_dN1_3, 
//                                                                                    0);
//              rec_PM.d2I0_Plus_dnorm_f_2 = Cartesian_to_spherical_Coordinates_Jacobian_Second_Derivatives(rec_PM.N1, 
//                                                                                                          rec_PM.x_val, 
//                                                                                                          rec_PM.y_val, 
//                                                                                                          rec_PM.z_val, 
//                                                                                                          0.0, 
//                                                                                                          rec_PM.d2I0_Plus_dN1_2, 
//                                                                                                          rec_PM.d2I0_Plus_dN1_3, 
//                                                                                                          0);
//              break;
//          default:
//              cout << "N1_1 Domain Type not Specified for dI0_Plus_dnorm_f !!!!!!!!!!!!!!!!!!!!!" << endl;
//              exit(0);
//              break;
//          
//      }
     
     M1_State.Index_Partial_Moments_Derivatives = 1;
     Partial_Quadrature_3D_Func(dIn_plus_d_lam, M1_State.NVARS, Angular_Moments_First_Derivatives, &M1_State);
     
     rec_PM.dN1_1_Plus_dN1_1 = 0.0;
     rec_PM.dN1_1_Plus_dN1_2 = 0.0;
     rec_PM.dN1_1_Plus_dN1_3 = 0.0;
     
     for (int i = 0; i < M1_State.NVARS; i++) {
         rec_PM.dN1_1_Plus_dN1_1 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 1];
         rec_PM.dN1_1_Plus_dN1_2 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 2];
         rec_PM.dN1_1_Plus_dN1_3 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 3];
     }
     
     rec_PM.d2N1_1_Plus_dN1_1 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 1, 1, M1_State);
     rec_PM.d2N1_1_Plus_dN1_2 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 2, 2, M1_State);
     rec_PM.d2N1_1_Plus_dN1_3 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 3, 3, M1_State);
     
//      switch (M1_State.N1_1_VAL_DOMAIN) {
//          case FULL_3D:
             rec_PM.dN1_1_Plus_dnorm_f = Cartesian_to_spherical_Coordinates_Jacobian(rec_PM.N1, 
                                                                                     rec_PM.x_val, rec_PM.y_val, rec_PM.z_val, 
                                                                                     rec_PM.dN1_1_Plus_dN1_1, rec_PM.dN1_1_Plus_dN1_2, rec_PM.dN1_1_Plus_dN1_3, 0);
             
             rec_PM.d2N1_1_Plus_dnorm_f_2 = Cartesian_to_spherical_Coordinates_Jacobian_Second_Derivatives(rec_PM.N1, 
                                                                                                           rec_PM.x_val, 
                                                                                                           rec_PM.y_val, 
                                                                                                           rec_PM.z_val, 
                                                                                                           rec_PM.d2N1_1_Plus_dN1_1, 
                                                                                                           rec_PM.d2N1_1_Plus_dN1_2, 
                                                                                                           rec_PM.d2N1_1_Plus_dN1_3, 
                                                                                                           0);
//              break;
//          case N1_1_EQUAL_ZERO:
//              // Then N1_1 not included
//              rec_PM.dN1_1_Plus_dnorm_f = Cartesian_to_spherical_Coordinates_Jacobian(rec_PM.N1, 
//                                                                                      rec_PM.x_val, rec_PM.y_val, rec_PM.z_val, 
//                                                                                      0.0, rec_PM.dN1_1_Plus_dN1_2, rec_PM.dN1_1_Plus_dN1_3, 0);
//              
//              rec_PM.d2N1_1_Plus_dnorm_f_2 = Cartesian_to_spherical_Coordinates_Jacobian_Second_Derivatives(rec_PM.N1, 
//                                                                                                            rec_PM.x_val, 
//                                                                                                            rec_PM.y_val, 
//                                                                                                            rec_PM.z_val, 
//                                                                                                            0.0, 
//                                                                                                            rec_PM.d2N1_1_Plus_dN1_2, 
//                                                                                                            rec_PM.d2N1_1_Plus_dN1_3, 
//                                                                                                            0);
//              break;
//          default:
//              cout << "N1_1 Domain Type not Specified for dI0_Plus_dnorm_f !!!!!!!!!!!!!!!!!!!!!" << endl;
//              exit(0);
//              break;
//          
//      }
     
     if (M1_State.display) {
         cout << "dN1_1_Plus_dN1_1 = " << rec_PM.dN1_1_Plus_dN1_1 << " " << "dN1_1_Plus_dN1_2 = " << rec_PM.dN1_1_Plus_dN1_2 << " " << "dN1_1_Plus_dN1_3 = " << rec_PM.dN1_1_Plus_dN1_3 << endl;
         cout << "********************************************************************" << endl;
         cout << endl;
     }
     
     M1_State.Index_Partial_Moments_Derivatives = 2;
     Partial_Quadrature_3D_Func(dIn_plus_d_lam, M1_State.NVARS, Angular_Moments_First_Derivatives, &M1_State);
     
     rec_PM.dN1_2_Plus_dN1_1 = 0.0;
     rec_PM.dN1_2_Plus_dN1_2 = 0.0;
     rec_PM.dN1_2_Plus_dN1_3 = 0.0;
     
     for (int i = 0; i < M1_State.NVARS; i++) {
         rec_PM.dN1_2_Plus_dN1_1 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 1];
         rec_PM.dN1_2_Plus_dN1_2 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 2];
         rec_PM.dN1_2_Plus_dN1_3 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 3];
     }
     
     rec_PM.d2N1_2_Plus_dN1_1 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 1, 1, M1_State);
     rec_PM.d2N1_2_Plus_dN1_2 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 2, 2, M1_State);
     rec_PM.d2N1_2_Plus_dN1_3 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 3, 3, M1_State);
     
//      switch (M1_State.N1_1_VAL_DOMAIN) {
//          case FULL_3D:
             rec_PM.dN1_2_Plus_dnorm_f = Cartesian_to_spherical_Coordinates_Jacobian(rec_PM.N1, 
                                                                                     rec_PM.x_val, rec_PM.y_val, rec_PM.z_val, 
                                                                                     rec_PM.dN1_2_Plus_dN1_1, rec_PM.dN1_2_Plus_dN1_2, rec_PM.dN1_2_Plus_dN1_3, 
                                                                                     0);
             
             rec_PM.d2N1_2_Plus_dnorm_f_2 = Cartesian_to_spherical_Coordinates_Jacobian_Second_Derivatives(rec_PM.N1, 
                                                                                                           rec_PM.x_val, 
                                                                                                           rec_PM.y_val, 
                                                                                                           rec_PM.z_val, 
                                                                                                           rec_PM.d2N1_2_Plus_dN1_1, 
                                                                                                           rec_PM.d2N1_2_Plus_dN1_2, 
                                                                                                           rec_PM.d2N1_2_Plus_dN1_3, 
                                                                                                           0);
//              break;
//          case N1_1_EQUAL_ZERO:
//              // Then N1_1 not included
//              rec_PM.dN1_2_Plus_dnorm_f = Cartesian_to_spherical_Coordinates_Jacobian(rec_PM.N1, 
//                                                                                      rec_PM.x_val, rec_PM.y_val, rec_PM.z_val, 
//                                                                                      0.0, rec_PM.dN1_2_Plus_dN1_2, rec_PM.dN1_2_Plus_dN1_3, 
//                                                                                      0);
//              
//              rec_PM.d2N1_2_Plus_dnorm_f_2 = Cartesian_to_spherical_Coordinates_Jacobian_Second_Derivatives(rec_PM.N1, 
//                                                                                                            rec_PM.x_val, 
//                                                                                                            rec_PM.y_val, 
//                                                                                                            rec_PM.z_val, 
//                                                                                                            0.0, 
//                                                                                                            rec_PM.d2N1_2_Plus_dN1_2, 
//                                                                                                            rec_PM.d2N1_2_Plus_dN1_3, 
//                                                                                                            0);
//              break;
//          default:
//              cout << "N1_1 Domain Type not Specified for dI0_Plus_dnorm_f !!!!!!!!!!!!!!!!!!!!!" << endl;
//              exit(0);
//              break;
//          
//      }
     
     if (M1_State.display) {
         cout << "dN1_2_Plus_dN1_1 = " << rec_PM.dN1_2_Plus_dN1_1 << endl;
         cout << "********************************************************************" << endl;
         cout << endl;
     }
     
     M1_State.Index_Partial_Moments_Derivatives = 3;
     Partial_Quadrature_3D_Func(dIn_plus_d_lam, M1_State.NVARS, Angular_Moments_First_Derivatives, &M1_State);
     
     rec_PM.dN2_11_Plus_dN1_1 = 0.0;
     rec_PM.dN2_11_Plus_dN1_2 = 0.0;
     rec_PM.dN2_11_Plus_dN1_3 = 0.0;
     
     for (int i = 0; i < M1_State.NVARS; i++) {
         rec_PM.dN2_11_Plus_dN1_1 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 1];
         rec_PM.dN2_11_Plus_dN1_2 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 2];
         rec_PM.dN2_11_Plus_dN1_3 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 3];
     }
     
     rec_PM.d2N2_11_Plus_dN1_1 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 1, 1, M1_State);
     rec_PM.d2N2_11_Plus_dN1_2 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 2, 2, M1_State);
     rec_PM.d2N2_11_Plus_dN1_3 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 3, 3, M1_State);
     
//      switch (M1_State.N1_1_VAL_DOMAIN) {
//          case FULL_3D:
             rec_PM.dN2_11_Plus_dnorm_f = Cartesian_to_spherical_Coordinates_Jacobian(rec_PM.N1, 
                                                                                      rec_PM.x_val, rec_PM.y_val, rec_PM.z_val, 
                                                                                      rec_PM.dN2_11_Plus_dN1_1, 
                                                                                      rec_PM.dN2_11_Plus_dN1_2, 
                                                                                      rec_PM.dN2_11_Plus_dN1_3, 
                                                                                      0);
             
             rec_PM.d2N2_11_Plus_dnorm_f_2 = Cartesian_to_spherical_Coordinates_Jacobian_Second_Derivatives(rec_PM.N1, 
                                                                                                            rec_PM.x_val, 
                                                                                                            rec_PM.y_val, 
                                                                                                            rec_PM.z_val, 
                                                                                                            rec_PM.d2N2_11_Plus_dN1_1, 
                                                                                                            rec_PM.d2N2_11_Plus_dN1_2, 
                                                                                                            rec_PM.d2N2_11_Plus_dN1_3, 
                                                                                                            0);
//              break;
//          case N1_1_EQUAL_ZERO:
//              // Then N1_1 not included
//              rec_PM.dN2_11_Plus_dnorm_f = Cartesian_to_spherical_Coordinates_Jacobian(rec_PM.N1, 
//                                                                                       rec_PM.x_val, rec_PM.y_val, rec_PM.z_val, 
//                                                                                       0.0, 
//                                                                                       rec_PM.dN2_11_Plus_dN1_2, 
//                                                                                       rec_PM.dN2_11_Plus_dN1_3, 
//                                                                                       0);
//              
//              rec_PM.d2N2_11_Plus_dnorm_f_2 = Cartesian_to_spherical_Coordinates_Jacobian_Second_Derivatives(rec_PM.N1, 
//                                                                                                             rec_PM.x_val, 
//                                                                                                             rec_PM.y_val, 
//                                                                                                             rec_PM.z_val, 
//                                                                                                             0.0, 
//                                                                                                             rec_PM.d2N2_11_Plus_dN1_2, 
//                                                                                                             rec_PM.d2N2_11_Plus_dN1_3, 
//                                                                                                             0);
//              break;
//          default:
//              cout << "N1_1 Domain Type not Specified for dI0_Plus_dnorm_f !!!!!!!!!!!!!!!!!!!!!" << endl;
//              exit(0);
//              break;
//          
//      }
     
     if (M1_State.display) {
         cout << "dN2_11_Plus_dN1_1 = " << rec_PM.dN2_11_Plus_dN1_1 << endl;
         cout << "********************************************************************" << endl;
         cout << endl;
     }
     
     M1_State.Index_Partial_Moments_Derivatives = 4;
     Partial_Quadrature_3D_Func(dIn_plus_d_lam, M1_State.NVARS, Angular_Moments_First_Derivatives, &M1_State);
     
     rec_PM.dN2_12_Plus_dN1_1 = 0.0;
     rec_PM.dN2_12_Plus_dN1_2 = 0.0;
     rec_PM.dN2_12_Plus_dN1_3 = 0.0;
     
     for (int i = 0; i < M1_State.NVARS; i++) {
         rec_PM.dN2_12_Plus_dN1_1 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 1];
         rec_PM.dN2_12_Plus_dN1_2 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 2];
         rec_PM.dN2_12_Plus_dN1_3 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 3];
     }
     
     rec_PM.d2N2_12_Plus_dN1_1 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 1, 1, M1_State);
     rec_PM.d2N2_12_Plus_dN1_2 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 2, 2, M1_State);
     rec_PM.d2N2_12_Plus_dN1_3 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 3, 3, M1_State);
     
//      switch (M1_State.N1_1_VAL_DOMAIN) {
//          case FULL_3D:
             rec_PM.dN2_12_Plus_dnorm_f = Cartesian_to_spherical_Coordinates_Jacobian(rec_PM.N1, 
                                                                                      rec_PM.x_val, rec_PM.y_val, rec_PM.z_val, 
                                                                                      rec_PM.dN2_12_Plus_dN1_1, rec_PM.dN2_12_Plus_dN1_2, rec_PM.dN2_12_Plus_dN1_3, 
                                                                                      0);
             
             rec_PM.d2N2_12_Plus_dnorm_f_2 = Cartesian_to_spherical_Coordinates_Jacobian_Second_Derivatives(rec_PM.N1,
                                                                                                            rec_PM.x_val, 
                                                                                                            rec_PM.y_val, 
                                                                                                            rec_PM.z_val, 
                                                                                                            rec_PM.d2N2_12_Plus_dN1_1, 
                                                                                                            rec_PM.d2N2_12_Plus_dN1_2, 
                                                                                                            rec_PM.d2N2_12_Plus_dN1_3, 
                                                                                                            0);
//              break;
//          case N1_1_EQUAL_ZERO:
//              // Then N1_1 not included
//              rec_PM.dN2_12_Plus_dnorm_f = Cartesian_to_spherical_Coordinates_Jacobian(rec_PM.N1, 
//                                                                                       rec_PM.x_val, rec_PM.y_val, rec_PM.z_val, 
//                                                                                       0.0, rec_PM.dN2_12_Plus_dN1_2, rec_PM.dN2_12_Plus_dN1_3, 
//                                                                                       0);
//              
//              rec_PM.d2N2_12_Plus_dnorm_f_2 = Cartesian_to_spherical_Coordinates_Jacobian_Second_Derivatives(rec_PM.N1, 
//                                                                                                             rec_PM.x_val, 
//                                                                                                             rec_PM.y_val, 
//                                                                                                             rec_PM.z_val, 
//                                                                                                             0.0, 
//                                                                                                             rec_PM.d2N2_12_Plus_dN1_2, 
//                                                                                                             rec_PM.d2N2_12_Plus_dN1_3, 
//                                                                                                             0);
//              break;
//          default:
//              cout << "N1_1 Domain Type not Specified for dI0_Plus_dnorm_f !!!!!!!!!!!!!!!!!!!!!" << endl;
//              exit(0);
//              break;
//          
//      }
     
     if (M1_State.display) {
         cout << "dN2_12_Plus_dN1_1 = " << rec_PM.dN2_12_Plus_dN1_1 << endl;
         cout << "********************************************************************" << endl;
         cout << endl;
     }
     
     M1_State.Index_Partial_Moments_Derivatives = 5;
     Partial_Quadrature_3D_Func(dIn_plus_d_lam, M1_State.NVARS, Angular_Moments_First_Derivatives, &M1_State);
     
     rec_PM.dN2_22_Plus_dN1_1 = 0.0;
     rec_PM.dN2_22_Plus_dN1_2 = 0.0;
     rec_PM.dN2_22_Plus_dN1_3 = 0.0;
     
     for (int i = 0; i < M1_State.NVARS; i++) {
         rec_PM.dN2_22_Plus_dN1_1 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 1];
         rec_PM.dN2_22_Plus_dN1_2 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 2];
         rec_PM.dN2_22_Plus_dN1_3 += dIn_plus_d_lam[i]*A_ONE[i*M1_State.NVARS + 3];
     }
     
     rec_PM.d2N2_22_Plus_dN1_1 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 1, 1, M1_State);
     rec_PM.d2N2_22_Plus_dN1_2 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 2, 2, M1_State);
     rec_PM.d2N2_22_Plus_dN1_3 = Calculate_Angular_Moments_Second_Derivatives(dIn_plus_d_lam, A_ONE, 3, 3, M1_State);
     
//      switch (M1_State.N1_1_VAL_DOMAIN) {
//          case FULL_3D:
             rec_PM.dN2_22_Plus_dnorm_f = Cartesian_to_spherical_Coordinates_Jacobian(rec_PM.N1, 
                                                                                      rec_PM.x_val, rec_PM.y_val, rec_PM.z_val, 
                                                                                      rec_PM.dN2_22_Plus_dN1_1, rec_PM.dN2_22_Plus_dN1_2, rec_PM.dN2_22_Plus_dN1_3, 
                                                                                      0);
             
             rec_PM.d2N2_22_Plus_dnorm_f_2 = Cartesian_to_spherical_Coordinates_Jacobian_Second_Derivatives(rec_PM.N1, 
                                                                                                            rec_PM.x_val, 
                                                                                                            rec_PM.y_val, 
                                                                                                            rec_PM.z_val, 
                                                                                                            rec_PM.d2N2_22_Plus_dN1_1, 
                                                                                                            rec_PM.d2N2_22_Plus_dN1_2, 
                                                                                                            rec_PM.d2N2_22_Plus_dN1_3, 
                                                                                                            0);
//              break;
//          case N1_1_EQUAL_ZERO:
//              // Then N1_1 not included
//              rec_PM.dN2_22_Plus_dnorm_f = Cartesian_to_spherical_Coordinates_Jacobian(rec_PM.N1, 
//                                                                                       rec_PM.x_val, rec_PM.y_val, rec_PM.z_val, 
//                                                                                       0.0, rec_PM.dN2_22_Plus_dN1_2, rec_PM.dN2_22_Plus_dN1_3, 
//                                                                                       0);
//              
//              rec_PM.d2N2_22_Plus_dnorm_f_2 = Cartesian_to_spherical_Coordinates_Jacobian_Second_Derivatives(rec_PM.N1, 
//                                                                                                             rec_PM.x_val, 
//                                                                                                             rec_PM.y_val, 
//                                                                                                             rec_PM.z_val, 
//                                                                                                             0.0, 
//                                                                                                             rec_PM.d2N2_22_Plus_dN1_2, 
//                                                                                                             rec_PM.d2N2_22_Plus_dN1_3, 
//                                                                                                             0);
//              break;
//          default:
//              cout << "N1_1 Domain Type not Specified for dI0_Plus_dnorm_f !!!!!!!!!!!!!!!!!!!!!" << endl;
//              exit(0);
//              break;
//          
//      }
     
     if (M1_State.display) {
         cout << "dN2_22_Plus_dN1_1 = " << rec_PM.dN2_22_Plus_dN1_1 << endl;
         cout << "********************************************************************" << endl;
         cout << endl;
     }
     
     delete[] A_ONE;
     delete[] H_ONE;
     delete[] dIn_plus_d_lam;
 }
 
 void Partial_Moments_M1_NON_GRAY_3D(long double *PARTIAL_MOMS, const int &NFUN, const long double *x, const long double *Sk, M1_State_Param &M1_State) {
    M1_State.set_x(x);
    M1_State.set_Sk(Sk);
    
    Partial_Quadrature_3D_Func(PARTIAL_MOMS, NFUN, Partial_Moments_3D, &M1_State);
    
    for (int i = 0; i < NFUN; i++) {
        PARTIAL_MOMS[i] /= M1_State.I0;
    }
    if (PARTIAL_MOMS[0] < 0.0 || fabs(PARTIAL_MOMS[0]) > 1.0+1.0e-4) {
        cout << endl;
        cout << "***************************" << "M1_State.I0 = " << M1_State.I0 << "      " << "E_Plus = " << PARTIAL_MOMS[0] << endl;
        cout << endl;
        exit(0);
    }
}

int Partial_Moments_3D(long double *Partial_Moms, const int &NFUN, const int &Id_angle, void *fdata) {
     long double coeff_1;
     M1_State_Param *M1_State = (M1_State_Param *) fdata;
     
     coeff_1 = 1.0;
     
     long double *poly, *poly_Sk, *poly_PM_basis;
     poly = new long double[M1_State->NVARS];
     poly_Sk = new long double[M1_State->NVARS];
     poly_PM_basis = new long double[M1_State->NFUNS_PM_3D];
     
     generate_polynomials_3D(poly, M1_State->Omega1[Id_angle], M1_State->Omega2[Id_angle], M1_State->Omega3[Id_angle]);
     
     for (int i = 0; i < M1_State->NVARS; i++) {
         poly_Sk[i] = 0.0;
         for (int j = 0; j < M1_State->NVARS; j++) {
             poly_Sk[i] += M1_State->Sk[i*M1_State->NVARS+j]*poly[j];
        }
    }
    
    for (int i = 0; i < M1_State->NFUNS_PM_3D; i++) {
        M1_State->Index_Partial_Moments_Derivatives = i;
        poly_PM_basis[i] = generate_partial_moments_monomials_basis(Id_angle, *M1_State);
    }
    
    switch (M1_State->Regime) {
        case GRAY:
            for (int i = 0; i < M1_State->NFUNS_PM_3D; i++) {
                Partial_Moms[i] = poly_PM_basis[i]/(pow(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 4));
            }
            break;
        case BOSE_EINSTEIN:
            for (int i = 0; i < M1_State->NFUNS_PM_3D; i++) {
                Partial_Moms[i] = coeff_1*poly_PM_basis[i]/(exp(-sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS)) - 1.0);
            }
            break;
        case HYPERBOLIC_LIMIT:
            for (int i = 0; i < M1_State->NFUNS_PM_3D; i++) {
                Partial_Moms[i] = coeff_1*poly_PM_basis[i]*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
            }
            break;
        case LOGARITHMIC_LIMIT:
            for (int i = 0; i < M1_State->NFUNS_PM_3D; i++) {
                Partial_Moms[i] = coeff_1*poly_PM_basis[i]/(-sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
            }
            break;
    }
     
     delete[] poly;
     delete[] poly_Sk;
     delete[] poly_PM_basis;
     return 0;
 }
