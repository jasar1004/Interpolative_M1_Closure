#ifndef _M1_DATA_1D_CHEBY_H_INCLUDED
#include "M1_Data_1D_Cheby.h"
#endif // _M1_DATA_1D_CHEBY_H_INCLUDED

// ******************************************************************************
// This routine sets up the least-squares algorithm for computing optimal values 
// of the length scale L_Chi2 of the exponential mapping for the zeroth-order 
// angular moment I^{(0)}
// ******************************************************************************
void M1_Data_1D_Cheby :: Nested_Least_Squares_Optimization_L_Chi2(M1_1D_Data_Pointer &M1_1D_Data) {
    Least_Squares_L_Chi2_Mode = LEAST_SQUARES_L_CHI2_OFF;
    int NVARS = N_pts_L_Chi2_LS;
    long double x[NVARS];
    long double tol_x[NVARS];
    M1_Data_1D_Cheby Chi2_M1_1D_Current(N_Points_E, N_Points_f, N_pts_Mob_Scale, N_pts_L_Chi2_LS);
    
    long double minf = 1.0e6, min_grad = 1.0e6; /* the minimum objective value, upon return */
    nlopt_opt opt = NULL;
    
    for (int i = 0; i < NVARS; i++) {
        x[i] = 0.0;
        tol_x[i] = 1.0e-12;
    }
    
    x[0] = 1.0;
    
//     x[0] = 4.018234;
//     x[1] = 0.704705;
//     x[2] = 0.334077;
//     x[3] = -0.240986;
//     x[4] = 0.652468;
    
    opt = nlopt_create(NLOPT_LN_COBYLA, NVARS); /* algorithm and dimensionality */
    
    Copy_to(&Chi2_M1_1D_Current);
    
    M1_1D_Data.M1_Data_BE = &Chi2_M1_1D_Current;
    
    nlopt_set_maxeval(opt, 40);
    
    nlopt_set_min_objective(opt, myfunc_Nested_Least_Squares_L_Chi2, &M1_1D_Data);
    
    nlopt_set_xtol_abs(opt, tol_x);
    nlopt_set_ftol_abs(opt, 1.0e-12);
    
    my_constraint_data data[1];
    add_constraints_Nested_Least_Squares_L_Chi2(opt, data, &M1_1D_Data);
                                
    if (nlopt_optimize(opt, x, &minf, &min_grad) < 0) {
        printf("....................nlopt failed!.....................\n");
        exit(0);
    } else {
        if (id_proc == 0) {
            printf("***********found minimum at f(");
            for (int index_vars = 0; index_vars < NVARS; index_vars++) {
                if (index_vars < NVARS - 1) {
                    printf("%Lf,", x[index_vars]);
                } else {
                    printf("%Lf", x[index_vars]);
                }
            }
            printf(") = %0.10Le***********\n ", minf);
        }
    }
    
    // Update coefficients for the approximation of optimal length scale L_Chi2
    for (int i = 0; i < N_pts_L_Chi2_LS; i++) {
        Coefficient_Matrix_Fit_Mob_Scale[i] = x[i]; //(M1_1D_Data.M1_Data_BE)->Coefficient_Matrix_Fit_Mob_Scale[i];
    }
    
//     // Rewrite Coefficients in Original Basis
//     Chebyshev_First_Kind_to_Monomial_Basis_Norm_f(Coefficient_Matrix_Fit_Mob_Scale, 1, N_pts_L_Chi2_LS, 1);
//     printf("*********** Coefficients in Monomial Basis (");
//     for (int index_vars = 0; index_vars < NVARS; index_vars++) {
//         if (index_vars < NVARS - 1) {
//             printf("%Lf,", Coefficient_Matrix_Fit_Mob_Scale[index_vars]);
//         } else {
//             printf("%Lf", Coefficient_Matrix_Fit_Mob_Scale[index_vars]);
//         }
//     }
//     printf(")***********\n ");
    Least_Squares_L_Chi2_Mode = LEAST_SQUARES_L_CHI2_OFF;
}
 long double myfunc_Nested_Least_Squares_L_Chi2(unsigned n, const long double *x, long double *grad, void *my_func_data) {
    ofstream out_L_inf_Norm_temp;
    M1_1D_Data_Pointer *Chi2_M1_1D_Data_Pointer = (M1_1D_Data_Pointer *) my_func_data;
    
    M1_Data_1D_Cheby *Chi2_M1_1D_Cheby = Chi2_M1_1D_Data_Pointer->M1_Data_BE;
    M1_Data_1D_Cheby *Chi2_M1_1D_HL = Chi2_M1_1D_Data_Pointer->M1_Data_HL;
    M1_Data_1D_Cheby *Chi2_M1_1D_LL = Chi2_M1_1D_Data_Pointer->M1_Data_LL;
    M1_Data_1D_Cheby *Chi2_M1_1D_Uniform = Chi2_M1_1D_Data_Pointer->M1_Data_Uniform_BE;
    M1_Data_1D_Cheby *Chi2_M1_1D_Uniform_HL = Chi2_M1_1D_Data_Pointer->M1_Data_Uniform_HL;
    M1_Data_1D_Cheby *Chi2_M1_1D_Uniform_LL = Chi2_M1_1D_Data_Pointer->M1_Data_Uniform_LL;
    
    //**********************************************
    // Compute Objective Function
    //**********************************************
    long double objective_function;
    int NF_obj = 1;
    
    for (int i = 0; i < Chi2_M1_1D_Cheby->N_pts_L_Chi2_LS; i++) {
        Chi2_M1_1D_Cheby->Coefficient_Matrix_Fit_Mob_Scale[i] = x[i];
        cout << "i = " << i << "  " << "Coefficient_Matrix_Fit_Mob_Scale = " << x[i] << endl;
    }
    
    // Now perform polynomial interpolation
    Chi2_M1_1D_Cheby->Polynomial_Interpolation_BE(*Chi2_M1_1D_HL, *Chi2_M1_1D_LL);
    
    // Chi2_M1_1D_Uniform->SetupInterpolant_Values_BE(*Chi2_M1_1D_Uniform_HL, *Chi2_M1_1D_Uniform_LL, true, Chi2_M1_1D_Cheby);
    
    // Now assess the error of our interpolative-based approximation of the Eddington factor based on the curremt
    // iterate for the coefficients of the interpolative-based approximation of the optimal length scale for the 
    // exponential mapping of the radiative energy density
    Chi2_M1_1D_Cheby->Compute_L_ONE_L_TWO_Errors_Nested_Least_Squares_L_Chi2(Chi2_M1_1D_Uniform, 0, out_L_inf_Norm_temp, 0);
    
    objective_function = Chi2_M1_1D_Uniform->L2_Norm_Chi2;
    objective_function = objective_function*objective_function;
    // cout << "objective_function = " << objective_function << endl;
    
    return objective_function;   
}

void add_constraints_Nested_Least_Squares_L_Chi2(nlopt_opt &opt, my_constraint_data *data, M1_1D_Data_Pointer *data_realiz) {
    nlopt_add_inequality_constraint(opt, myconstraint_Nested_Least_Squares_L_Chi2_realizability_Upper_Bound, data_realiz, 0.0);
    nlopt_add_inequality_constraint(opt, myconstraint_Nested_Least_Squares_L_Chi2_realizability_Lower_Bound, data_realiz, 0.0);
    
//     nlopt_add_inequality_constraint(opt, myconstraint_Nested_Least_Squares_L_Chi2_Hyperbolicity_Lambda_s, data_realiz, 0.0);
//     // nlopt_add_inequality_constraint(opt, myconstraint_Nested_Least_Squares_L_Chi2_Speed_Limit, data_realiz, 0.0);
//     
    nlopt_add_inequality_constraint(opt, myconstraint_Nested_Least_Squares_L_Chi2_Hyperbolicity_Lambda_s_2D, data_realiz, 0.0);
}
 
long double myconstraint_Nested_Least_Squares_L_Chi2_realizability_Upper_Bound(unsigned n, const long double *x, long double *grad, void *data) {
    long double constr_val;
    M1_1D_Data_Pointer *Chi2_M1_1D_Data_Pointer = (M1_1D_Data_Pointer *) data;
    
    M1_Data_1D_Cheby *Chi2_M1_1D_Current = Chi2_M1_1D_Data_Pointer->M1_Data_BE;
    M1_Data_1D_Cheby *Chi2_M1_1D_HL = Chi2_M1_1D_Data_Pointer->M1_Data_HL;
    M1_Data_1D_Cheby *Chi2_M1_1D_LL = Chi2_M1_1D_Data_Pointer->M1_Data_LL;
    M1_Data_1D_Cheby *Chi2_M1_1D_Uniform = Chi2_M1_1D_Data_Pointer->M1_Data_Uniform_BE;
    
//     for (int i = 0; i < Chi2_M1_1D_Current->N_pts_L_Chi2_LS; i++) {
//         Chi2_M1_1D_Current->Coefficient_Matrix_Fit_Mob_Scale[i] = x[i];
//     }
    
//     Chi2_M1_1D_Current->Polynomial_Interpolation_BE(*Chi2_M1_1D_HL, *Chi2_M1_1D_LL);
    
    // Compute the nodal distribution for the length scale L_Chi2
    double Mobius_Scale_Actual, Mobius_Scale_Fit;
    int N_pts_Mob_Scale = 100;
    long double alpha_Lag;
    alpha_Lag = 0;
    long double x_L_Chi2[N_pts_Mob_Scale], w_L_Chi2[N_pts_Mob_Scale];
    gen_laguerre_ek_compute ( N_pts_Mob_Scale, alpha_Lag, x_L_Chi2, w_L_Chi2 );
    
    int num_pts_temp;
    long double ratio_E_orig, ratio_E, r_N1;
    long double norm_f, Chi2_Fit, Chi2_Fit_max;
    Chi2_Fit_max = 0.0;
    num_pts_temp = 100;
    
    for (int id_Mobius = 0; id_Mobius < N_pts_Mob_Scale; id_Mobius++) {
        Mobius_Scale_Actual = x_L_Chi2[id_Mobius];
        for (int incr_r_I0 = 0 + 1; incr_r_I0 < num_pts_temp - 1; incr_r_I0++) {
            ratio_E_orig = zeros_shifted(incr_r_I0, num_pts_temp, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
            for (int incr_N1 = 0; incr_N1 < num_pts_temp; incr_N1++) {
                norm_f = zeros_shifted(incr_N1 + num_pts_temp - 1, 2*(num_pts_temp - 1) + 1, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
                
                switch (Chi2_M1_1D_Current->Least_Squares_L_Chi2_Mode) {
                    case LEAST_SQUARES_L_CHI2_ON:
                        Mobius_Scale_Fit = Chi2_M1_1D_Current->Evaluate_Length_Scale_Least_Squares(norm_f);
                        break;
                    case LEAST_SQUARES_L_CHI2_OFF:
                        Mobius_Scale_Fit = Chi2_M1_1D_Current->Evaluate_Length_Scale(norm_f);
                        break;
                    default:
                        cout << "Invalid value for Least_Squares_L_Chi2_Mode !!!!!!!!!!!" << endl;
                        exit(0);
                        break;
                };
                
                ratio_E = Chi2_M1_1D_Current->Recompute_I0_Mapping(ratio_E_orig, Mobius_Scale_Actual, Mobius_Scale_Fit);
                
                Chi2_Fit = Chi2_M1_1D_Current->Evaluate_Chi2(ratio_E, norm_f);
                
                // cout << "ratio_E init = " << ratio_E << "   " << "ratio_E = " << ratio_E << "  " << "norm_f = " << norm_f << "Chi2_Fit = " << Chi2_Fit << endl;
                
                Chi2_Fit_max = max(Chi2_Fit_max, Chi2_Fit);
            }
        }
    }
    
    constr_val = -(1.0 - Chi2_Fit_max);
    constr_val *= 1.0e4;
            
    return constr_val;
 }
 
long double myconstraint_Nested_Least_Squares_L_Chi2_realizability_Lower_Bound(unsigned n, const long double *x, long double *grad, void *data) {
    long double constr_val;
    M1_1D_Data_Pointer *Chi2_M1_1D_Data_Pointer = (M1_1D_Data_Pointer *) data;
    
    M1_Data_1D_Cheby *Chi2_M1_1D_Current = Chi2_M1_1D_Data_Pointer->M1_Data_BE;
    M1_Data_1D_Cheby *Chi2_M1_1D_HL = Chi2_M1_1D_Data_Pointer->M1_Data_HL;
    M1_Data_1D_Cheby *Chi2_M1_1D_LL = Chi2_M1_1D_Data_Pointer->M1_Data_LL;
    M1_Data_1D_Cheby *Chi2_M1_1D_Uniform = Chi2_M1_1D_Data_Pointer->M1_Data_Uniform_BE;
    
//     for (int i = 0; i < Chi2_M1_1D_Current->N_pts_L_Chi2_LS; i++) {
//         Chi2_M1_1D_Current->Coefficient_Matrix_Fit_Mob_Scale[i] = x[i];
//     }
    
//     Chi2_M1_1D_Current->Polynomial_Interpolation_BE(*Chi2_M1_1D_HL, *Chi2_M1_1D_LL);
    
    // Compute the nodal distribution for the length scale L_Chi2
    double Mobius_Scale_Actual, Mobius_Scale_Fit;
    int N_pts_Mob_Scale = 100;
    long double alpha_Lag;
    alpha_Lag = 0;
    long double x_L_Chi2[N_pts_Mob_Scale], w_L_Chi2[N_pts_Mob_Scale];
    gen_laguerre_ek_compute ( N_pts_Mob_Scale, alpha_Lag, x_L_Chi2, w_L_Chi2 );
    
    int num_pts_temp;
    long double ratio_E_orig, ratio_E, r_N1;
    long double norm_f, norm_f_2, Chi2_Fit, Chi2_Fit_min;
    Chi2_Fit_min = 1.0;
    num_pts_temp = 100;
    norm_f_2 = norm_f * norm_f;
    
    for (int id_Mobius = 0; id_Mobius < N_pts_Mob_Scale; id_Mobius++) {
        Mobius_Scale_Actual = x_L_Chi2[id_Mobius];
        for (int incr_r_I0 = 0 + 1; incr_r_I0 < num_pts_temp - 1; incr_r_I0++) {
            ratio_E_orig = zeros_shifted(incr_r_I0, num_pts_temp, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
            for (int incr_N1 = 0; incr_N1 < num_pts_temp; incr_N1++) {
                norm_f = zeros_shifted(incr_N1 + num_pts_temp - 1, 2*(num_pts_temp - 1) + 1, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
                
                switch (Chi2_M1_1D_Current->Least_Squares_L_Chi2_Mode) {
                    case LEAST_SQUARES_L_CHI2_ON:
                        Mobius_Scale_Fit = Chi2_M1_1D_Current->Evaluate_Length_Scale_Least_Squares(norm_f);
                        break;
                    case LEAST_SQUARES_L_CHI2_OFF:
                        Mobius_Scale_Fit = Chi2_M1_1D_Current->Evaluate_Length_Scale(norm_f);
                        break;
                    default:
                        cout << "Invalid value for Least_Squares_L_Chi2_Mode !!!!!!!!!!!" << endl;
                        exit(0);
                        break;
                };
                
                ratio_E = Chi2_M1_1D_Current->Recompute_I0_Mapping(ratio_E_orig, Mobius_Scale_Actual, Mobius_Scale_Fit);
                
                Chi2_Fit = Chi2_M1_1D_Current->Evaluate_Chi2(ratio_E, norm_f);
                // cout << "ratio_E init = " << ratio_E << "   " << "ratio_E = " << ratio_E << "  " << "norm_f = " << norm_f << "Chi2_Fit = " << Chi2_Fit << endl;
                
                Chi2_Fit_min = min(Chi2_Fit_min, Chi2_Fit);
            }
        }
    }
    
    constr_val = -(Chi2_Fit_min - norm_f_2);
    constr_val *= 1.0e4;
            
    return constr_val;
 }
 
long double myconstraint_Nested_Least_Squares_L_Chi2_Hyperbolicity_Lambda_s(unsigned n, const long double *x, long double *grad, void *data) {
    long double constr_val;
    M1_1D_Data_Pointer *Chi2_M1_1D_Data_Pointer = (M1_1D_Data_Pointer *) data;
    
    M1_Data_1D_Cheby *Chi2_M1_1D_Current = Chi2_M1_1D_Data_Pointer->M1_Data_BE;
    M1_Data_1D_Cheby *Chi2_M1_1D_HL = Chi2_M1_1D_Data_Pointer->M1_Data_HL;
    M1_Data_1D_Cheby *Chi2_M1_1D_LL = Chi2_M1_1D_Data_Pointer->M1_Data_LL;
    M1_Data_1D_Cheby *Chi2_M1_1D_Uniform = Chi2_M1_1D_Data_Pointer->M1_Data_Uniform_BE;
    
//     for (int i = 0; i < Chi2_M1_1D_Current->N_pts_L_Chi2_LS; i++) {
//         Chi2_M1_1D_Current->Coefficient_Matrix_Fit_Mob_Scale[i] = x[i];
//     }
    
//     Chi2_M1_1D_Current->Polynomial_Interpolation_BE(*Chi2_M1_1D_HL, *Chi2_M1_1D_LL);
    
    // Compute the nodal distribution for the length scale L_Chi2
    double Mobius_Scale_Actual, Mobius_Scale_Fit;
    int N_pts_Mob_Scale = 100;
    long double alpha_Lag;
    alpha_Lag = 0;
    long double x_L_Chi2[N_pts_Mob_Scale], w_L_Chi2[N_pts_Mob_Scale];
    gen_laguerre_ek_compute ( N_pts_Mob_Scale, alpha_Lag, x_L_Chi2, w_L_Chi2 );
    
    long double ratio_E_orig, ratio_E;
    long double r_N1;
    long double norm_f;
    long double Chi2, dChi2_dN1, E_dChi2_dE;
    int num_pts_temp_E, num_pts_temp_f;
    
    long double temp_val, min_temp_val;
    min_temp_val = 1.0e12;
    num_pts_temp_E = 100;
    num_pts_temp_f = 100;
    
    long double E_dChi2_dE_min, dChi2_dN1_min, r_I0_min, norm_f_min, Chi2_min;
    
    for (int id_Mobius = 0; id_Mobius < N_pts_Mob_Scale; id_Mobius++) {
        Mobius_Scale_Actual = x_L_Chi2[id_Mobius];
        for (int incr_r_I0 = 0; incr_r_I0 < num_pts_temp_E; incr_r_I0++) {
            ratio_E_orig = zeros_shifted(incr_r_I0, num_pts_temp_E, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
            for (int incr_N1 = 0; incr_N1 < num_pts_temp_f; incr_N1++) {
                norm_f = zeros_shifted(incr_N1 + num_pts_temp_f - 1, 2*(num_pts_temp_f - 1) + 1, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
                
                switch (Chi2_M1_1D_Current->Least_Squares_L_Chi2_Mode) {
                    case LEAST_SQUARES_L_CHI2_ON:
                        Mobius_Scale_Fit = Chi2_M1_1D_Current->Evaluate_Length_Scale_Least_Squares(norm_f);
                        break;
                    case LEAST_SQUARES_L_CHI2_OFF:
                        Mobius_Scale_Fit = Chi2_M1_1D_Current->Evaluate_Length_Scale(norm_f);
                        break;
                    default:
                        cout << "Invalid value for Least_Squares_L_Chi2_Mode !!!!!!!!!!!" << endl;
                        exit(0);
                        break;
                };
                
                ratio_E = Chi2_M1_1D_Current->Recompute_I0_Mapping(ratio_E_orig, Mobius_Scale_Actual, Mobius_Scale_Fit);
                
                Chi2 = Chi2_M1_1D_Current->Evaluate_Chi2(ratio_E, norm_f);
                Chi2_M1_1D_Current->Evaluate_Chi2_derivatives(E_dChi2_dE, dChi2_dN1, ratio_E, norm_f);
                
                // cout << "ratio_E init = " << ratio_E << "   " << "ratio_E = " << ratio_E << "  " << "norm_f = " << norm_f << "Chi2_Fit = " << Chi2_Fit << endl;
                
                temp_val = pow(dChi2_dN1, 2) + 4.0*(Chi2 + E_dChi2_dE);
                
                min_temp_val = min(min_temp_val, temp_val);
                
                if (min_temp_val == temp_val) {
                    Chi2_min = Chi2;
                    E_dChi2_dE_min = E_dChi2_dE;
                    dChi2_dN1_min = dChi2_dN1;
                    r_I0_min = ratio_E;
                    norm_f_min = norm_f;
                }
            }
        }
    }
    
    constr_val = -min_temp_val;
    constr_val = 1.0e4*constr_val;
    
    if (Chi2_M1_1D_Current->id_proc == 0) {
        // for (int i = 0; i < Chi2_M1_1D_Current->N_Points_f; i++) {
        //     cout << "i = " << i << "   " << "x = " << x[i] << endl;
        // }
        cout << "min_temp_val = " << min_temp_val << "  " << "r_I0_min = " << r_I0_min << "  " << "norm_f_min = " << norm_f_min << "  " << "Chi2_min = " << Chi2_min << "  " << "E_dChi2_dE_min = " << E_dChi2_dE_min << "  " << "dChi2_dN1_min = " << dChi2_dN1_min << endl;
        
    }
            
    return constr_val;
}

long double myconstraint_Nested_Least_Squares_L_Chi2_Speed_Limit(unsigned n, const long double *x, long double *grad, void *data) {
    long double constr_val;
    M1_1D_Data_Pointer *Chi2_M1_1D_Data_Pointer = (M1_1D_Data_Pointer *) data;
    
    M1_Data_1D_Cheby *Chi2_M1_1D_Current = Chi2_M1_1D_Data_Pointer->M1_Data_BE;
    M1_Data_1D_Cheby *Chi2_M1_1D_HL = Chi2_M1_1D_Data_Pointer->M1_Data_HL;
    M1_Data_1D_Cheby *Chi2_M1_1D_LL = Chi2_M1_1D_Data_Pointer->M1_Data_LL;
    M1_Data_1D_Cheby *Chi2_M1_1D_Uniform = Chi2_M1_1D_Data_Pointer->M1_Data_Uniform_BE;
    
    long double ratio_E_orig, ratio_E;
    long double r_N1;
    long double norm_f;
    long double lambda_1, lambda_2, max_lambda;
    long double Chi2, dChi2_dN1, E_dChi2_dE;
    int num_pts_temp;
    
    long double temp_val, temp_val_temp;
    max_lambda = 0.0;
    num_pts_temp = 100;
    
    long double dChi2_dN1_max, r_I0_max, norm_f_max;
    
    for (int incr_r_I0 = 0 + 1; incr_r_I0 < num_pts_temp - 1; incr_r_I0++) {
        ratio_E_orig = zeros_shifted(incr_r_I0, num_pts_temp, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
        for (int incr_N1 = 0; incr_N1 < num_pts_temp; incr_N1++) {
            norm_f = zeros_shifted(incr_N1 + num_pts_temp - 1, 2*(num_pts_temp - 1) + 1, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
            Chi2 = Chi2_M1_1D_Current->Evaluate_Chi2(ratio_E, norm_f);
            // cout << "Chi2 = " << Chi2 << "  " << "norm_f = " << norm_f << "   " << "ratio_E = " << ratio_E << endl;
            Chi2_M1_1D_Current->Evaluate_Chi2_derivatives(E_dChi2_dE, dChi2_dN1, ratio_E, norm_f);
            
            temp_val_temp = pow(dChi2_dN1, 2) + 4.0*(Chi2 + E_dChi2_dE);
            temp_val = sqrt(temp_val_temp);
            
            lambda_1 = (1.0/2.0)*(dChi2_dN1 + temp_val);
            lambda_2 = (1.0/2.0)*(dChi2_dN1 - temp_val);
            
            max_lambda = max(max_lambda, fabs(lambda_1));
            max_lambda = max(max_lambda, fabs(lambda_2));
            
            if (max_lambda == fabs(lambda_1) || max_lambda == fabs(lambda_2)) {
                dChi2_dN1_max = dChi2_dN1;
                r_I0_max = ratio_E;
                norm_f_max = norm_f;
            }
            
        }
    }
    
    constr_val = -(1.0 - max_lambda);
    constr_val = 1.0e2*constr_val;
    
    if (Chi2_M1_1D_Current->id_proc == 0) {
        cout << "diff_lambdas_max = " << constr_val << "  " << "max_lambda = " << max_lambda << "  " << "r_I0_max = " << r_I0_max << "  " << "norm_f_max = " << norm_f_max << "  " << "dChi2_dN1_max = " << dChi2_dN1_max << "  " << "lambda_1 = " << lambda_1 << "  " << "lambda_2 = " << lambda_2 << endl;
    }
            
    return constr_val;
}

long double myconstraint_Nested_Least_Squares_L_Chi2_Hyperbolicity_Lambda_s_2D(unsigned n, const long double *x, long double *grad, void *data) {
    long double constr_val;
    M1_1D_Data_Pointer *Chi2_M1_1D_Data_Pointer = (M1_1D_Data_Pointer *) data;
    
    M1_Data_1D_Cheby *Chi2_M1_1D_Current = Chi2_M1_1D_Data_Pointer->M1_Data_BE;
    M1_Data_1D_Cheby *Chi2_M1_1D_HL = Chi2_M1_1D_Data_Pointer->M1_Data_HL;
    M1_Data_1D_Cheby *Chi2_M1_1D_LL = Chi2_M1_1D_Data_Pointer->M1_Data_LL;
    M1_Data_1D_Cheby *Chi2_M1_1D_Uniform = Chi2_M1_1D_Data_Pointer->M1_Data_Uniform_BE;
    
//     for (int i = 0; i < Chi2_M1_1D_Current->N_pts_L_Chi2_LS; i++) {
//         Chi2_M1_1D_Current->Coefficient_Matrix_Fit_Mob_Scale[i] = x[i];
//     }
    
//     Chi2_M1_1D_Current->Polynomial_Interpolation_BE(*Chi2_M1_1D_HL, *Chi2_M1_1D_LL);
    
    long double aa, bb, cc, dd, delta;
    long double dFdU_12, dFdU_21, dFdU_22, dFdU_23, dFdU_31, dFdU_32, dFdU_33;
    
    long double ratio_E_orig, ratio_E, norm_f, phi;
    long double Chi2, dChi2_dN1, E_dChi2_dE;
    int num_pts_temp, num_pts_temp_phi;
    
    long double min_det;
    min_det = 1.0e12;
    num_pts_temp = 20; //100;
    num_pts_temp_phi = 20;
    
    long double E_val;
    long double E_val_min, r_I0_min, norm_f_min, phi_min;
    long double Chi2_min, E_dChi2_dE_min, dChi2_dN1_min;
    
    // Chi2_M1_1D_Current->Problem_Type = GRAY;
    
    // Compute the nodal distribution for the length scale L_Chi2
    double Mobius_Scale_Actual, Mobius_Scale_Fit;
    int N_pts_Mob_Scale = 100;
    long double alpha_Lag;
    alpha_Lag = 0;
    long double x_L_Chi2[N_pts_Mob_Scale], w_L_Chi2[N_pts_Mob_Scale];
    gen_laguerre_ek_compute ( N_pts_Mob_Scale, alpha_Lag, x_L_Chi2, w_L_Chi2 );
    
    for (int id_Mobius = 0; id_Mobius < N_pts_Mob_Scale; id_Mobius++) {
        Mobius_Scale_Actual = x_L_Chi2[id_Mobius];
        for (int incr_r_I0 = 0; incr_r_I0 < num_pts_temp; incr_r_I0++) {
            ratio_E_orig = zeros_shifted(incr_r_I0, num_pts_temp, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
            // ratio_E_orig = 1.0;
            // cout << "incr_r_I0 = " << incr_r_I0 << "   " << "ratio_E_orig = " << ratio_E_orig << endl;
            for (int incr_N1 = 0; incr_N1 < num_pts_temp; incr_N1++) {
                norm_f = zeros_shifted(incr_N1 + num_pts_temp - 1, 2*(num_pts_temp - 1) + 1, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
                
                switch (Chi2_M1_1D_Current->Least_Squares_L_Chi2_Mode) {
                    case LEAST_SQUARES_L_CHI2_ON:
                        Mobius_Scale_Fit = Chi2_M1_1D_Current->Evaluate_Length_Scale_Least_Squares(norm_f);
                        break;
                    case LEAST_SQUARES_L_CHI2_OFF:
                        Mobius_Scale_Fit = Chi2_M1_1D_Current->Evaluate_Length_Scale(norm_f);
                        break;
                    default:
                        cout << "Invalid value for Least_Squares_L_Chi2_Mode !!!!!!!!!!!" << endl;
                        exit(0);
                        break;
                };
                
                E_val = Inverse_Mobius_Transformation(ratio_E, Mobius_Scale_Fit);
                
                ratio_E = Chi2_M1_1D_Current->Recompute_I0_Mapping(ratio_E_orig, Mobius_Scale_Actual, Mobius_Scale_Fit);
                
                Chi2 = Chi2_M1_1D_Current->Evaluate_Chi2(ratio_E, norm_f);
                Chi2_M1_1D_Current->Evaluate_Chi2_derivatives(E_dChi2_dE, dChi2_dN1, ratio_E, norm_f);
                
                // if (fabs(1.0 - ratio_E) > 1.0e-6)
                for (int incr_phi = 0; incr_phi < num_pts_temp_phi; incr_phi++) {
                    phi = Uniform_Distribution(incr_phi, num_pts_temp_phi, 0.0, 2.0*PI);
                    // phi = 0.0;
                    
                    dFdU_12 = 1.0;
                    dFdU_21 = dI2ij_dI0(1, 1, ratio_E, norm_f, phi, Chi2, E_dChi2_dE);
                    dFdU_22 = dI2ij_dIi(1, 1, 1, ratio_E, norm_f, phi, Chi2, dChi2_dN1);
                    dFdU_23 = dI2ij_dIi(1, 1, 2, ratio_E, norm_f, phi, Chi2, dChi2_dN1);
                    dFdU_31 = dI2ij_dI0(1, 2, ratio_E, norm_f, phi, Chi2, E_dChi2_dE);
                    dFdU_32 = dI2ij_dIi(1, 2, 1, ratio_E, norm_f, phi, Chi2, dChi2_dN1);
                    dFdU_33 = dI2ij_dIi(1, 2, 2, ratio_E, norm_f, phi, Chi2, dChi2_dN1);
                    
                    aa = -1.0;
                    bb = dFdU_22 + dFdU_33;
                    cc = dFdU_12*dFdU_21 - dFdU_22*dFdU_33 + dFdU_23*dFdU_32;
                    dd = -dFdU_12*dFdU_21*dFdU_33 + dFdU_12*dFdU_23*dFdU_31;
                    
                    // Discriminant of the general cubic function
                    delta = 18.0*aa*bb*cc*dd - 4.0*bb*bb*bb*dd + bb*bb*cc*cc - 4.0*aa*cc*cc*cc - 27.0*aa*aa*dd*dd;
                    
                    min_det = min(min_det, delta);
                    
                    if (min_det == delta) {
                        r_I0_min = ratio_E;
                        E_val_min = E_val;
                        norm_f_min = norm_f;
                        phi_min = phi;
                        
                        Chi2_min = Chi2;
                        E_dChi2_dE_min = E_dChi2_dE;
                        dChi2_dN1_min = dChi2_dN1;
                    }
                }
            }
        }
    }
    
    // The cubic equation has real roots if and only if delta >= 0 ==> min_det > 0
    // min_det > 0 ==> -min_det < 0
    constr_val = -min_det; // don't use pow(min_det, 1.0/3.0) in case min_det is negative;
    constr_val *= 1.0e4;
    
    if (Chi2_M1_1D_Current->id_proc == 0) {
        cout << "min_det in 2D = " << min_det << "  " << "E_val = " << E_val_min << "  " << "r_I0 = " << r_I0_min << "  " << "norm_f = " << norm_f_min << "  " << "phi = " << phi_min << "  " << "Chi2 = " << Chi2_min << "  " << "E_dChi2_dE = " << E_dChi2_dE_min << "  " << "dChi2_dN1= " << dChi2_dN1_min << endl;
    }
            
    // Chi2_M1_1D_Current->Problem_Type = NON_GRAY;
    
    return constr_val;
}


void M1_Data_1D_Cheby :: Compute_dI2ij_dIn_Finite_Difference(const long double &r_I0, const long double &N1, const long double &phi) {
    long double dFdU_11, dFdU_12, dFdU_13, dFdU_21, dFdU_22, dFdU_23, dFdU_31, dFdU_32, dFdU_33;
    long double E_dChi2_dE, dChi2_df;
    long double epsilon = 1.0e-6;
    long double Mobius_Scale_Fit, E;
    
    switch (Problem_Type) {
        case GRAY:
            dChi2_df = Evaluate_dChi2_dN1(r_I0, N1);
            E_dChi2_dE = -dChi2_df*N1;
            break;
        case NON_GRAY:
            switch (Least_Squares_L_Chi2_Mode) {
                case LEAST_SQUARES_L_CHI2_ON:
                    Mobius_Scale_Fit = Evaluate_Length_Scale_Least_Squares(N1);
                    break;
                case LEAST_SQUARES_L_CHI2_OFF:
                    Mobius_Scale_Fit = Evaluate_Length_Scale(N1);
                    break;
            };
            E = Inverse_Mobius_Transformation(r_I0, Mobius_Scale_Fit);
            break;
        default:
            cout << "Problem Type not specified" << endl;
            exit(0);
            break;
    }
    
    long double Chi2, Chi2_L, Chi2_R;
    long double I2_11_L, I2_11_R, I2_12_L, I2_12_R;
    long double dI2_11_dI0, dI2_12_dI0, dI2_11_dI1_1, dI2_12_dI1_1, dI2_11_dI1_2, dI2_12_dI1_2;
    
    // Finite differencing with respect to I0
    I2_11_L = I2ij(1, 1, r_I0, N1, phi, Chi2);
    I2_11_R = I2ij(1, 1, r_I0, N1, phi, Chi2);
    dI2_11_dI0 = (I2_11_R - I2_11_L)/(2.0*epsilon);
    
    I2_12_L = I2ij(1, 2, r_I0, N1, phi, Chi2);
    I2_12_R = I2ij(1, 2, r_I0, N1, phi, Chi2);
    dI2_12_dI0 = (I2_12_R - I2_12_L)/(2.0*epsilon);
    
    // Finite differencing with respect to I1_1
    I2_11_L = I2ij(1, 1, r_I0, N1, phi, Chi2);
    I2_11_R = I2ij(1, 1, r_I0, N1, phi, Chi2);
    dI2_11_dI1_1 = (I2_11_R - I2_11_L)/(2.0*epsilon);
    
    I2_12_L = I2ij(1, 2, r_I0, N1, phi, Chi2);
    I2_12_R = I2ij(1, 2, r_I0, N1, phi, Chi2);
    dI2_12_dI1_1 = (I2_12_R - I2_12_L)/(2.0*epsilon);
    
    // Finite differencing with respect to I1_2
    I2_11_L = I2ij(1, 1, r_I0, N1, phi, Chi2);
    I2_11_R = I2ij(1, 1, r_I0, N1, phi, Chi2);
    dI2_11_dI1_2 = (I2_11_R - I2_11_L)/(2.0*epsilon);
    
    I2_11_L = I2ij(1, 1, r_I0, N1, phi, Chi2);
    I2_11_R = I2ij(1, 1, r_I0, N1, phi, Chi2);
    dI2_12_dI1_2 = (I2_11_R - I2_11_L)/(2.0*epsilon);
    
}

// ********************************************************************************
// This routine computes the L^{1} and L^{2} errors of the polynomial interpolation
// for the least squares algorithm for computing optimal values of L_Chi2
// The polynomial interpolation consists of performing interpolation with respect
// to both the exponential mapping for the radiative energy density and the first-
// order normalized angular moment
// ********************************************************************************
void M1_Data_1D_Cheby :: Compute_L_ONE_L_TWO_Errors_Nested_Least_Squares_L_Chi2(M1_Data_1D_Cheby *Chi2_M1_1D_Uniform, const int &flag_Write_Output, ofstream &out_L_inf_Norm, const int &VAR_NUM) {
    long double max_err_Chi2_E, max_err_Chi2_N1_1;
    long double L2_Norm_Chi2, L_inf_Norm_Chi2;
    long double Chi2_Fit, Chi2_Numerical;
    long double ratio_E, norm_f;
    long double Mobius_Scale_Actual, Mobius_Scale_Fit;
    long double error_fit;
    long double weight_total;
    long double weight_L_Chi2_val, weight_r_I0_val, weight_N1_val;
    long double E_dChi2_dE, dChi2_dN1;
    int index;
    L2_Norm_Chi2 = 0.0;
    L_inf_Norm_Chi2 = 0.0;
    
    long double x_L_Chi2[Chi2_M1_1D_Uniform->N_pts_Mob_Scale], weight_L_Chi2[Chi2_M1_1D_Uniform->N_pts_Mob_Scale];
    long double x_r_I0[Chi2_M1_1D_Uniform->N_Points_E], weight_r_I0[Chi2_M1_1D_Uniform->N_Points_E];
    long double x_N1[2*(Chi2_M1_1D_Uniform->N_Points_f - 1) + 1], weight_N1[2*(Chi2_M1_1D_Uniform->N_Points_f - 1) + 1];
    
    // chebyshev_quadrature ( weight_r_I0, x_r_I0, Chi2_M1_1D_Uniform->N_Points_E, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
    // chebyshev_quadrature ( weight_N1, x_N1, 2*(Chi2_M1_1D_Uniform->N_Points_f - 1) + 1, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
    
    lobatto_compute ( Chi2_M1_1D_Uniform->N_Points_E, x_r_I0, weight_r_I0 );
    lobatto_compute ( 2*(Chi2_M1_1D_Uniform->N_Points_f - 1) + 1, x_N1, weight_N1 );
    
    gen_laguerre_ek_compute ( Chi2_M1_1D_Uniform->N_pts_Mob_Scale, 0.0, x_L_Chi2, weight_L_Chi2 );
    
    for (int id_Mobius = 0; id_Mobius < Chi2_M1_1D_Uniform->N_pts_Mob_Scale; id_Mobius++) {
        Mobius_Scale_Actual = x_L_Chi2[id_Mobius];
        for (int i_E = 0; i_E < Chi2_M1_1D_Uniform->N_Points_E; i_E++) {
            for (int i_f = 0+1; i_f < Chi2_M1_1D_Uniform->N_Points_f-1; i_f++) {
                index = (id_Mobius*Chi2_M1_1D_Uniform->N_Points_E + i_E)*Chi2_M1_1D_Uniform->N_Points_f + i_f;
                
                ratio_E = Chi2_M1_1D_Uniform->ratio_I0_NON_GRAY[index];
                norm_f = Chi2_M1_1D_Uniform->N1_1_NON_GRAY[index];
                
                switch (Least_Squares_L_Chi2_Mode) {
                    case LEAST_SQUARES_L_CHI2_ON:
                        Mobius_Scale_Fit = Evaluate_Length_Scale_Least_Squares(norm_f);
                        break;
                    case LEAST_SQUARES_L_CHI2_OFF:
                        Mobius_Scale_Fit = Evaluate_Length_Scale(norm_f);
                        break;
                    default:
                        cout << "Invalid value for Least_Squares_L_Chi2_Mode !!!!!!!!!!!" << endl;
                        exit(0);
                        break;
                };
                
                ratio_E = Recompute_I0_Mapping(ratio_E, Mobius_Scale_Actual, Mobius_Scale_Fit);
                
                // Chi2_Numerical = Chi2_M1_1D_Uniform->Chi2_NON_GRAY[index];
                // Chi2_Fit = Evaluate_Chi2(ratio_E, norm_f);
                
                Chi2_Numerical = Chi2_M1_1D_Uniform->f_Chi2_NON_GRAY[index];
                Chi2_Fit = Evaluate_h_Chi2(ratio_E, norm_f);
                
//                 Chi2_Numerical = Chi2_M1_1D_Uniform->dChi2_dI0[index];
                
//                 Evaluate_Chi2_derivatives(E_dChi2_dE, dChi2_dN1, ratio_E, norm_f);
//                 Chi2_Fit = Evaluate_Chi2(ratio_E, norm_f);
//                 Chi2_Fit = Chi2_Fit + E_dChi2_dE;
                
                // cout << "ratio_E = " << ratio_E << "   " << "norm_f = " << norm_f << "   " << "Chi2_Fit = " << Chi2_Fit << "    " << "Chi2_Numerical = " << Chi2_Numerical << endl;
                
                error_fit = Chi2_Fit - Chi2_Numerical;
                
//                 if (x_r_I0[i_E] != Chi2_M1_1D_Uniform->ratio_I0_NON_GRAY[index] && Chi2_M1_1D_Uniform->E_NON_GRAY[index] < 1.0e7) {
//                     cout << "Problem with gauss-Lobatto-Chebyshev quadrature for r_I0" << endl;
//                     cout << "i_E = " << i_E << "  " << "i_f = " << i_f << "  " << "x_r_I0 = " << x_r_I0[i_E] << "  " << "ratio_I0_NON_GRAY = " << Chi2_M1_1D_Uniform->ratio_I0_NON_GRAY[index] << endl;
//                     exit(0);
//                 }
//                 
//                 if (x_N1[i_f + (Chi2_M1_1D_Uniform->N_Points_f - 1)] != norm_f) {
//                     cout << "Problem with gauss-Lobatto-Chebyshev quadrature for N1" << endl;
//                     cout << "i_E = " << i_E << "  " << "i_f = " << i_f << "  " << "x_N1 = " << x_N1[i_f + (Chi2_M1_1D_Uniform->N_Points_f - 1)] << "  " << "norm_f = " << norm_f << endl;
//                     exit(0);
//                 }
                
                // cout << "i_f = " << i_f << "  " << "x_N1 = " << x_N1[i_f + (Chi2_M1_1D_Uniform->N_Points_f - 1)] << "  " << "norm_f = " << norm_f << endl;
                
                weight_L_Chi2_val = exp(x_L_Chi2[id_Mobius]) * weight_L_Chi2[id_Mobius];
                
//                 if (i_E == 0 && i_f == 0) { 
//                     cout << "weight_L_Chi2_val = " << weight_L_Chi2_val << endl;
//                 }
                
                if (isinf(weight_L_Chi2_val) || isnan(weight_L_Chi2_val)) {
                    weight_L_Chi2_val = 0.0;
                    cout << "weight_L_Chi2_val = " << weight_L_Chi2_val << endl;
                    cout << "id_Mobius = " << id_Mobius << "  " << "x_L_Chi2 = " << x_L_Chi2[id_Mobius] << "  " << "weight_L_Chi2 = " << weight_L_Chi2[id_Mobius] << endl;
                    exit(0);
                }
                
                weight_r_I0_val = weight_r_I0[i_E];
                weight_N1_val = weight_N1[i_f + (Chi2_M1_1D_Uniform->N_Points_f - 1)];
                if (i_f != 0) {
                    weight_N1_val *= 2.0;
                }
                
                weight_total = weight_L_Chi2_val * weight_r_I0_val * weight_N1_val;
                
                L2_Norm_Chi2 += /*weight_total **/ pow(error_fit, 2);
                L_inf_Norm_Chi2 = max(L_inf_Norm_Chi2, fabs(error_fit));
                
                if (L_inf_Norm_Chi2 == fabs(error_fit)) {
                    max_err_Chi2_E = Chi2_M1_1D_Uniform->ratio_I0_NON_GRAY[index];
                    max_err_Chi2_N1_1 = Chi2_M1_1D_Uniform->N1_1_NON_GRAY[index];
                } // end if statement
            } // end for i_f
        } // end for i_E
    } // end for id_Mobius
    
    L2_Norm_Chi2 = sqrt(L2_Norm_Chi2);
    
    Chi2_M1_1D_Uniform->L_inf_Norm = L_inf_Norm_Chi2;
    
    Chi2_M1_1D_Uniform->L2_Norm_Chi2 = L2_Norm_Chi2;
    
    if (flag_Write_Output) {
        out_L_inf_Norm  << L_inf_Norm_Chi2 << setw(16) << L2_Norm_Chi2 << endl;
    }
    
    if (id_proc == 0) {
        cout << "Convergence Stats ........." << "L2_Norm_Chi2 = " << L2_Norm_Chi2 << "     "  << "L_inf_Norm_Chi2 = " << L_inf_Norm_Chi2 << "    " << "max_err_Chi2_E = " << max_err_Chi2_E << "    " << "max_err_Chi2_N1_1 = " << max_err_Chi2_N1_1 << endl;
    }
}

long double M1_Data_1D_Cheby :: Evaluate_diff_Length_Scale_dN1_Least_Squares(const long double &norm_f) {
    long double dLength_Scale, Length_Scale;
    long double norm_f_2 = norm_f*norm_f;
    int index = N_pts_L_Chi2_LS - 1;
    
    long double dpoly_val;
    Length_Scale = Evaluate_Length_Scale_Least_Squares(norm_f);
    dLength_Scale = 0.0;
    for (int i_fit_f = 0; i_fit_f < N_pts_L_Chi2_LS; i_fit_f++) {
        dpoly_val = 2*i_fit_f*Chebyshev_Second_Kind_Polynomial_Basis(norm_f, 2*i_fit_f-1);
        dLength_Scale += Coefficient_Matrix_Fit_Mob_Scale[i_fit_f] * dpoly_val;
    }
    dLength_Scale *= Length_Scale;
    
    return dLength_Scale;
}

long double M1_Data_1D_Cheby :: Evaluate_Length_Scale_Least_Squares(const long double &norm_f) {
    long double Length_Scale;
    int index = N_pts_L_Chi2_LS - 1;
    
    long double poly_val;
    Length_Scale = 0.0;
    for (int i_fit_f = 0; i_fit_f < N_pts_L_Chi2_LS; i_fit_f++) {
        poly_val = Chebyshev_Polynomial_Basis(norm_f, 2*i_fit_f);
        Length_Scale += Coefficient_Matrix_Fit_Mob_Scale[i_fit_f] * poly_val;
    }
    
    Length_Scale = exp(Length_Scale);
    
    // cout << "norm_f = " << norm_f << "  " << "Length_Scale = " << Length_Scale << endl;
    
    return Length_Scale;
}
