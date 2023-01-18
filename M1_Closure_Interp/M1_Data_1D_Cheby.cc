#ifndef _M1_DATA_1D_CHEBY_H_INCLUDED
#include "M1_Data_1D_Cheby.h"
#endif // _M1_DATA_1D_CHEBY_H_INCLUDED

int M1_Data_1D_Cheby :: Problem_Type = 0;
bool M1_Data_1D_Cheby :: flag_test_Implementation = false;
int M1_Data_1D_Cheby :: Implementation_type = 0;
int M1_Data_1D_Cheby :: Least_Squares_L_Chi2_Mode = LEAST_SQUARES_L_CHI2_OFF;
int M1_Data_1D_Cheby :: Weighting_Function_fchi2_Interpolation_Type = POLYNOMIAL_TYPE; //RATIONAL_TYPE; //

void M1_Data_1D_Cheby :: allocate() {
    // Deallocate first
    deallocate();
    
    // Allocate now
    int N_pts_total = N_Points_E*N_Points_f;
    
    Coefficients_Vander_Matrix_2vars = new long double[N_pts_total*N_pts_total];
    
    Coefficient_Matrix_Fit_Chi2 = new long double [N_pts_total];
    Coefficient_Matrix_Fit_Mob_Scale = new long double [N_pts_L_Chi2_LS];
    
    Coefficient_Matrix_Fit_Chi2_HL = new long double [N_Points_f];
    Coefficient_Matrix_Fit_Chi2_LL = new long double [N_Points_f];
    
    N_pts_total = N_Points_E*N_Points_f*N_pts_Mob_Scale;
    
    E_NON_GRAY = new long double [N_pts_total];
    ratio_I0_NON_GRAY = new long double [N_pts_total];
    N1_1_NON_GRAY = new long double [N_pts_total];
        
    Chi2_NON_GRAY = new long double [N_pts_total];
        
    f_Chi2_NON_GRAY = new long double [N_pts_total];
    error_fit = new long double [N_pts_total];
    
    dChi2_drI0 = new long double [N_pts_total];
    dChi2_dI0 = new long double [N_pts_total];
    dChi2_dN1 = new long double [N_pts_total];
    d2Chi2_dN12 = new long double [N_pts_total];
    d3Chi2_dN13 = new long double [N_pts_total];
    d2Chi2_drI0_dN1 = new long double [N_pts_total];
    d3Chi2_drI0_dN1_2 = new long double [N_pts_total];
    d2Chi2_dI0_dN1 = new long double [N_pts_total];
    d3Chi2_dI0_dN1_2 = new long double [N_pts_total];
}

void M1_Data_1D_Cheby :: deallocate() {
    if (Coefficients_Vander_Matrix_2vars != NULL) {
        delete[] Coefficients_Vander_Matrix_2vars; Coefficients_Vander_Matrix_2vars = NULL;
    }
    
    if (E_NON_GRAY != NULL) {
        delete[] E_NON_GRAY; E_NON_GRAY = NULL;
    }
    
    if (N1_1_NON_GRAY != NULL) {
        delete[] N1_1_NON_GRAY; N1_1_NON_GRAY = NULL;
    }
    
    if (dChi2_drI0 != NULL) {
        delete[] dChi2_drI0; dChi2_drI0 = NULL;
    }
    
    if (dChi2_dI0 != NULL) {
        delete[] dChi2_dI0; dChi2_dI0 = NULL;
    }
    
    if (dChi2_dN1 != NULL) {
        delete[] dChi2_dN1; dChi2_dN1 = NULL;
    }
    
    if (d2Chi2_dN12 != NULL) {
        delete[] d2Chi2_dN12; d2Chi2_dN12 = NULL;
    }
    
    if (d3Chi2_dN13 != NULL) {
        delete[] d3Chi2_dN13; d3Chi2_dN13 = NULL;
    }
    
    if (d2Chi2_drI0_dN1 != NULL) {
        delete[] d2Chi2_drI0_dN1; d2Chi2_drI0_dN1 = NULL;
    }
    
    if (d3Chi2_drI0_dN1_2 != NULL) {
        delete[] d3Chi2_drI0_dN1_2; d3Chi2_drI0_dN1_2 = NULL;
    }
    
    if (d2Chi2_dI0_dN1 != NULL) {
        delete[] d2Chi2_dI0_dN1; d2Chi2_dI0_dN1 = NULL;
    }
    
    if (d3Chi2_dI0_dN1_2 != NULL) {
        delete[] d3Chi2_dI0_dN1_2; d3Chi2_dI0_dN1_2 = NULL;
    }
    
    if (Chi2_NON_GRAY != NULL) {
        delete[] Chi2_NON_GRAY; Chi2_NON_GRAY = NULL;
    }
    
    if (f_Chi2_NON_GRAY != NULL) {
        delete[] f_Chi2_NON_GRAY; f_Chi2_NON_GRAY = NULL;
    }
    
    if (error_fit != NULL) {
        delete[] error_fit; error_fit = NULL;
    }
    
    if (Coefficient_Matrix_Fit_Chi2 != NULL) {
        delete[] Coefficient_Matrix_Fit_Chi2; Coefficient_Matrix_Fit_Chi2 = NULL;
    }
    
    if (Coefficient_Matrix_Fit_Mob_Scale != NULL) {
        delete[] Coefficient_Matrix_Fit_Mob_Scale; Coefficient_Matrix_Fit_Mob_Scale = NULL;
    }
}

void M1_Data_1D_Cheby :: Allocate_Coefficients_L_Chi2() {
    Deallocate_Coefficients_L_Chi2();
    Coefficient_Matrix_Fit_Mob_Scale = new long double [N_pts_L_Chi2_LS];
}

void M1_Data_1D_Cheby :: Deallocate_Coefficients_L_Chi2() {
    if (Coefficient_Matrix_Fit_Mob_Scale != NULL) {
        delete[] Coefficient_Matrix_Fit_Mob_Scale; Coefficient_Matrix_Fit_Mob_Scale = NULL;
    }
}

void M1_Data_1D_Cheby :: OpenInputFile(char *filename) {
    strcpy(path_out, getenv(PATHVAR));
    strcat(path_out, "/M1_Model/Non_Gray_M1_Closure/");
    strcat(path_out, filename);
    sprintf(extension, "_%.6d", 0);
    strcat(extension, ".dat");
    strcat(path_out, extension);
        
//     fstream in_out;
    in_out_Optim.open(path_out, ios::in|ios::binary);
    
    if (!in_out_Optim) {
        cout << path_out << "could not be accessed!" << endl;    
    }
}

void M1_Data_1D_Cheby :: CloseInputFile() {
    in_out_Optim.close();
}

void M1_Data_1D_Cheby :: ReadInputData() {
    int index;
    if (in_out_Optim.good()) {
        for (int id_Mob_Scale = 0 ; id_Mob_Scale < N_pts_Mob_Scale; id_Mob_Scale++) {
            for (int index_e = 0 ; index_e < N_Points_E; index_e++) {
                for (int i_f = 0; i_f < N_Points_f; i_f++) {
                    index = (id_Mob_Scale*N_Points_E + index_e)*N_Points_f + i_f;
                    
                    read<record_Chi2>(in_out_Optim, rec_Chi2, index);
                    
                    ratio_I0_NON_GRAY[index] = rec_Chi2.ratio_I0;
                    E_NON_GRAY[index] = rec_Chi2.I0;
                    N1_1_NON_GRAY[index] = rec_Chi2.N1;
                    Chi2_NON_GRAY[index] = rec_Chi2.Chi2;
                    
                    dChi2_drI0[index] = rec_Chi2.dChi2_drI0;
                    dChi2_dI0[index] = rec_Chi2.dChi2_dI0;
                    dChi2_dN1[index] = rec_Chi2.dChi2_dN1;
                    d2Chi2_dN12[index] = rec_Chi2.d2Chi2_dN12;
                    d3Chi2_dN13[index] = rec_Chi2.d3Chi2_dN13;
                    d2Chi2_drI0_dN1[index] = rec_Chi2.d2_Chi2_drI0_dN1;
                    d3Chi2_drI0_dN1_2[index] = rec_Chi2.d3Chi2_drI0_dN1_2;
                    d2Chi2_dI0_dN1[index] = rec_Chi2.d2_Chi2_dI0_dN1;
                    d3Chi2_dI0_dN1_2[index] = rec_Chi2.d3Chi2_dI0_dN1_2;
                    
                    // cout << "index = " << index << "  " << "E_NON_GRAY = " << E_NON_GRAY[index] << "     " << "N1_1_NON_GRAY = " << N1_1_NON_GRAY[index] << "     " << "Chi2_NON_GRAY = " << Chi2_NON_GRAY[index] << endl;
                    
                    if (Chi2_NON_GRAY[index] < pow(N1_1_NON_GRAY[index], 2)) {
                        cout << "Chi2_NON_GRAY = " << Chi2_NON_GRAY[index] << "     " << "N1_1_NON_GRAY = " << N1_1_NON_GRAY[index] << "     " << "pow(N1_1_NON_GRAY[index], 2) = " << pow(N1_1_NON_GRAY[index], 2) << endl;
                        cout << "Non Realizable moment" << endl;
                        exit(0);
                    }
                }
            }
        }
        
//         cout << "Checking hyperbolicity for Max Entropy solutions !!!!!!!!!!!!!!!" << endl;
//         Check_Hyperbolicity_Max_Ent_Solutions(0);
    }
}

void M1_Data_1D_Cheby :: Check_Hyperbolicity_Max_Ent_Solutions(const int &Type) {
    int index;
    long double r_I0_min, I0_min, norm_f_min, phi_min, Chi2_min, E_dChi2_dE_min, dChi2_dnorm_f_min;
    long double ratio_E, I0, norm_f, Chi2, phi;
    long double E_dChi2_dE, dChi2_dnorm_f;
    long double delta_min, delta, aa, bb, cc, dd;
    long double dFdU_12, dFdU_21, dFdU_22, dFdU_23, dFdU_31, dFdU_32, dFdU_33;
    int num_pts_temp_phi = 20;
    delta_min = 1.0e12;
    for (int id_Mob_Scale = 0 ; id_Mob_Scale < N_pts_Mob_Scale; id_Mob_Scale++) {
        for (int index_e = 0 ; index_e < N_Points_E; index_e++) {
            for (int i_f = 0; i_f < N_Points_f; i_f++) {
                index = (id_Mob_Scale*N_Points_E + index_e)*N_Points_f + i_f;
                
                I0 = E_NON_GRAY[index];
                ratio_E = ratio_I0_NON_GRAY[index];
                norm_f = N1_1_NON_GRAY[index];
                Chi2 = Chi2_NON_GRAY[index];
                
                switch (Type) {
                    case 0:
                        E_dChi2_dE = dChi2_dI0[index] - Chi2;
                        dChi2_dnorm_f = dChi2_dN1[index];
                        break;
                    case 1:
                        Evaluate_Chi2_derivatives(E_dChi2_dE, dChi2_dnorm_f, ratio_E, norm_f);
                        break;
                    default:
                        cout << "Invalid value for Type !!!!!!!!!!!!!!!!!!" << endl;
                        exit(0);
                        break;
                }
                
                for (int incr_phi = 0; incr_phi < num_pts_temp_phi; incr_phi++) {
                    phi = Uniform_Distribution(incr_phi, num_pts_temp_phi, 0.0, 2.0*PI);
                    // phi = 0.0;
                    
                    dFdU_12 = 1.0;
                    dFdU_21 = dI2ij_dI0(1, 1, ratio_E, norm_f, phi, Chi2, E_dChi2_dE);
                    dFdU_22 = dI2ij_dIi(1, 1, 1, ratio_E, norm_f, phi, Chi2, dChi2_dnorm_f);
                    dFdU_23 = dI2ij_dIi(1, 1, 2, ratio_E, norm_f, phi, Chi2, dChi2_dnorm_f);
                    dFdU_31 = dI2ij_dI0(1, 2, ratio_E, norm_f, phi, Chi2, E_dChi2_dE);
                    dFdU_32 = dI2ij_dIi(1, 2, 1, ratio_E, norm_f, phi, Chi2, dChi2_dnorm_f);
                    dFdU_33 = dI2ij_dIi(1, 2, 2, ratio_E, norm_f, phi, Chi2, dChi2_dnorm_f);
                    
                    aa = -1.0;
                    bb = dFdU_22 + dFdU_33;
                    cc = dFdU_12*dFdU_21 - dFdU_22*dFdU_33 + dFdU_23*dFdU_32;
                    dd = -dFdU_12*dFdU_21*dFdU_33 + dFdU_12*dFdU_23*dFdU_31;
                    
                    // Discriminant of the general cubic function
                    delta = 18.0*aa*bb*cc*dd - 4.0*bb*bb*bb*dd + bb*bb*cc*cc - 4.0*aa*cc*cc*cc - 27.0*aa*aa*dd*dd;
                    
                    delta_min = min(delta_min, delta);
                    
                    if (delta_min == delta) {
                        r_I0_min = ratio_E;
                        I0_min = I0;
                        norm_f_min = norm_f;
                        phi_min = phi;
                        Chi2_min = Chi2;
                        E_dChi2_dE_min = E_dChi2_dE;
                        dChi2_dnorm_f_min = dChi2_dnorm_f;
                    }
                }
            }
        }
    }
    
    if (delta_min < 0.0) {
        cout << "Problem with hyperbolicity, delta_min = " << delta_min << endl;
        
        cout << "delta_min = " << delta_min << "  " << "I0 = " << I0_min << "  " << "norm_f = " << norm_f_min << "  " << "phi = " << phi_min << "  " << "Chi2 = " << Chi2_min << "  " << "E_dChi2_dE = " << E_dChi2_dE_min << "  " << "dChi2_dnorm_f = " << dChi2_dnorm_f_min << endl;
        exit(0);
    }
}

//***********************************************************************************************
// This routine sets up interpolant values in either the Hyperbolic or the Logarithmic limit
//***********************************************************************************************
void M1_Data_1D_Cheby :: SetupInterpolant_Values_HL_LL() {
    int index;
    long double norm_f_2;
    
    if (N_Points_E != 1) {
        cout << "Problem with number of points in HL or LL" << endl;
        exit(0);
    }
    
    for (int i_Cheby_f = 0; i_Cheby_f < N_Points_f; i_Cheby_f++) {
        index = i_Cheby_f;
        
        norm_f_2 = pow(N1_1_NON_GRAY[index],2);
        
        f_Chi2_NON_GRAY[index] = (3.0*Chi2_NON_GRAY[index] - 1.0)/2.0;
        
        f_Chi2_NON_GRAY[index] = (f_Chi2_NON_GRAY[index] - norm_f_2)/(norm_f_2*(1.0 - norm_f_2));
        
        if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
            // Apply l'Hopital's rule twice since the first derivative is also zero
            // in the isotropic regime
            f_Chi2_NON_GRAY[index] = 3.0*d2Chi2_dN12[index]/2.0;
            f_Chi2_NON_GRAY[index] = f_Chi2_NON_GRAY[index]/2.0 - 1.0;
            f_Chi2_NON_GRAY[index] /= 1.0 - norm_f_2;
        } else if (fabs(N1_1_NON_GRAY[index] - 1.0) < 1.0e-8) {
            // Then free-streaming limit
            f_Chi2_NON_GRAY[index] = 3.0*dChi2_dN1[index]/2.0;
            f_Chi2_NON_GRAY[index] = 1.0 - f_Chi2_NON_GRAY[index]/2.0;
        }
        
        if (f_Chi2_NON_GRAY[index] != f_Chi2_NON_GRAY[index]) {
            cout << "In Hyperbolic or Logarithmic limit" << endl;
            cout << "E_NON_GRAY = " << E_NON_GRAY[index] << "     " << "N1_1_NON_GRAY = " << N1_1_NON_GRAY[index] << "     " << "Chi2_NON_GRAY = " << Chi2_NON_GRAY[index] << "     " << "f_Chi2_NON_GRAY = " << f_Chi2_NON_GRAY[index] << "     " << "dChi2_dN1 = " << dChi2_dN1[index] << endl;
            exit(0);
        }
    } // end for i_Cheby_f
} 


void M1_Data_1D_Cheby :: SetupInterpolant_Values_BE(M1_Data_1D_Cheby &M1_Data_Chi2_HL, M1_Data_1D_Cheby &M1_Data_Chi2_LL) {
    switch (Implementation_type) {
        case IMPLEMENTATION_DERIVATIVES_N1:
            SetupInterpolant_Values_BE_Impl_Derivs_N1(M1_Data_Chi2_HL, M1_Data_Chi2_LL);
            break;
        case IMPLEMENTATION_DERIVATIVES_N1_DERIVATIVES_R_I0:
            SetupInterpolant_Values_BE_Impl_Derivs_N1_Derivs_r_I0(M1_Data_Chi2_HL, M1_Data_Chi2_LL);
            break;
        default:
            cout << "Implementation type not specified, value is " << Implementation_type << endl;
            exit(0);
            break;
    }
}

void M1_Data_1D_Cheby :: SetupInterpolant_Values_BE_Impl_Derivs_N1(M1_Data_1D_Cheby &M1_Data_Chi2_HL, M1_Data_1D_Cheby &M1_Data_Chi2_LL) {
    int index;
    long double norm_f_2, norm_f_10;
    
    long double diff_f_Chi2_NON_GRAY, g_Chi2_HL, g_Chi2_LL;
    
    // In the case of the optimization problem for the M1 closure: I0 and I1
    // are the independent variables, which allows us to write Chi2 = Chi2(I0, I1). 
    // Defining N1 = I1/I0, we can write Chi2 = Chi2(I0, N1) where N1 = N1(I0, I1)
    // In the case where we consider the space (I0, N1), i.e., assume I0 and N1 to 
    // be independent, we have
    // {d Chi2}{d N1} = {\partial Chi2}{\partial N1}
    // In the space (r_I0, N1) where r_I0 = r_I0(I0, N1), as is the case for our interpolation
    // of the Eddington factor for the M1 closure, we have
    // {d Chi2}{d N1} = {\partial Chi2}{\partial N1} + {\partial Chi2}{\partial r_I0}.{\partial r_I0}{\partial N1}
    // We can then write in the space (r_I0, N1) where r_I0 = r_I0(I0, N1):
    // {\partial Chi2}{\partial N1} = {d Chi2}{d N1} - {\partial Chi2}{\partial r_I0}.{\partial r_I0}{\partial N1}
    // This is how we compute the partial first derivatives of Chi2 in the space (r_I0, N1) where r_I0 = r_I0(I0, N1) from the space (I0, N1).
    // We can also write:
    // {d^2 Chi2}{d N1^2} = {\partial^2 Chi2}{\partial N1^2} + {\partial Chi2}{\partial r_I0}.{\partial r_I0}{\partial N1^2}
    
    for (int id_Mobius = 0; id_Mobius < N_pts_Mob_Scale; id_Mobius++) {
        for (int i_Cheby_E = 0; i_Cheby_E < N_Points_E; i_Cheby_E++) {
            for (int i_Cheby_f = 0; i_Cheby_f < N_Points_f; i_Cheby_f++) {
                index = (id_Mobius * N_Points_E + i_Cheby_E)*N_Points_f + i_Cheby_f;
                
                norm_f_2 = pow(N1_1_NON_GRAY[index],2);
                
                norm_f_10 = pow(norm_f_2, 1);
                
                switch (Weighting_Function_fchi2_Interpolation_Type) {
                    case POLYNOMIAL_TYPE:
                        f_Chi2_NON_GRAY[index] = (3.0*Chi2_NON_GRAY[index] - 1.0)/2.0;
                        f_Chi2_NON_GRAY[index] = (f_Chi2_NON_GRAY[index] - norm_f_2)/(norm_f_2*(1.0 - norm_f_10));
                        
                        if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
                            // Apply l'Hopital's rule twice since the first derivative is also zero
                            // in the isotropic regime
                            f_Chi2_NON_GRAY[index] = 3.0*d2Chi2_dN12[index]/2.0;
                            f_Chi2_NON_GRAY[index] = f_Chi2_NON_GRAY[index]/2.0 - 1.0;
                            f_Chi2_NON_GRAY[index] /= 1.0 - norm_f_10;
                        } else if (fabs(N1_1_NON_GRAY[index] - 1.0) < 1.0e-8) {
                            // Then free-streaming limit
                            f_Chi2_NON_GRAY[index] = 3.0*dChi2_dN1[index]/2.0;
                            f_Chi2_NON_GRAY[index] = (1.0/2.0)*(2.0 - f_Chi2_NON_GRAY[index]);
                            // cout << "id_Mobius = " << id_Mobius << "   " << "E = " << E_NON_GRAY[index] << " " << "N1 = " << N1_1_NON_GRAY[index] << " " << "Chi2 = " << Chi2_NON_GRAY[index] << " " << "f_Chi2 = " << f_Chi2_NON_GRAY[index] << " " << "dChi2_dN1 = " << dChi2_dN1[index] << endl;
                        }
                        break;
                    case RATIONAL_TYPE:
                        // f = N1^2 / (N1^2 + (1 - N1^2) g) ==> g = N1^2 (1 - f) / [f (1 - N1^2)]
                        f_Chi2_NON_GRAY[index] = (3.0*Chi2_NON_GRAY[index] - 1.0)/2.0;
                        f_Chi2_NON_GRAY[index] = norm_f_2*(1.0 - f_Chi2_NON_GRAY[index])/(f_Chi2_NON_GRAY[index]*(1.0 - norm_f_2));
                        if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
                            // Apply l'Hopital's rule twice since the first derivative is also zero
                            // in the isotropic regime
                            f_Chi2_NON_GRAY[index] = 3.0*d2Chi2_dN12[index]/2.0;
                            f_Chi2_NON_GRAY[index] = 2.0/f_Chi2_NON_GRAY[index];
                        } else if (fabs(N1_1_NON_GRAY[index] - 1.0) < 1.0e-8) {
                            // Then free-streaming limit
                            f_Chi2_NON_GRAY[index] = 3.0*dChi2_dN1[index]/2.0;
                            f_Chi2_NON_GRAY[index] = f_Chi2_NON_GRAY[index]/2.0;
                        }
                        break;
                    default:
                        cout << "Invalid Value for !!!!!!!!!!!!!!!!!" << endl;
                        exit(0);
                        break;
                }
                
//                 f_Chi2_NON_GRAY[index] = f_Chi2_NON_GRAY[index]/norm_f_2;
//                 
//                 if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) {
//                     // Apply l'Hopital's rule twice since the first derivative is also zero
//                     // in the isotropic regime
//                     f_Chi2_NON_GRAY[index] = 3.0*d2Chi2_dN12[index]/2.0;
//                     f_Chi2_NON_GRAY[index] = f_Chi2_NON_GRAY[index]/2.0;
//                 }
                
                if (f_Chi2_NON_GRAY[index] != f_Chi2_NON_GRAY[index]) {
                    cout << "In Bose-Einstein Regime" << endl;
                    cout << "E = " << E_NON_GRAY[index] << " " << "N1 = " << N1_1_NON_GRAY[index] << " " << "Chi2 = " << Chi2_NON_GRAY[index] << " " << "f_Chi2 = " << f_Chi2_NON_GRAY[index] << " " << "dChi2_dN1 = " << dChi2_dN1[index] /*<< "     " << "dChi2_dI0 = " << dChi2_dI0[index] << " " << "d2Chi2_dI0_dN1 = " << d2Chi2_dI0_dN1[index] << " " << "d3Chi2_drI0_dN1_2 = " << d3Chi2_drI0_dN1_2[index] */<< endl;
                    
                    cout << "dChi2_dN1 = " << dChi2_dN1[index] << "     " << "d2Chi2_dN12 = " << d2Chi2_dN12[index] << " " << "dChi2_drI0 = " << dChi2_drI0[index] << " " << " " << "d2Chi2_drI0_dN1 = " << d2Chi2_drI0_dN1[index] << " " << "d3Chi2_drI0_dN1_2 = " << d3Chi2_drI0_dN1_2[index] << endl;
                    exit(0);
                }
            } // end for i_Cheby_f
        } // end for i_Cheby_E
    } // end for id_Mobius
//     exit(0);
} 

void M1_Data_1D_Cheby :: SetupInterpolant_Values_BE_Impl_Derivs_N1_Derivs_r_I0(M1_Data_1D_Cheby &M1_Data_Chi2_HL, M1_Data_1D_Cheby &M1_Data_Chi2_LL) {
    int index, index_HL_LL;
    long double norm_f_2;
    long double Chi2_HL, Chi2_LL;
    long double dChi2_dN1_HL, dChi2_dN1_LL;
    long double d2_Chi2_dN12_HL, d2_Chi2_dN12_LL;
    long double d3_Chi2_dN13_HL, d3_Chi2_dN13_LL;
    
    long double diff_f_Chi2_NON_GRAY, g_Chi2_HL, g_Chi2_LL, g_Chi2_BE;
    long double theta_BE, dtheta_BE_dr;
    
    for (int id_Mobius = 0; id_Mobius < N_pts_Mob_Scale; id_Mobius++) {
        for (int i_Cheby_E = 0; i_Cheby_E < N_Points_E; i_Cheby_E++) {
            for (int i_Cheby_f = 0; i_Cheby_f < N_Points_f; i_Cheby_f++) {
                index = (id_Mobius * N_Points_E + i_Cheby_E)*N_Points_f + i_Cheby_f;
                index_HL_LL = i_Cheby_f;
                
                Chi2_HL = M1_Data_Chi2_HL.Chi2_NON_GRAY[index_HL_LL];
                Chi2_LL = M1_Data_Chi2_LL.Chi2_NON_GRAY[index_HL_LL];
                
                dChi2_dN1_HL = M1_Data_Chi2_HL.dChi2_dN1[index_HL_LL];
                dChi2_dN1_LL = M1_Data_Chi2_LL.dChi2_dN1[index_HL_LL];
                
                d2_Chi2_dN12_HL = M1_Data_Chi2_HL.d2Chi2_dN12[index_HL_LL];
                d2_Chi2_dN12_LL = M1_Data_Chi2_LL.d2Chi2_dN12[index_HL_LL];
                
                d3_Chi2_dN13_HL = M1_Data_Chi2_HL.d3Chi2_dN13[index_HL_LL];
                d3_Chi2_dN13_LL = M1_Data_Chi2_LL.d3Chi2_dN13[index_HL_LL];
                
                norm_f_2 = pow(N1_1_NON_GRAY[index],2);
                
                if ((fabs(1.0 + ratio_I0_NON_GRAY[index]) < 1.0e-4) || (fabs(1.0 - ratio_I0_NON_GRAY[index]) < 1.0e-4)) {
                    // Then Hyperbolic or Logarithmic limit
                    if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) { // Then isotropic limit
                        // Apply l'Hopital's rule twice since the first derivative is also zero
                        // in the isotropic regime
                        diff_f_Chi2_NON_GRAY = d3Chi2_drI0_dN1_2[index];
                        
                        g_Chi2_BE = d2Chi2_dN12[index];
                        
                        g_Chi2_HL = d2_Chi2_dN12_HL;
                        g_Chi2_LL = d2_Chi2_dN12_LL;
                    } else if (fabs(N1_1_NON_GRAY[index] - 1.0) < 1.0e-8) { 
                        // Then free-streaming limit
                        diff_f_Chi2_NON_GRAY = d2Chi2_drI0_dN1[index];
                        
                        g_Chi2_BE = dChi2_dN1[index];
                        
                        g_Chi2_HL = dChi2_dN1_HL;
                        g_Chi2_LL = dChi2_dN1_LL;
                    } else { // Then neither isotropic nor free-streaming limit
                        diff_f_Chi2_NON_GRAY = dChi2_drI0[index];
                        
                        g_Chi2_BE = Chi2_NON_GRAY[index];
                        
                        g_Chi2_HL = Chi2_HL;
                        g_Chi2_LL = Chi2_LL;
                    }
                    
                    // Now compute theta_BE
                    theta_BE = g_Chi2_BE/(g_Chi2_LL - g_Chi2_HL);
                    
                    // Compute partial theta_BE / partial r_{I^{(0)}}
                    dtheta_BE_dr = diff_f_Chi2_NON_GRAY/(g_Chi2_LL - g_Chi2_HL);
                        
                    //                         if (N_pts_f_L_Chi2 == 1) {
                    //                             cout << "Check if this is zero everywhere for LL f_Chi2_NON_GRAY = " << f_Chi2_NON_GRAY[index] << "  " << "diff_f_Chi2_NON_GRAY = " << diff_f_Chi2_NON_GRAY << "  " << "ratio_I0_NON_GRAY = " << ratio_I0_NON_GRAY[index] << "  " << "N1 = " << N1_1_NON_GRAY[index] << endl;
                    //                         }
                    
                    // Now compute h_BE
                    //                         if (fabs(1.0 - ratio_I0_NON_GRAY[index]) < 1.0e-4) {
                    //                             f_Chi2_NON_GRAY[index] = -dtheta_BE_dr/theta_BE;
                    //                             f_Chi2_NON_GRAY[index] = exp(f_Chi2_NON_GRAY[index]);
                    //                             f_Chi2_NON_GRAY[index] /= 1.0 + ratio_I0_NON_GRAY[index];
                    //                             
                    //                             cout << "f_Chi2 = " << f_Chi2_NON_GRAY[index] << endl;
                    //                         } else {
                    // //                             f_Chi2_NON_GRAY[index] = 2.0*dtheta_BE_dr/theta_BE;
                    // //                             f_Chi2_NON_GRAY[index] -= -log(theta_BE);
                    // //                             f_Chi2_NON_GRAY[index] /= pow(1.0 - ratio_I0_NON_GRAY[index], 2);
                    // //                             f_Chi2_NON_GRAY[index] *= exp(log(theta_BE)/2.0);
                    //                         }
                    f_Chi2_NON_GRAY[index] = (0.5 - dtheta_BE_dr)/(2.0 * ratio_I0_NON_GRAY[index]);
                } else {
                    // Then Bose-Einstein Regime
                    if (fabs(N1_1_NON_GRAY[index]) < 1.0e-8) { // Then isotropic limit
                        // Apply l'Hopital's rule twice since the first derivative is also zero
                        // in the isotropic regime
                        // Compute theta_BE
                        f_Chi2_NON_GRAY[index] = (d2Chi2_dN12[index] - d2_Chi2_dN12_HL)/(d2_Chi2_dN12_LL - d2_Chi2_dN12_HL);
                    } else if (fabs(N1_1_NON_GRAY[index] - 1.0) < 1.0e-8) { // Then free-streaming limit
                        // Compute theta_BE
                        f_Chi2_NON_GRAY[index] = (dChi2_dN1[index] - dChi2_dN1_HL)/(dChi2_dN1_LL - dChi2_dN1_HL);
                    } else { // Then neither isotropic nor free-streaming limit
                        // Compute theta_BE
                        f_Chi2_NON_GRAY[index] = (Chi2_NON_GRAY[index] - Chi2_HL)/(Chi2_LL - Chi2_HL);
                    }
                    
                    // Now compute h_BE
                    //                         f_Chi2_NON_GRAY[index] = log(f_Chi2_NON_GRAY[index])/(1.0 - ratio_I0_NON_GRAY[index]);
                    //                         f_Chi2_NON_GRAY[index] = exp(f_Chi2_NON_GRAY[index]);
                    //                         f_Chi2_NON_GRAY[index] /= 1.0 + ratio_I0_NON_GRAY[index];
                    
                    f_Chi2_NON_GRAY[index] = (f_Chi2_NON_GRAY[index] - 0.5*(1.0 + ratio_I0_NON_GRAY[index]))/(1.0 - pow(ratio_I0_NON_GRAY[index], 2));
                }
                
                if (f_Chi2_NON_GRAY[index] != f_Chi2_NON_GRAY[index]) {
                    cout << "In Bose-Einstein Regime" << endl;
                    cout << "E = " << E_NON_GRAY[index] << " " << "N1 = " << N1_1_NON_GRAY[index] << " " << "Chi2 = " << Chi2_NON_GRAY[index] << " " << "f_Chi2 = " << f_Chi2_NON_GRAY[index] << " " << "dChi2_dN1 = " << dChi2_dN1[index] /*<< "     " << "dChi2_dI0 = " << dChi2_dI0[index] << " " << "d2Chi2_dI0_dN1 = " << d2Chi2_dI0_dN1[index] << " " << "d3Chi2_drI0_dN1_2 = " << d3Chi2_drI0_dN1_2[index] */<< endl;
                    
                    cout << "dChi2_dN1 = " << dChi2_dN1[index] << "     " << "d2Chi2_dN12 = " << d2Chi2_dN12[index] << " " << "dChi2_drI0 = " << dChi2_drI0[index] << " " << " " << "d2Chi2_drI0_dN1 = " << d2Chi2_drI0_dN1[index] << " " << "d3Chi2_drI0_dN1_2 = " << d3Chi2_drI0_dN1_2[index] << endl;
                    exit(0);
                }
            } // end for i_Cheby_f
        } // end for i_Cheby_E
    } // end for id_Mobius
} 

//***********************************************************************************************
// This routine sets up and solves the Vandermonde system of equations for the purpose of our 
// proposed interpolative-based approximation of the Eddington factor
//***********************************************************************************************
void M1_Data_1D_Cheby :: Setup_Vandermonde_Matrix_2vars() {
    int index_Cheby, index_Cheby_lin, index_Cheby_col;
    int N_Points_Non_GRAY = N_Points_E*N_Points_f;
    long double ratio_E, norm_f;
    long double *VanderMonde_Matrix, *VanderMonde_Vector; 
    VanderMonde_Matrix = new long double[N_Points_Non_GRAY*N_Points_Non_GRAY];
    VanderMonde_Vector = new long double[N_Points_Non_GRAY]; 
    // Chebyshev quadrature
    for (int i_lin_Cheby_E = 0; i_lin_Cheby_E < N_Points_E; i_lin_Cheby_E++) {
        for (int i_lin_Cheby_f = 0; i_lin_Cheby_f < N_Points_f; i_lin_Cheby_f++) {
            index_Cheby_lin = i_lin_Cheby_E*N_Points_f + i_lin_Cheby_f;
            ratio_E = ratio_I0_NON_GRAY[index_Cheby_lin];
            norm_f = N1_1_NON_GRAY[index_Cheby_lin];
            
            for (int i_col_Cheby_E = 0; i_col_Cheby_E < N_Points_E; i_col_Cheby_E++) {
                for (int i_col_Cheby_f = 0; i_col_Cheby_f < N_Points_f; i_col_Cheby_f++) {
                    index_Cheby_col = i_col_Cheby_E*N_Points_f + i_col_Cheby_f;
                    
                    VanderMonde_Matrix[index_Cheby_lin*N_Points_Non_GRAY+index_Cheby_col] = VanderMonde_Matrix_2_vars(ratio_E, norm_f, i_col_Cheby_E, 2*i_col_Cheby_f);
                } // end for i_col_Cheby_f
            } // end for i_col_Cheby_E
        } // end for i_lin_Cheby_f
    } // end for i_lin_Cheby_E
    
    for (int i_Cheby_E = 0; i_Cheby_E < N_Points_E; i_Cheby_E++) {
        for (int i_Cheby_f = 0; i_Cheby_f < N_Points_f; i_Cheby_f++) {
            index_Cheby = i_Cheby_E*N_Points_f + i_Cheby_f;
            
            for (int i_Coeff_Cheby = 0; i_Coeff_Cheby < N_Points_Non_GRAY; i_Coeff_Cheby++) {
                VanderMonde_Vector[i_Coeff_Cheby] = VanderMonde_Vector_N_vars(i_Coeff_Cheby, index_Cheby);
            } // end for i_Coeff_Cheby
            
            // cout  << "*************************** solving the Vandermonde System ************************"<< endl;
            Solve_A_x_b(VanderMonde_Matrix, Coefficients_Vander_Matrix_2vars, VanderMonde_Vector, N_Points_Non_GRAY, index_Cheby);
        } // end for i_Cheby_f
    } // end for i_Cheby_E
    
    // cout << "Checking Vandermonde System !!!!!!!!!!!!!!!!!!!!" << endl;
    // Check_A_x_b(VanderMonde_Matrix, Coefficients_Vander_Matrix_2vars, N_Points_Non_GRAY);
    
    delete[] VanderMonde_Matrix;
    delete[] VanderMonde_Vector;
}

void M1_Data_1D_Cheby :: Vandermonde_Interpolation_2vars() {
    int N_Points_NON_GRAY = N_Points_E*N_Points_f;
    for (int i_Coeff_Cheby = 0; i_Coeff_Cheby < N_Points_NON_GRAY; i_Coeff_Cheby++) {
        Coefficient_Matrix_Fit_Chi2[i_Coeff_Cheby] = 0.0;
        
        for (int j_Coeff_Cheby = 0; j_Coeff_Cheby < N_Points_NON_GRAY; j_Coeff_Cheby++) {
            Coefficient_Matrix_Fit_Chi2[i_Coeff_Cheby] += Coefficients_Vander_Matrix_2vars[i_Coeff_Cheby*N_Points_NON_GRAY+j_Coeff_Cheby]*f_Chi2_NON_GRAY[j_Coeff_Cheby];
        }
        
        if (Coefficient_Matrix_Fit_Chi2[i_Coeff_Cheby] != Coefficient_Matrix_Fit_Chi2[i_Coeff_Cheby]) {
            cout << "Coefficient_Matrix_Fit_Chi2[i_Coeff_Cheby] = " <<Coefficient_Matrix_Fit_Chi2[i_Coeff_Cheby] << endl;
            
            for (int j_Coeff_Cheby = 0; j_Coeff_Cheby < N_Points_NON_GRAY; j_Coeff_Cheby++) {
                cout << "Coefficients_Vander_Matrix_2vars = " << Coefficients_Vander_Matrix_2vars[i_Coeff_Cheby*N_Points_NON_GRAY+j_Coeff_Cheby] << endl;
            }
            exit(0);
        }
    }
    
    Chebyshev_First_Kind_to_Monomial_Basis_ratio_E(Coefficient_Matrix_Fit_Chi2, N_Points_E, N_Points_f, 1);
    Chebyshev_First_Kind_to_Monomial_Basis_Norm_f(Coefficient_Matrix_Fit_Chi2, N_Points_E, N_Points_f, 1);
    
//     cout << "Checking hyperbolicity for Max Entropy Fit !!!!!!!!!!!!!!!" << endl;
//     Check_Hyperbolicity_Max_Ent_Solutions(1);
}

// long double M1_Data_1D_Cheby :: Chebyshev_First_Kind_Basis(const long double &x, const int &Index) {
//     long double Cheby, Cheby_n_m_1, Cheby_n;
//     Cheby_n_m_1 = 1.0;
//     Cheby_n = x;
//     
//     if (Index == 0) {
//         Cheby = Cheby_n_m_1;
//     } else if (Index == 1) {
//         Cheby = Cheby_n;
//     } else {
//         for (int i = 0; i < Index-1; i++) {
//             Cheby = 2.0*x*Cheby_n - Cheby_n_m_1;
//             Cheby_n_m_1 = Cheby_n;
//             Cheby_n = Cheby;
//         }
//     }
//     
//     return Cheby;
// }
// 
// long double M1_Data_1D_Cheby :: Chebyshev_First_Kind_Derivatives_Basis(const long double &x, const int &Index) {
//     long double dCheby;
//     
//     dCheby = Index*Chebyshev_Second_Kind_Basis(x, Index - 1);
//     
//     return dCheby;
// }
// 
// long double M1_Data_1D_Cheby :: Chebyshev_Second_Kind_Basis(const long double &x, const int &Index) {
//     long double Cheby, Cheby_n_m_1, Cheby_n;
//     Cheby_n_m_1 = 1.0;
//     Cheby_n = 2.0*x;
//     
//     if (Index == 0) {
//         Cheby = Cheby_n_m_1;
//     } else if (Index == 1) {
//         Cheby = Cheby_n;
//     } else {
//         for (int i = 0; i < Index-1; i++) {
//             Cheby = 2.0*x*Cheby_n - Cheby_n_m_1;
//             Cheby_n_m_1 = Cheby_n;
//             Cheby_n = Cheby;
//         }
//     }
//     
//     return Cheby;
// }

long double M1_Data_1D_Cheby :: Evaluate_diff_Length_Scale_d_sqr_N1(const long double &norm_f) {
    long double dLength_Scale, Length_Scale;
    long double norm_f_2 = norm_f*norm_f;
    int index = N_pts_L_Chi2_LS - 1;
    
    if (!flag_test_Implementation) {
        Length_Scale = Evaluate_Length_Scale(norm_f);
        
        for (int i_fit_f = N_pts_L_Chi2_LS - 1; i_fit_f >= 1; i_fit_f--) {
            if (i_fit_f == N_pts_L_Chi2_LS - 1) {
                dLength_Scale = i_fit_f*Coefficient_Matrix_Fit_Mob_Scale[i_fit_f];
            } else {
                dLength_Scale = i_fit_f*Coefficient_Matrix_Fit_Mob_Scale[i_fit_f] + dLength_Scale*norm_f_2;
            }
            index--;
        }
        
//         dLength_Scale = 0.0;
//         long double dpoly_val;
//         
//         for (int i_fit_f = 0; i_fit_f < N_pts_L_Chi2_LS; i_fit_f++) {
//             dpoly_val = Chebyshev_First_Kind_Derivatives_Basis(norm_f, 2*i_fit_f);
//             dLength_Scale += Coefficient_Matrix_Fit_Mob_Scale[i_fit_f] * dpoly_val;
//         }
        dLength_Scale *= Length_Scale;
    } else {
        dLength_Scale = 0.0;
    }
    
    return dLength_Scale;
}

long double M1_Data_1D_Cheby :: Evaluate_diff_Length_Scale_dN1(const long double &norm_f) {
    long double dLength_Scale;
    long double norm_f_2 = norm_f*norm_f;
    
    long double dL_Chi2_dnorm_f_2;
    
    // Compute \frac{\partial L_Chi2}{\partial  N1^{2}}
    dL_Chi2_dnorm_f_2 = Evaluate_diff_Length_Scale_d_sqr_N1(norm_f);
    
    // Compute \frac{\partial L_Chi2}{\partial N1} = \frac{\partial L_Chi2}{\partial N1^{2}} \frac{\partial N1^{2}}{\partial N1}
    // Compute \frac{\partial L_Chi2}{\partial N1} = 2 N1 \frac{\partial L_Chi2}{\partial N1^{2}}
    dLength_Scale = 2.0*norm_f*dL_Chi2_dnorm_f_2;
    
    return dLength_Scale;
}

long double M1_Data_1D_Cheby :: Evaluate_Length_Scale(const long double &norm_f) {
    long double Length_Scale;
    int index = N_pts_L_Chi2_LS - 1;
    
    long double norm_f_2 = norm_f*norm_f;
    
    if (!flag_test_Implementation) {
        for (int i_fit_f = N_pts_L_Chi2_LS - 1; i_fit_f >= 0; i_fit_f--) {
            if (i_fit_f == N_pts_L_Chi2_LS - 1) {
                Length_Scale = Coefficient_Matrix_Fit_Mob_Scale[i_fit_f];
            } else {
                Length_Scale = Coefficient_Matrix_Fit_Mob_Scale[i_fit_f] + Length_Scale*norm_f_2;
            }
            index--;
        }
        
//         long double poly_val;
//         Length_Scale = 0.0;
//         for (int i_fit_f = 0; i_fit_f < N_pts_L_Chi2_LS; i_fit_f++) {
//             poly_val = Chebyshev_First_Kind_Basis(norm_f, 2*i_fit_f);
//             Length_Scale += Coefficient_Matrix_Fit_Mob_Scale[i_fit_f] * poly_val;
//         }
        
        Length_Scale = exp(Length_Scale);
    } else {
        Length_Scale = 1.0;
    }
    
    // cout << "norm_f = " << norm_f << "  " << "Length_Scale = " << Length_Scale << endl;
    
    return Length_Scale;
}

// long double M1_Data_1D_Cheby :: Evaluate_dratio_E_dLength_Scale(const long double &ratio_E, const long double &Length_Scale) {
//     long double dratio_dL;
//     long double Moment;
//     
//     Moment = Inverse_Mobius_Transformation(ratio_E, Length_Scale);
//     
//     // [-(Moment + Length_Scale) - (Moment - Length_Scale)]/(Moment + Length_Scale)^2
//     // (-2*Moment)/(Moment + Length_Scale)^2
//     
//     if (fabs(fabs(ratio_E) - 1) < 1.0e-6) {
//         dratio_dL = 0.0;
//     } else {
// //         dratio_dL = -2.0*(Moment/pow(Length_Scale, 2))*exp(-Moment/Length_Scale);
//         dratio_dL = -2.0*Moment/pow(Moment + Length_Scale, 2);
//     }
//     
//     if (dratio_dL != dratio_dL) {
//         cout << "dratio_dL = " << dratio_dL << "   " << "Length_Scale = " << Length_Scale << "   " << "Moment = " << Moment << endl;
//     }
//     
//     return dratio_dL;
// }

long double M1_Data_1D_Cheby :: Evaluate_g_Chi2_HL_LL(const long double &ratio_E, const long double &norm_f, const int &Regime) {
    long double f_Chi2;
    int index = N_Points_f - 1;
    long double norm_f_2 = norm_f*norm_f;
    long double coeff_val;
    
    for (int i_fit_f = N_Points_f - 1; i_fit_f >= 0; i_fit_f--) {
        switch (Regime) {
            case HYPERBOLIC_LIMIT:
                coeff_val = Coefficient_Matrix_Fit_Chi2_HL[index];
                break;
            case LOGARITHMIC_LIMIT:
                coeff_val = Coefficient_Matrix_Fit_Chi2_LL[index];
                break;
            default:
                cout << "Invalid Regime for computing Chi2 in HL or LL" << endl;
                exit(0);
                break;
        }
        if (i_fit_f == N_Points_f - 1) {
            f_Chi2 = coeff_val;
        } else {
            f_Chi2 = coeff_val + f_Chi2*norm_f_2;
        }
        index--;
    }
    
    return f_Chi2;
}

long double M1_Data_1D_Cheby :: Evaluate_dg_Chi2_HL_LL_d_sqr_N1(const long double &ratio_E, const long double &norm_f, const int &Regime) {
    long double f_Chi2;
    int index = N_Points_f - 1;
    long double norm_f_2 = norm_f*norm_f;
    long double coeff_val;
    
    for (int i_fit_f = N_Points_f - 1; i_fit_f >= 0; i_fit_f--) {
        switch (Regime) {
            case HYPERBOLIC_LIMIT:
                coeff_val = Coefficient_Matrix_Fit_Chi2_HL[index];
                break;
            case LOGARITHMIC_LIMIT:
                coeff_val = Coefficient_Matrix_Fit_Chi2_LL[index];
                break;
            default:
                cout << "Invalid Regime for computing Chi2 in HL or LL" << endl;
                exit(0);
                break;
        }
        if (i_fit_f == N_Points_f - 1) {
            f_Chi2 = i_fit_f*coeff_val;
        } else if (i_fit_f >= 1) {
            f_Chi2 = i_fit_f*coeff_val + f_Chi2*norm_f_2;
        }
        index--;
    }
    
    return f_Chi2;
}

// long double M1_Data_1D_Cheby :: Evaluate_h_Chi2(const long double &ratio_E, const long double &norm_f) {
//     long double f_Chi2;
//     int index;
//     long double poly_map_I0_star, poly_norm_f;
//     
//     f_Chi2 = 0.0;
//     for (int i_fit_E = 0; i_fit_E < N_Points_E; i_fit_E++) {
//         for (int i_fit_f = 0; i_fit_f < N_Points_f; i_fit_f++) {
//             index = i_fit_E*N_Points_f + i_fit_f;
//             poly_map_I0_star = Chebyshev_Polynomial_Basis(ratio_E, i_fit_E);
//             poly_norm_f = Chebyshev_Polynomial_Basis(norm_f, 2*i_fit_f);
//             
//             f_Chi2 += Coefficient_Matrix_Fit_Chi2[index]*poly_map_I0_star*poly_norm_f; 
//         }
//     }
//     
//     return f_Chi2;
// }

long double M1_Data_1D_Cheby :: Evaluate_h_Chi2(const long double &ratio_E, const long double &norm_f) {
    long double f_Chi2, f_Chi2_temp_f;
    int index = N_Points_E*N_Points_f - 1;
    long double norm_f_2 = norm_f*norm_f;
    
    for (int i_fit_E = N_Points_E - 1; i_fit_E >= 0; i_fit_E--) {
        for (int i_fit_f = N_Points_f - 1; i_fit_f >= 0; i_fit_f--) {
            if (i_fit_f == N_Points_f - 1) {
                f_Chi2_temp_f = Coefficient_Matrix_Fit_Chi2[index];
            } else {
                f_Chi2_temp_f = Coefficient_Matrix_Fit_Chi2[index] + f_Chi2_temp_f*norm_f_2;
            }
            index--;
        }
        if (i_fit_E == N_Points_E - 1) {
            f_Chi2 = f_Chi2_temp_f;
        } else {
            f_Chi2 = f_Chi2_temp_f + f_Chi2*ratio_E;   
        }
    }
    
    if (f_Chi2 != f_Chi2) {
        index = N_Points_E*N_Points_f - 1;
        for (int i_fit_E = N_Points_E - 1; i_fit_E >= 0; i_fit_E--) {
            for (int i_fit_f = N_Points_f - 1; i_fit_f >= 0; i_fit_f--) {
                cout << "Coefficient_Matrix_Fit_Chi2 = " << Coefficient_Matrix_Fit_Chi2[index] << endl;
            }
        }
        exit(0);
    }
    
    return f_Chi2;
}

long double M1_Data_1D_Cheby :: Evaluate_dh_Chi2_dratio_E(const long double &ratio_E, const long double &norm_f) {
    long double Chi2, Chi2_temp_f;
    int index = N_Points_E*N_Points_f - 1;
    long double norm_f_2 = norm_f*norm_f;
    
    for (int i_fit_E = N_Points_E - 1; i_fit_E >= 0; i_fit_E--) {
        for (int i_fit_f = N_Points_f - 1; i_fit_f >= 0; i_fit_f--) {
            if (i_fit_f == N_Points_f - 1) {
                Chi2_temp_f = Coefficient_Matrix_Fit_Chi2[index];
            } else {
                Chi2_temp_f = Coefficient_Matrix_Fit_Chi2[index] + Chi2_temp_f*norm_f_2;
            }
            index--;
        }
        if (i_fit_E == N_Points_E - 1) {
            Chi2 = i_fit_E*Chi2_temp_f;
        } else if (i_fit_E >= 1) {
            Chi2 = i_fit_E*Chi2_temp_f + Chi2*ratio_E;   
        }
    }
    
    return Chi2;
}

long double M1_Data_1D_Cheby :: Evaluate_dh_Chi2_d_sqr_N1(const long double &ratio_E, const long double &norm_f) {
    long double Chi2, Chi2_temp_f;
    int index = N_Points_E*N_Points_f - 1;
    long double norm_f_2 = norm_f*norm_f;
    
    index = N_Points_E*N_Points_f - 1;
    for (int i_fit_E = N_Points_E - 1; i_fit_E >= 0; i_fit_E--) {
        for (int i_fit_f = N_Points_f - 1; i_fit_f >= 0; i_fit_f--) {
            if (i_fit_f == N_Points_f - 1) {
                Chi2_temp_f = i_fit_f*Coefficient_Matrix_Fit_Chi2[index];
            } else if (i_fit_f >= 1) {
                Chi2_temp_f = i_fit_f*Coefficient_Matrix_Fit_Chi2[index] + Chi2_temp_f*norm_f_2;
            }
            index--;
        }
        if (i_fit_E == N_Points_E - 1) {
            Chi2 = Chi2_temp_f;
        } else {
            Chi2 = Chi2_temp_f + Chi2*ratio_E;   
        }
    }
    
    return Chi2;
}

long double M1_Data_1D_Cheby :: Evaluate_theta_Chi2(const long double &ratio_E, const long double &norm_f) {
    long double h_BE, theta_BE;
    
    if (!flag_test_Implementation) {
        h_BE = Evaluate_h_Chi2(ratio_E, norm_f);
        switch (Implementation_type) {
            case IMPLEMENTATION_DERIVATIVES_N1:
                theta_BE = h_BE;
                break;
            case IMPLEMENTATION_DERIVATIVES_N1_DERIVATIVES_R_I0:
                // h_BE = log(h_BE);
                
                // (1 + r) [ 0.5 + (1 - r) theta_BE]
                theta_BE = (1.0 + ratio_E) * (0.5 + (1.0 - ratio_E) * h_BE);
                break;
            default:
                cout << "Implementation type not specified, value is " << Implementation_type << endl;
                exit(0);
                break;
        }
    } else {
        switch (Implementation_type) {
            case IMPLEMENTATION_DERIVATIVES_N1:
                theta_BE = 0.0;
                break;
            case IMPLEMENTATION_DERIVATIVES_N1_DERIVATIVES_R_I0:
                theta_BE = ratio_E;
                break;
            default:
                cout << "Implementation type not specified, value is " << Implementation_type << endl;
                exit(0);
                break;
        }
    }
    
    return theta_BE;
}

long double M1_Data_1D_Cheby :: Evaluate_dtheta_Chi2_dratio_E(const long double &ratio_E, const long double &norm_f) {
    long double h_BE, dh_BE, dtheta_BE;
    
    if (!flag_test_Implementation) {
        dh_BE = Evaluate_dh_Chi2_dratio_E(ratio_E, norm_f);
        switch (Implementation_type) {
            case IMPLEMENTATION_DERIVATIVES_N1:
                dtheta_BE = dh_BE;
                break;
            case IMPLEMENTATION_DERIVATIVES_N1_DERIVATIVES_R_I0:
                h_BE = Evaluate_h_Chi2(ratio_E, norm_f);
                
                // dh_BE = dh_BE/h_BE;
                
                dtheta_BE = 0.5 - 2.0*ratio_E * h_BE + (1.0 - pow(ratio_E,2)) * dh_BE;
                break;
            default:
                cout << "Implementation type not specified, value is " << Implementation_type << endl;
                exit(0);
                break;
        }
    } else {
        switch (Implementation_type) {
            case IMPLEMENTATION_DERIVATIVES_N1:
                dtheta_BE = 0.0;
                break;
            case IMPLEMENTATION_DERIVATIVES_N1_DERIVATIVES_R_I0:
                dtheta_BE = 1.0;
                break;
            default:
                cout << "Implementation type not specified, value is " << Implementation_type << endl;
                exit(0);
                break;
        }
    }
    
    return dtheta_BE;
}

long double M1_Data_1D_Cheby :: Evaluate_dtheta_Chi2_d_sqr_N1(const long double &ratio_E, const long double &norm_f) {
    long double h_BE, dh_BE, dtheta_BE;
    
    if (!flag_test_Implementation) {
        dh_BE = Evaluate_dh_Chi2_d_sqr_N1(ratio_E, norm_f);
        switch (Implementation_type) {
            case IMPLEMENTATION_DERIVATIVES_N1:
                dtheta_BE = dh_BE;
                break;
            case IMPLEMENTATION_DERIVATIVES_N1_DERIVATIVES_R_I0:
                // h_BE = Evaluate_h_Chi2(ratio_E, norm_f);
                
                // dh_BE = dh_BE/h_BE;
                
                // (1 + r) [ 0.5 + (1 - r) theta_BE]
                dtheta_BE = (1.0 - pow(ratio_E,2)) * dh_BE;
                break;
            default:
                cout << "Implementation type not specified, value is " << Implementation_type << endl;
                exit(0);
                break;
        }
    } else {
        switch (Implementation_type) {
            case IMPLEMENTATION_DERIVATIVES_N1:
                dtheta_BE = 0.0;
                break;
            case IMPLEMENTATION_DERIVATIVES_N1_DERIVATIVES_R_I0:
                dtheta_BE = 0.0;
                break;
            default:
                cout << "Implementation type not specified, value is " << Implementation_type << endl;
                exit(0);
                break;
        }
    }
    
    return dtheta_BE;
}

long double M1_Data_1D_Cheby :: Evaluate_g_Chi2(const long double &ratio_E, const long double &norm_f) {
    long double g_Chi2, g_Chi2_HL, g_Chi2_LL, theta_BE;
    
    theta_BE = Evaluate_theta_Chi2(ratio_E, norm_f);
            
    switch (Implementation_type) {
        case IMPLEMENTATION_DERIVATIVES_N1:
            g_Chi2 = theta_BE;
            break;
        case IMPLEMENTATION_DERIVATIVES_N1_DERIVATIVES_R_I0:
            g_Chi2_HL = Evaluate_g_Chi2_HL_LL(ratio_E, norm_f, HYPERBOLIC_LIMIT);
            g_Chi2_LL = Evaluate_g_Chi2_HL_LL(ratio_E, norm_f, LOGARITHMIC_LIMIT);
            g_Chi2 = g_Chi2_HL + (g_Chi2_LL - g_Chi2_HL)*theta_BE;
            break;
        default:
            cout << "Implementation type not specified, value is " << Implementation_type << endl;
            exit(0);
            break;
    }
    
    return g_Chi2;
}

long double M1_Data_1D_Cheby :: Evaluate_dg_Chi2_dratio_E(const long double &ratio_E, const long double &norm_f) {
    long double dg_Chi2, g_Chi2_HL, g_Chi2_LL, dtheta_BE_d_rI0;
    
    dtheta_BE_d_rI0 = Evaluate_dtheta_Chi2_dratio_E(ratio_E, norm_f);
    
    switch (Implementation_type) {
        case IMPLEMENTATION_DERIVATIVES_N1:
            dg_Chi2 = dtheta_BE_d_rI0;
            break;
        case IMPLEMENTATION_DERIVATIVES_N1_DERIVATIVES_R_I0:
            g_Chi2_HL = Evaluate_g_Chi2_HL_LL(ratio_E, norm_f, HYPERBOLIC_LIMIT);
            g_Chi2_LL = Evaluate_g_Chi2_HL_LL(ratio_E, norm_f, LOGARITHMIC_LIMIT);
            dg_Chi2 = (g_Chi2_LL - g_Chi2_HL)*dtheta_BE_d_rI0;
            break;
        default:
            cout << "Implementation type not specified, value is " << Implementation_type << endl;
            exit(0);
            break;
    }
    
    return dg_Chi2;
}

long double M1_Data_1D_Cheby :: Evaluate_dg_Chi2_d_sqr_N1(const long double &ratio_E, const long double &norm_f) {
    long double dg_Chi2, dg_Chi2_HL, dg_Chi2_LL, g_Chi2_HL, g_Chi2_LL, theta_BE, dtheta_BE_d_sqr_N1;
    
    dtheta_BE_d_sqr_N1 = Evaluate_dtheta_Chi2_d_sqr_N1(ratio_E, norm_f);
    
    switch (Implementation_type) {
        case IMPLEMENTATION_DERIVATIVES_N1:
            dg_Chi2 = dtheta_BE_d_sqr_N1;
            break;
        case IMPLEMENTATION_DERIVATIVES_N1_DERIVATIVES_R_I0:
            g_Chi2_HL = Evaluate_g_Chi2_HL_LL(ratio_E, norm_f, HYPERBOLIC_LIMIT);
            g_Chi2_LL = Evaluate_g_Chi2_HL_LL(ratio_E, norm_f, LOGARITHMIC_LIMIT);
            
            dg_Chi2_HL = Evaluate_dg_Chi2_HL_LL_d_sqr_N1(ratio_E, norm_f, HYPERBOLIC_LIMIT);
            dg_Chi2_LL = Evaluate_dg_Chi2_HL_LL_d_sqr_N1(ratio_E, norm_f, LOGARITHMIC_LIMIT);
            theta_BE = Evaluate_theta_Chi2(ratio_E, norm_f);
            
            dg_Chi2 = dg_Chi2_HL + (dg_Chi2_LL - dg_Chi2_HL)*theta_BE + (g_Chi2_LL - g_Chi2_HL)*dtheta_BE_d_sqr_N1;
            break;
        default:
            cout << "Implementation type not specified, value is " << Implementation_type << endl;
            exit(0);
            break;
    }
    
    return dg_Chi2;
}

long double M1_Data_1D_Cheby :: Evaluate_f_Chi2(const long double &ratio_E, const long double &norm_f) {
    long double f_Chi2, g_Chi2;
    long double norm_f_2 = norm_f*norm_f;
    long double norm_f_10 = pow(norm_f_2, 1);
    
    g_Chi2 = Evaluate_g_Chi2(ratio_E, norm_f);
    
    switch (Weighting_Function_fchi2_Interpolation_Type) {
        case POLYNOMIAL_TYPE:
            f_Chi2 = norm_f_2*(1.0 + (1.0 - norm_f_10)*g_Chi2);
            break;
        case RATIONAL_TYPE:
            // f = N1^2 / (N1^2 + (1 - N1^2) g)
            f_Chi2 = norm_f_2/(norm_f_2 + (1.0 - norm_f_2)*g_Chi2);
            break;
        default:
            cout << "Invalid Value for !!!!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    }
    
    // f_Chi2 = norm_f_2*g_Chi2;
    
    return f_Chi2;
}

long double M1_Data_1D_Cheby :: Evaluate_df_Chi2_dratio_E(const long double &ratio_E, const long double &norm_f) {
    long double d_f_Chi2, d_g_Chi2, g_Chi2;
    long double norm_f_2 = norm_f*norm_f;
    long double norm_f_10 = pow(norm_f_2, 1);
    
    g_Chi2 = Evaluate_g_Chi2(ratio_E, norm_f);
    d_g_Chi2 = Evaluate_dg_Chi2_dratio_E(ratio_E, norm_f);
    
    // f_chi2 = (N1)^2 * [1 + (1 - (N1)^2) g_chi2]
    // f_chi2' = (N1)^2 (1 - (N1)^2) g_chi2'
    // where f_chi2' = \frac{\partial f_chi2}{\partial ratio_E}
    // and g_chi2' = \frac{\partial g_chi2}{\partial ratio_E}
    
    switch (Weighting_Function_fchi2_Interpolation_Type) {
        case POLYNOMIAL_TYPE:
            d_f_Chi2 = norm_f_2*(1.0 - norm_f_10)*d_g_Chi2;
            break;
        case RATIONAL_TYPE:
            // f = N1^2 / (N1^2 + (1 - N1^2) g)
            d_f_Chi2 = -norm_f_2*(1.0 - norm_f_2)*d_g_Chi2/pow((norm_f_2 + (1.0 - norm_f_2)*g_Chi2), 2);
            break;
        default:
            cout << "Invalid Value for !!!!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    }
    
    // d_f_Chi2 = norm_f_2*d_g_Chi2;
    
    return d_f_Chi2;
}

long double M1_Data_1D_Cheby :: Evaluate_df_Chi2_d_sqr_N1(const long double &ratio_E, const long double &norm_f) {
    long double d_f_Chi2, g_Chi2, d_g_Chi2;
    long double norm_f_2 = norm_f*norm_f;
    long double norm_f_10 = pow(norm_f_2, 1);
    
    g_Chi2 = Evaluate_g_Chi2(ratio_E, norm_f);
    
    d_g_Chi2 = Evaluate_dg_Chi2_d_sqr_N1(ratio_E, norm_f);
    
    // f_chi2 = (N1)^2 * [1 + (1 - (N1)^2) g_chi2]
    // f_chi2' = (1 + (1 - 2 (N1)^2 ) g_chi2 + (N1)^2 (1 - (N1)^2) g_chi2')
    // where f_chi2' = \frac{\partial f_chi2}{\partial N1^{2}}
    // and g_chi2' = \frac{\partial g_chi2}{\partial N1^{2}}
    
    switch (Weighting_Function_fchi2_Interpolation_Type) {
        case POLYNOMIAL_TYPE:
            d_f_Chi2 = 1.0 + (1.0 - 2.0*norm_f_10)*g_Chi2 + norm_f_2*(1.0 - norm_f_10)*d_g_Chi2;
            break;
        case RATIONAL_TYPE:
            // f = N1^2 / (N1^2 + (1 - N1^2) g)
            d_f_Chi2 = (g_Chi2 - norm_f_2*(1.0 - norm_f_2)*d_g_Chi2)/pow((norm_f_2 + (1.0 - norm_f_2)*g_Chi2), 2);
            break;
        default:
            cout << "Invalid Value for !!!!!!!!!!!!!!!!!" << endl;
            exit(0);
            break;
    }
    
    // d_f_Chi2 = g_Chi2 + norm_f_2*d_g_Chi2;
    
    return d_f_Chi2;
}

//*************************************************************************************
// This routine computes the interpolative-based approximation of the Eddington factor
// in the case of the non-gray M1 closure
//*************************************************************************************
long double M1_Data_1D_Cheby :: Evaluate_Chi2(const long double &ratio_E, const long double &norm_f) {
    long double Chi2, f_Chi2;
    
    long double zeta;
    switch (Problem_Type) {
        case GRAY:
            zeta = sqrt(4.0 - 3.0*pow(norm_f, 2));
            Chi2 = 5.0 - 2.0*zeta;
            Chi2 /= 3.0;
            break;
        case NON_GRAY:
            f_Chi2 = Evaluate_f_Chi2(ratio_E, norm_f);
            Chi2 = (1.0 + 2.0*f_Chi2)/3.0;
            break;
        default:
            cout << "Problem Type not specified" << endl;
            exit(0);
            break;
    }
    
    return Chi2;
}

long double M1_Data_1D_Cheby :: Evaluate_dChi2_dratio_E(const long double &ratio_E, const long double &norm_f) {
    long double d_chi2, d_f_chi2;
    
    switch (Problem_Type) {
        case GRAY:
            d_f_chi2 = 0.0;
            break;
        case NON_GRAY:
            // ratio_E and norm_f_2 are treated as independent variables for now
            d_f_chi2 = Evaluate_df_Chi2_dratio_E(ratio_E, norm_f);
            break;
        default:
            cout << "Problem Type not specified" << endl;
            exit(0);
            break;
    }
    
    d_chi2 = 2.0*d_f_chi2/3.0;
    
    return d_chi2;
}

long double M1_Data_1D_Cheby :: Evaluate_dChi2_dN1(const long double &ratio_E, const long double &norm_f) {
    long double d_chi2, d_f_chi2;
    
    // ratio_E and norm_f_2 are treated as independent variables for now
    // Now we compute d_f_Chi2_dnorm_f_2
    
    long double zeta;
    switch (Problem_Type) {
        case GRAY:
            zeta = sqrt(4.0 - 3.0*pow(norm_f, 2));
            d_chi2 = 2.0*norm_f/zeta;
            break;
        case NON_GRAY:
            d_f_chi2 = Evaluate_df_Chi2_d_sqr_N1(ratio_E, norm_f);
            
            // dChi2_dnorm_f = dChi2_dnorm_f_2*dnorm_f_2_dnorm_f
            // dChi2_dnorm_f = 2.0*norm_f*dChi2_dnorm_f_2
            d_f_chi2 *= 2.0*norm_f;
            
            d_chi2 = 2.0*d_f_chi2/3.0;
            break;
        default:
            cout << "Problem Type not specified" << endl;
            exit(0);
            break;
    }
    
    return d_chi2;
}


void M1_Data_1D_Cheby :: Evaluate_Chi2_derivatives(long double &E_dChi2_dE, long double &dChi2_df, const long double &r, const long double &N1) {
    long double Mobius_Scale_Fit;
    long double dMobius_Scale_Fit_df, dmap_E_dE;
    long double dmap_E_df;
    long double dChi2_dmapE;
    long double E;
    
    switch (Problem_Type) {
        case GRAY:
            dChi2_df = Evaluate_dChi2_dN1(r, N1);
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
            E = Inverse_Mobius_Transformation(r, Mobius_Scale_Fit);
            
            if (r == 1.0) {
                E = 1.0e32;
            }
            
            // Compute \frac{\partial L_Chi2}{\partial N1}
            switch (Least_Squares_L_Chi2_Mode) {
                case LEAST_SQUARES_L_CHI2_ON:
                    dMobius_Scale_Fit_df = Evaluate_diff_Length_Scale_dN1_Least_Squares(N1);
                    break;
                case LEAST_SQUARES_L_CHI2_OFF:
                    dMobius_Scale_Fit_df = Evaluate_diff_Length_Scale_dN1(N1);
                    break;
            };
            
            // Compute \frac{\partial ratio_E}{\partial I0}
            dmap_E_dE = d_ratio_E_d_I0(E, Mobius_Scale_Fit);
            
            // Compute \frac{\partial ratio_E}{\partial N1} = \frac{\partial ratio_E}{\partial L_Chi2} \frac{\partial L_Chi2}{\partial N1}
            dmap_E_df = d_ratio_E_d_L_Chi2(E, Mobius_Scale_Fit);
            dmap_E_df *= dMobius_Scale_Fit_df;
            
            // Compute \frac{\partial chi2}{\partial ratio_E}
            dChi2_dmapE = Evaluate_dChi2_dratio_E(r, N1);
            
            // Compute \frac{\partial chi2}{\partial N1}
            dChi2_df = Evaluate_dChi2_dN1(r, N1);
            
            // Since ratio_E = ratio_E(I0, N1), it follows that
            // \frac{d chi2}{d N1} = \frac{\partial chi2}{\partial N1} + (\frac{\partial chi2}{\partial ratio_E}) (\frac{\partial ratio_E}{\partial N1})
            dChi2_df = dChi2_df + dChi2_dmapE*dmap_E_df;
            
            // Also, since ratio_E = ratio_E(I0, N1) and N1 = I1/I0, it follows that
            // \frac{d chi2}{d I0} = \frac{\partial chi2}{\partial N1} \frac{\partial N1}{\partial I0} + 
            //                       (\frac{\partial chi2}{\partial ratio_E}) (\frac{\partial ratio_E}{\partial I0})
            // ==> \frac{d chi2}{d I0} = - \frac{\partial chi2}{\partial N1} (N1/I0) + 
            //                           (\frac{\partial chi2}{\partial ratio_E}) (\frac{\partial ratio_E}{\partial I0})
            // ==> I0 \frac{d chi2}{d I0} = - N1 \frac{\partial chi2}{\partial N1} + 
            //                             I0 (\frac{\partial chi2}{\partial ratio_E}) (\frac{\partial ratio_E}{\partial I0})
            E_dChi2_dE = E*dChi2_dmapE*dmap_E_dE - dChi2_df*N1;
            break;
        default:
            cout << "Problem Type not specified" << endl;
            exit(0);
            break;
    }
    
    // cout << "E = " << E << "    " << "dmap_E_dE = " << dmap_E_dE << "    " << "dChi2_dmapE = " << dChi2_dmapE << endl;
    
//     long double E_dChi2_dE_temp, dChi2_df_temp;
//     Evaluate_Chi2_derivatives_Finite_Difference(E_dChi2_dE_temp, dChi2_df_temp, r, N1);
//     E_dChi2_dE = E_dChi2_dE_temp;
//     dChi2_df = dChi2_df_temp;
//     if (fabs(E_dChi2_dE_temp - E_dChi2_dE) > 1.0e-4 || fabs(dChi2_df_temp - dChi2_df) > 1.0e-4) {
//         cout << "E_dChi2_dE = " << E_dChi2_dE << "  " << "E_dChi2_dE_temp = " << E_dChi2_dE_temp << "  " << "diff = " << E_dChi2_dE_temp - E_dChi2_dE << endl;
//         cout << "dChi2_df = " << dChi2_df << "  " << "dChi2_df_temp = " << dChi2_df_temp << "  " << "diff = " << dChi2_df_temp - dChi2_df << endl;
//         // exit(0);
//     }
}

void M1_Data_1D_Cheby :: Evaluate_Chi2_derivatives_Finite_Difference(long double &E_dChi2_dE, long double &dChi2_df, const long double &r, const long double &N1) {
    long double Mobius_Scale_Fit;
    long double dMobius_Scale_Fit_df, dmap_E_dE;
    long double dmap_E_df;
    long double dChi2_dmapE;
    long double E;
    long double L_Chi2_L, L_Chi2_R;
    
    long double Chi2_L, Chi2_R;
    long double map_E_L, map_E_R;
    long double epsilon = 1.0e-6;
    
    switch (Problem_Type) {
        case GRAY:
            // Compute \frac{\partial chi2}{\partial N1}
            Chi2_L = Evaluate_Chi2(r, N1 - epsilon);
            Chi2_R = Evaluate_Chi2(r, N1 + epsilon);
            dChi2_df = (Chi2_R - Chi2_L)/(2.0*epsilon);
            
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
            E = Inverse_Mobius_Transformation(r, Mobius_Scale_Fit);
            
            if (r == 1.0) {
                E = 1.0e32;
            }
            
            // Compute \frac{\partial L_Chi2}{\partial N1}
            switch (Least_Squares_L_Chi2_Mode) {
                case LEAST_SQUARES_L_CHI2_ON:
                    L_Chi2_L = Evaluate_Length_Scale_Least_Squares(N1 - epsilon);
                    L_Chi2_R = Evaluate_Length_Scale_Least_Squares(N1 + epsilon);
                    break;
                case LEAST_SQUARES_L_CHI2_OFF:
                    L_Chi2_L = Evaluate_Length_Scale(N1 - epsilon);
                    L_Chi2_R = Evaluate_Length_Scale(N1 + epsilon);
                    break;
            };
            dMobius_Scale_Fit_df = (L_Chi2_R - L_Chi2_L)/(2.0*epsilon);
            
            // Compute \frac{\partial ratio_E}{\partial I0}
            map_E_L = Mobius_Transformation(E - epsilon, Mobius_Scale_Fit);
            map_E_R = Mobius_Transformation(E + epsilon, Mobius_Scale_Fit);
            dmap_E_dE = (map_E_R - map_E_L)/(2.0*epsilon);
            
            // Compute \frac{\partial ratio_E}{\partial N1} = \frac{\partial ratio_E}{\partial L_Chi2} \frac{\partial L_Chi2}{\partial N1}
            map_E_L = Mobius_Transformation(E, Mobius_Scale_Fit - epsilon);
            map_E_R = Mobius_Transformation(E, Mobius_Scale_Fit + epsilon);
            dmap_E_df = (map_E_R - map_E_L)/(2.0*epsilon);
            dmap_E_df *= dMobius_Scale_Fit_df;
            
            // Compute \frac{\partial chi2}{\partial ratio_E}
            Chi2_L = Evaluate_Chi2(r - epsilon, N1);
            Chi2_R = Evaluate_Chi2(r + epsilon, N1);
            dChi2_dmapE = (Chi2_R - Chi2_L)/(2.0*epsilon);
            
            // Compute \frac{\partial chi2}{\partial N1}
            Chi2_L = Evaluate_Chi2(r, N1 - epsilon);
            Chi2_R = Evaluate_Chi2(r, N1 + epsilon);
            dChi2_df = (Chi2_R - Chi2_L)/(2.0*epsilon);
            
            // Since ratio_E = ratio_E(I0, N1), it follows that
            // \frac{d chi2}{d N1} = \frac{\partial chi2}{\partial N1} + (\frac{\partial chi2}{\partial ratio_E}) (\frac{\partial ratio_E}{\partial N1})
            dChi2_df = dChi2_df + dChi2_dmapE*dmap_E_df;
            
            // Also, since ratio_E = ratio_E(I0, N1) and N1 = I1/I0, it follows that
            // \frac{d chi2}{d I0} = \frac{\partial chi2}{\partial N1} \frac{\partial N1}{\partial I0} + 
            //                       (\frac{\partial chi2}{\partial ratio_E}) (\frac{\partial ratio_E}{\partial I0})
            // ==> \frac{d chi2}{d I0} = - \frac{\partial chi2}{\partial N1} (N1/I0) + 
            //                           (\frac{\partial chi2}{\partial ratio_E}) (\frac{\partial ratio_E}{\partial I0})
            // ==> I0 \frac{d chi2}{d I0} = - N1 \frac{\partial chi2}{\partial N1} + 
            //                             I0 (\frac{\partial chi2}{\partial ratio_E}) (\frac{\partial ratio_E}{\partial I0})
            E_dChi2_dE = E*dChi2_dmapE*dmap_E_dE - dChi2_df*N1;
            break;
        default:
            cout << "Problem Type not specified" << endl;
            exit(0);
            break;
    }
}

// ******************************************************************************
// This routine performs the polynomial interpolation for the Eddington factor
// in either the Hyperbolic or the Logarithmic limit
// ******************************************************************************
void M1_Data_1D_Cheby :: Polynomial_Interpolation_HL_LL(const int &Maximum_Entropy_Solution_Regime) {
    // First Precompute maximum entropy solutions
    Precompute_Final_Max_Ent_Solution(Maximum_Entropy_Solution_Regime);
    
    // Setup interpolant values for the purpose of polynomial interpolation
    SetupInterpolant_Values_HL_LL();
    
//     int index;
//     for (int i_Cheby_f = 0; i_Cheby_f < N_Points_f; i_Cheby_f++) {
//         index = i_Cheby_f;
//         cout << "E = " << E_NON_GRAY[index] << "  " << "ratio_I0_NON_GRAY = " << ratio_I0_NON_GRAY[index] << "  " << "N1 = " << N1_1_NON_GRAY[index] << "  " << "f_Chi2 = " << f_Chi2_NON_GRAY[index] << endl;
//     }
    
    // Setup the Vandermonde matrix used to compute the interpolant
    Setup_Vandermonde_Matrix_2vars();
    
    // solve the Vandermonde system for the coefficients of the interpolant
    Vandermonde_Interpolation_2vars();
}

// ******************************************************************************
// This routine performs the polynomial interpolation for the Eddington factor
// in the Bose Einstein Regime
// ******************************************************************************
void M1_Data_1D_Cheby :: Polynomial_Interpolation_BE(M1_Data_1D_Cheby &M1_Data_1D_Cheby_HL, M1_Data_1D_Cheby &M1_Data_1D_Cheby_LL) {
    // First Precompute maximum entropy solutions
    Precompute_Final_Max_Ent_Solution();
    
    // Setup interpolant values for the purpose of polynomial interpolation
    SetupInterpolant_Values_BE(M1_Data_1D_Cheby_HL, M1_Data_1D_Cheby_LL);
    
//     int index;
//     for (int i_Cheby_E = 0; i_Cheby_E < N_Points_E; i_Cheby_E++) {
//         for (int i_Cheby_f = 0; i_Cheby_f < N_Points_f; i_Cheby_f++) {
//             index = i_Cheby_E*N_Points_f + i_Cheby_f;
//             cout << "E = " << E_NON_GRAY[index] << "  " << "ratio_I0_NON_GRAY = " << ratio_I0_NON_GRAY[index] << "  " << "N1 = " << N1_1_NON_GRAY[index] << "  " << "f_Chi2 = " << f_Chi2_NON_GRAY[index] << endl;
//         }
//     }
    
    // Setup coefficients for the interpolation in both the Hyperbolic and Logarithmic limits
    Setup_Coefficients_HL_LL(M1_Data_1D_Cheby_HL, M1_Data_1D_Cheby_LL);
    
    // Setup the Vandermonde matrix used to compute the interpolant
    Setup_Vandermonde_Matrix_2vars();
    
    // solve the Vandermonde system for the coefficients of the interpolant
    Vandermonde_Interpolation_2vars();
}

// ******************************************************************************
// This routine computes optimal values of the length scale f_L_Chi2 of the 
// algebraic mapping of the first-order normalized angular moment
// ******************************************************************************
void M1_Data_1D_Cheby :: Polynomial_Interpolation_Non_Gray_M1(M1_1D_Data_Pointer &M1_1D_Data, ofstream &output_Opt_Coefficients) {
    ofstream out_L_inf_Norm;
    
    // Testing finite difference
    // Test_Finite_Difference();
    
    // First perform polynomial interpolation in the Hyperbolic and Free-Streaming limits
    M1_Data_1D_Cheby M1_Data_1D_Cheby_HL(1, N_Points_f);
    M1_Data_1D_Cheby M1_Data_1D_Cheby_LL(1, N_Points_f);
    
    // Hyperbolic Limit
    M1_Data_1D_Cheby_HL.Polynomial_Interpolation_HL_LL(HYPERBOLIC_LIMIT);
    
    // Now Logarithmic Limit
    M1_Data_1D_Cheby_LL.Polynomial_Interpolation_HL_LL(LOGARITHMIC_LIMIT);
    
    M1_1D_Data.M1_Data_HL = &M1_Data_1D_Cheby_HL;
    M1_1D_Data.M1_Data_LL = &M1_Data_1D_Cheby_LL;
    
    // Compute optimal set of coefficients for the interpolative-based approximation of
    // the length of the exponential mapping of the radiative energy density
    Nested_Least_Squares_Optimization_L_Chi2(M1_1D_Data);
    
    // Perform polynomial interpolation for the Eddington factor based on the optimal set 
    // of coefficients for the interpolative-based approximation of the length of the 
    // exponential mapping of the radiative energy density
    Polynomial_Interpolation_BE(M1_Data_1D_Cheby_HL, M1_Data_1D_Cheby_LL);
    
    if (id_proc == 0) {
        // Test to make sure the interpolation is exact at the interpolation nodes
        // The error should be zero
        cout  << "************************ Testing fit in 1D based on selected node for interpolation "<< endl;
        Fit_Convergence_Check(0);
        cout  << "****************** Testing fit in 1D based on selected node for interpolation Completed "<< endl;
        
        cout << endl;
        
        // Test to assess the error of our interpolation for a set of evaluation points spanning
        // the full realizable space for angular moments up to first-order
        cout  << "************************ Testing accuracy of interpolation in 1D "<< endl;
        Fit_Convergence_Check(M1_1D_Data.M1_Data_Uniform_BE, 0, out_L_inf_Norm, 0);
        cout  << "************************ Testing accuracy of interpolation in 1D Completed "<< endl;
        
        // Write coefficients for our interpolative-based approximation of the 
        // non-gray M1 closure
        Write_Coefficients_M1_Closure_Interp(output_Opt_Coefficients);
        
        Write_Output_Data_Fit_Chi2_Matlab();
    }
}

void M1_Data_1D_Cheby :: Write_Coefficients_M1_Closure_Interp(ofstream &output_Opt_Coefficients) {
    int index_Coeffs_M1_Closure, index_Coeffs_M1_Closure_max;
    
    if (Implementation_type == IMPLEMENTATION_DERIVATIVES_N1_DERIVATIVES_R_I0) {
        // Output the coefficients for the polynomial interpolation approximation of the
        // Eddington factor in HL and LL
        for (int i_f = 0; i_f < N_Points_f; i_f++) {
            output_Opt_Coefficients << Coefficient_Matrix_Fit_Chi2_HL[i_f];
            if (i_f < N_Points_f - 1) {
                output_Opt_Coefficients << setw(16);
            }
        }
        output_Opt_Coefficients << endl;
        
        for (int i_f = 0; i_f < N_Points_f; i_f++) {
            output_Opt_Coefficients << Coefficient_Matrix_Fit_Chi2_LL[i_f];
            if (i_f < N_Points_f - 1) {
                output_Opt_Coefficients << setw(16);
            }
        }
        output_Opt_Coefficients << endl;
    }
    
    // Output the coefficients for the polynomial interpolation approximation of the
    // optimal length scale
    for (int i_f = 0; i_f < N_pts_L_Chi2_LS; i_f++) {
        output_Opt_Coefficients << setprecision(12) << Coefficient_Matrix_Fit_Mob_Scale[i_f];
        if (i_f < N_pts_L_Chi2_LS - 1) {
            output_Opt_Coefficients << setw(16);
        }
    }
    
    output_Opt_Coefficients << endl;
    
    // Output the coefficients for the polynomial interpolation approximation of the
    // Eddington factor in BE
    index_Coeffs_M1_Closure_max = N_Points_E*N_Points_f - 1;
    for (int i_E = 0; i_E < N_Points_E; i_E++) {
        for (int i_f = 0; i_f < N_Points_f; i_f++) {
            index_Coeffs_M1_Closure = i_E*N_Points_f + i_f;
            output_Opt_Coefficients << setprecision(12) << Coefficient_Matrix_Fit_Chi2[index_Coeffs_M1_Closure];
            if (index_Coeffs_M1_Closure < index_Coeffs_M1_Closure_max) {
                output_Opt_Coefficients << setw(16);
            }    
        }
    }
    output_Opt_Coefficients << endl;
}

void M1_Data_1D_Cheby :: Fit_Convergence_Check(const int &VAR_NUM) {
    long double max_err_Chi2_E, max_err_Chi2_N1_1;
    long double L2_Norm_Chi2, L_inf_Norm_Chi2;
    long double Chi2_Fit, Chi2_Numerical;
    long double ratio_E, norm_f;
    int index;
    L2_Norm_Chi2 = 0.0;
    L_inf_Norm_Chi2 = 0.0;
    
    for (int i_E = 0; i_E < N_Points_E; i_E++) {
        for (int i_f = 0; i_f < N_Points_f; i_f++) {
            index = i_E*N_Points_f + i_f;
            ratio_E = ratio_I0_NON_GRAY[index];
            norm_f = N1_1_NON_GRAY[index];
                
            Chi2_Numerical = Chi2_NON_GRAY[index];
            Chi2_Fit = Evaluate_Chi2(ratio_E, norm_f);
            
            // cout << "ratio_E = " << ratio_E << "   " << "norm_f = " << norm_f << "   " << "Chi2_Fit = " << Chi2_Fit << "    " << "Chi2_Numerical = " << Chi2_Numerical << endl;
                
            L2_Norm_Chi2 += pow(Chi2_Fit - Chi2_Numerical, 2);
            L_inf_Norm_Chi2 = max(L_inf_Norm_Chi2, fabs(Chi2_Fit - Chi2_Numerical));
                
            if (L_inf_Norm_Chi2 == fabs(Chi2_Fit - Chi2_Numerical)) {
                max_err_Chi2_E = ratio_I0_NON_GRAY[index];
                max_err_Chi2_N1_1 = N1_1_NON_GRAY[index];
            }
        }
    }
    L2_Norm_Chi2 /= N_Points_E*N_Points_f;
    L2_Norm_Chi2 = sqrt(L2_Norm_Chi2);
      
    cout << "Convergence Stats ........." << "L2_Norm_Chi2 = " << L2_Norm_Chi2 << "     "  << "L_inf_Norm_Chi2 = " << L_inf_Norm_Chi2 << "    " << "max_err_Chi2_E = " << max_err_Chi2_E << "    " << "max_err_Chi2_N1_1 = " << max_err_Chi2_N1_1 << endl;
}

void M1_Data_1D_Cheby :: Fit_Convergence_Check(M1_Data_1D_Cheby *Chi2_M1_1D_Uniform, const int &flag_Write_Output, ofstream &out_L_inf_Norm, const int &VAR_NUM) {
    long double max_err_Chi2_E, max_err_Chi2_N1_1;
    long double L2_Norm_Chi2, L_inf_Norm_Chi2;
    long double Chi2_Fit, Chi2_Numerical;
    long double ratio_E, norm_f;
    long double Mobius_Scale_Actual, Mobius_Scale_Fit;
    long double error_fit;
    long double weight_total;
    long double weight_L_Chi2_val, weight_r_I0_val, weight_N1_val;
    int index;
    L2_Norm_Chi2 = 0.0;
    L_inf_Norm_Chi2 = 0.0;
    
    long double x_L_Chi2[Chi2_M1_1D_Uniform->N_pts_Mob_Scale], weight_L_Chi2[Chi2_M1_1D_Uniform->N_pts_Mob_Scale];
    long double x_r_I0[Chi2_M1_1D_Uniform->N_Points_E], weight_r_I0[Chi2_M1_1D_Uniform->N_Points_E];
    long double x_N1[2*(Chi2_M1_1D_Uniform->N_Points_f - 1) + 1], weight_N1[2*(Chi2_M1_1D_Uniform->N_Points_f - 1) + 1];
    
    chebyshev_quadrature ( weight_r_I0, x_r_I0, Chi2_M1_1D_Uniform->N_Points_E, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
    chebyshev_quadrature ( weight_N1, x_N1, 2*(Chi2_M1_1D_Uniform->N_Points_f - 1) + 1, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
    gen_laguerre_ek_compute ( Chi2_M1_1D_Uniform->N_pts_Mob_Scale, 0.0, x_L_Chi2, weight_L_Chi2 );
    
    for (int id_Mobius = 0; id_Mobius < Chi2_M1_1D_Uniform->N_pts_Mob_Scale; id_Mobius++) {
        Mobius_Scale_Actual = x_L_Chi2[id_Mobius];
        for (int i_E = 0; i_E < Chi2_M1_1D_Uniform->N_Points_E; i_E++) {
            for (int i_f = 0; i_f < Chi2_M1_1D_Uniform->N_Points_f; i_f++) {
                index = (id_Mobius*Chi2_M1_1D_Uniform->N_Points_E + i_E)*Chi2_M1_1D_Uniform->N_Points_f + i_f;
                
                ratio_E = Chi2_M1_1D_Uniform->ratio_I0_NON_GRAY[index];
                norm_f = Chi2_M1_1D_Uniform->N1_1_NON_GRAY[index];
                
                Mobius_Scale_Fit = Evaluate_Length_Scale(norm_f);
                
                ratio_E = Recompute_I0_Mapping(ratio_E, Mobius_Scale_Actual, Mobius_Scale_Fit);
                
                Chi2_Numerical = Chi2_M1_1D_Uniform->Chi2_NON_GRAY[index];
                Chi2_Fit = Evaluate_Chi2(ratio_E, norm_f);
                
                // cout << "ratio_E = " << ratio_E << "   " << "norm_f = " << norm_f << "   " << "Chi2_Fit = " << Chi2_Fit << "    " << "Chi2_Numerical = " << Chi2_Numerical << endl;
                
                error_fit = Chi2_Fit - Chi2_Numerical;
                
                if (x_r_I0[i_E] != Chi2_M1_1D_Uniform->ratio_I0_NON_GRAY[index] && Chi2_M1_1D_Uniform->E_NON_GRAY[index] < 1.0e7) {
                    cout << "Problem with gauss-Lobatto-Chebyshev quadrature for r_I0" << endl;
                }
                
                if (x_N1[i_f + (Chi2_M1_1D_Uniform->N_Points_f - 1)] != norm_f) {
//                     cout << "Problem with gauss-Lobatto-Chebyshev quadrature for N1" << endl;
                }
                
                weight_L_Chi2_val = exp(x_L_Chi2[id_Mobius]) * weight_L_Chi2[id_Mobius];
                if (isinf(weight_L_Chi2_val) || isnan(weight_L_Chi2_val)) {
                    weight_L_Chi2_val = 0.0;
                }
                
                weight_r_I0_val = weight_r_I0[i_E];
                weight_N1_val = weight_N1[i_f + (Chi2_M1_1D_Uniform->N_Points_f - 1)];
                
                weight_total = weight_L_Chi2_val * weight_r_I0_val * weight_N1_val;
                
                L2_Norm_Chi2 += weight_total * pow(error_fit, 2);
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

void M1_Data_1D_Cheby :: Precompute_Final_Max_Ent_Solution(const int &Maximum_Entropy_Solution_Regime) {
    int Problem_Type;
    int Node_Distribution_E, Node_Distribution_f;
    int index;
    M1_Var_Num_Points Num_points;
    Algebraic_Mapping_N1 struct_map_N1;
    record_Chi2 *rec_Chi2_array;
    
    rec_Chi2_array = new record_Chi2[N_Points_E*N_Points_f];
    
    Problem_Type = NON_GRAY;
    Node_Distribution_E = UNIFORM_DISTRIBUTION; //CHEBYSHEV_FIRST_KIND_DISTRIBUTION; //CHEBYSHEV_SECOND_KIND_DISTRIBUTION; //
    Node_Distribution_f = CHEBYSHEV_FIRST_KIND_DISTRIBUTION; //CHEBYSHEV_SECOND_KIND_DISTRIBUTION; //
    
    Num_points.E = N_Points_E;
    Num_points.f = N_Points_f;
    Num_points.Maximum_Entropy_Solution_Regime = Maximum_Entropy_Solution_Regime;
    
    Num_points.Length_Scale_Dist_Type = LENGTH_SCALE_DIST_FIT;
    Num_points.Least_Squares_L_Chi2_Mode = Least_Squares_L_Chi2_Mode;
    
    // Setup parameters for algebraic mapping for the first-order normalized angular moment
    struct_map_N1.f_L_Chi2 = 0.0;
    struct_map_N1.flag_Algebraic_Map_N1 = false;
    
    OPTIM_NON_GRAY_M1_Array(rec_Chi2_array, NULL, CLOSING_FLUX, &Num_points, Problem_Type, ONE_DIMENSIONAL, Node_Distribution_E, Node_Distribution_f, 0, 1, N_pts_L_Chi2_LS, Coefficient_Matrix_Fit_Mob_Scale, &struct_map_N1, false, false);
    
    for (int index_e = 0 ; index_e < N_Points_E; index_e++) {
        for (int i_f = 0; i_f < N_Points_f; i_f++) {
            index = index_e*N_Points_f + i_f;
            
            ratio_I0_NON_GRAY[index] = rec_Chi2_array[index].ratio_I0;
            E_NON_GRAY[index] = rec_Chi2_array[index].I0;
            N1_1_NON_GRAY[index] = rec_Chi2_array[index].N1;
            Chi2_NON_GRAY[index] = rec_Chi2_array[index].Chi2;
            
            dChi2_drI0[index] = rec_Chi2_array[index].dChi2_drI0;
            dChi2_dI0[index] = rec_Chi2_array[index].dChi2_dI0;
            dChi2_dN1[index] = rec_Chi2_array[index].dChi2_dN1;
            d2Chi2_dN12[index] = rec_Chi2_array[index].d2Chi2_dN12;
            d3Chi2_dN13[index] = rec_Chi2_array[index].d3Chi2_dN13;
            d2Chi2_drI0_dN1[index] = rec_Chi2_array[index].d2_Chi2_drI0_dN1;
            d3Chi2_drI0_dN1_2[index] = rec_Chi2_array[index].d3Chi2_drI0_dN1_2;
            d2Chi2_dI0_dN1[index] = rec_Chi2_array[index].d2_Chi2_dI0_dN1;
            d3Chi2_dI0_dN1_2[index] = rec_Chi2_array[index].d3Chi2_dI0_dN1_2;
            
            if (id_proc == 0) {
                if (Chi2_NON_GRAY[index] < pow(N1_1_NON_GRAY[index], 2)) {
                    cout << "E = " << E_NON_GRAY[index] << "   "  << "N1_1 = " << N1_1_NON_GRAY[index] << "   " << "Chi2 = " << Chi2_NON_GRAY[index] << endl;
                    cout << "Non Realizable moment" << endl;
                    exit(0);
                }
            }
        }
    }
    
//     cout << "Checking hyperbolicity for Max Entropy solutions !!!!!!!!!!!!!!!" << endl;
//     Check_Hyperbolicity_Max_Ent_Solutions(0);
}

long double dmapping_L_Chi2_dnorm_f_2(const long double &L, const long double &norm_f) {
    long double dr_N1_dnorm_f_2;
    long double norm_f_2;
    
    norm_f_2 = norm_f*norm_f;
    
    dr_N1_dnorm_f_2 = 2.0*(1.0 - L)/pow((L*norm_f_2 - 1.0), 2);
    
    return dr_N1_dnorm_f_2;
}

long double dmapping_L_Chi2_d_f_L_Chi_2(const long double &L, const long double &norm_f) {
    long double dr_N1_d_f_L_Chi_2;
    long double norm_f_2;
    
    norm_f_2 = norm_f*norm_f;
    
    dr_N1_d_f_L_Chi_2 = -2.0*norm_f_2*(1.0 - norm_f_2)/pow((L*norm_f_2 - 1.0), 2);
    
    return dr_N1_d_f_L_Chi_2;
}

long double mapping_L_Chi2(const long double &L, const long double &norm_f) {
    long double r_N1;
    long double norm_f_2;
    
    norm_f_2 = norm_f*norm_f;
    
    r_N1 = 2.0*(1.0 - norm_f_2)/(L*norm_f_2 - 1.0) + 1.0;
    
    if (1.0 - L < 1.0e-8) {
        if (norm_f < 1.0e-6) {
            r_N1 = -1.0;
        } else {
            r_N1 = 1.0;
        }
    }
    
    return r_N1;
}

long double Inverse_mapping_L_Chi2(const long double &L, const long double &r_N1) {
    long double norm_f_2;
    
    norm_f_2 = (r_N1 + 1.0)/(L*r_N1 + (2.0 - L));
    
    if (1.0 - L < 1.0e-8) {
        if (r_N1 + 1.0 < 1.0e-6) {
            norm_f_2 = 0.0;
        } else {
            norm_f_2 = 1.0;
        }
    }
    
    return sqrt(norm_f_2);
}

long double M1_Data_1D_Cheby :: Recompute_I0_Mapping(const long double &ratio_E, const long double &L_N3_ijk_Actual, const long double &L_N3_ijk_Fit) {
    long double ratio_E_temp;
    
    if (fabs(ratio_E) == 1.0) {
        ratio_E_temp = ratio_E;
    } else {
        // Compute the energy density associated with the length scale associated
        // with the uniform distribution
        ratio_E_temp = Inverse_Mobius_Transformation(ratio_E, L_N3_ijk_Actual);
        
        // Compute the value of the mapping of I0 associated with the approximated
        // optimal length scale
        ratio_E_temp = Mobius_Transformation(ratio_E_temp, L_N3_ijk_Fit);
    }
    
    return ratio_E_temp;
}

long double M1_Data_1D_Cheby :: Compute_drI01_drI02(const long double &ratio_E, const long double &L_N3_ijk_Actual, const long double &L_N3_ijk_Fit) {
    long double drI01_drI02;
    long double temp_val;
    long double a_val, x_val;
    
    a_val = L_N3_ijk_Fit/L_N3_ijk_Actual;
    x_val = (1.0 - ratio_E)/2.0;
    
    drI01_drI02 = pow(x_val, a_val - 1.0);
    
    if (fabs(x_val) < 1.0e-8) {
        if (a_val > 1.0) {
            drI01_drI02 = 0.0;
        } else if (a_val < 1.0) {
            drI01_drI02 = 1.0e32; // infinity in this case
        } else {
            drI01_drI02 = 1.0;
        }
    }
    
    return drI01_drI02;
}
