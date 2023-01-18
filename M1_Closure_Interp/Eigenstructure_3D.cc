#ifndef _EIGENSTRUCTURE_3D_H_INCLUDED
#include "./Eigenstructure_3D.h"
#endif // _EIGENSTRUCTURE_3D_H_INCLUDED

long double c = 1.0;

//*************************************************************************
// This routine computes the lower bound of the realizability 
// domain for the normalized second-order angular moment
//*************************************************************************
long double Chi_2_Lower_bound_1D(const long double &N1) {
    long double b_minus;
    b_minus = N1*N1;
    return b_minus;
}

//*************************************************************************
// This routine computes the upper bound of the realizability 
// domain for the normalized second-order angular moment
//*************************************************************************
long double Chi_2_Upper_bound_1D(const long double &N1) {
    long double b_plus;
    b_plus = 1.0;
    return b_plus;
}

//*************************************************************************
// This routine check realizability of our approximation of the non-gray M1
// closure
//*************************************************************************
void Check_Realizability(const long double &N1, const long double &Chi2) {
    long double b_minus, b_plus, Mobius_Scale_Fit, E;
    
    b_minus = Chi_2_Lower_bound_1D(N1);
    b_plus = Chi_2_Upper_bound_1D(N1);
    
    // cout << "N1 = " << N1 << "   " << "b_minus = " << b_minus << "   " << "Chi2 = " << Chi2 << "   " << "b_plus = " << b_plus << endl;
    
    // Now assess realizability of our approximation of the Eddington factor
    if (Chi2 < b_minus) {
        cout << "Non-Realizable closure with Chi2 < b_minus" << "   " << "N1 = " << N1 << "   " << "Chi2 = " << Chi2 << "   " << "b_minus = " << b_minus << "   " << "Chi2 - b_minus = " << Chi2 - b_minus << endl;
        exit(0);
    } else if (Chi2 > b_plus) {
        cout << "Non-Realizable closure with Chi2 > b_plus" << "   " << "N1 = " << N1 << "   " << "Chi2 = " << Chi2 << "   " << "b_plus = " << b_plus << "   " << "Chi2 - b_plus = " << Chi2 - b_plus << endl;
        exit(0);
    }
}

long double unit_vec_N1(const int &i, const long double &N1, const long double &phi) {
    long double n_i;
    
    switch (i) {
        case 1:
            n_i = N1*cos(phi);
            break;
        case 2:
            n_i = N1*sin(phi);
            break;
    };
    n_i /= N1;
    return n_i;
}

long double I2ij(const int &i, const int &j, const long double &I0, const long double &N1, const long double &phi, const long double &Chi2) {
    long double I2_ij, N2_ij, chi_T;
    
    N2_ij = (1.0/2.0)*(3.0*Chi2 - 1.0)*unit_vec_N1(i, N1, phi)*unit_vec_N1(j, N1, phi);
    
    if (i == j){
        N2_ij += (1.0/2.0)*(1.0 - Chi2);
    }
    
    I2_ij = N2_ij*I0;
    
    // cout << "I2_ij = " << I2_ij << endl;
    
    return I2_ij;
}

long double dI2ij_dI0(const int &i, const int &j, const long double &r_E, const long double &N1, const long double &phi, const long double &Chi2, const long double &E_dChi_dE) {
    long double dN2_ij, chi_T;
    
    chi_T = Chi2 + E_dChi_dE;
    dN2_ij = (1.0/2.0)*(3.0*chi_T - 1.0)*unit_vec_N1(i, N1, phi)*unit_vec_N1(j, N1, phi);
    
    if (i == j){
        dN2_ij += (1.0/2.0)*(1.0 - chi_T);
    }
    
    // cout << "dI2ij_dI0 = " << dN2_ij << endl;
    
    return dN2_ij;
}
  
long double dI2ij_dIi(const int &i, const int &j, const int &index_N1, const long double &r_E, const long double &N1, const long double &phi, const long double &Chi2, const long double &dChi_dN1) {
    long double dN2_ij, N1_i;
    
    switch (index_N1) {
        case 1:
            N1_i = N1*cos(phi);
            break;
        case 2:
            N1_i = N1*sin(phi);
            break;
    };
    
    dN2_ij = (1.0/2.0)*(3.0*N1_i*dChi_dN1 - 2.0*(3.0*Chi2-1.0)*unit_vec_N1(index_N1,N1,phi))*unit_vec_N1(i,N1,phi)*unit_vec_N1(j,N1,phi);
    
    if (i == j){
        dN2_ij -= (1.0/2.0)*N1_i*dChi_dN1;
    }
    
    if (index_N1 == i) {
        dN2_ij += (1.0/2.0)*(3.0*Chi2-1.0)*unit_vec_N1(j,N1,phi);
    }
    
    if (index_N1 == j) {
        dN2_ij += (1.0/2.0)*(3.0*Chi2-1.0)*unit_vec_N1(i,N1,phi);
    }
    
    dN2_ij /= N1;
//     cout << "dI2ij_dI1 = " << dN2_ij << endl;
    return dN2_ij;
}

void Cubic_solver(const long double &aa, const long double &bb, const long double &cc, const long double &dd, long double &x1, long double &x2, long double &x3, bool flag_Imaginary_Part) {
     double TOLERANCE = 1.0e-10;
    // Solve: aa X^3 + bb X^2 + cc X + dd = 0
    
     // VARIABLE DECLARATION
     // ----------------------
     
     long double delta, delta0, delta1;
     std::complex<long double> C0, CC, CTEMP;
     long double x1temp, x2temp, x3temp;
     long double temp;
     std::complex<long double>  temp_Val;
     
     //  CUBIC SOLVER
     // ----------------------
     // Intermediate coefficients
     
     // Discriminant of the general cubic function
     delta = 18.0*aa*bb*cc*dd - 4.0*bb*bb*bb*dd + bb*bb*cc*cc - 4.0*aa*cc*cc*cc - 27.0*aa*aa*dd*dd;
     
     if (delta < -TOLERANCE) {
         cout << "Determinant of characteristic polynomial in 2D is negative : delta = " << delta << endl; 
         cout << "Complex eigenvalues might be encountered !!" << endl; 
         // exit(0);
     } else if (fabs(delta) < TOLERANCE) {
         delta = 0.0;
     }
     
     delta0 = bb*bb - 3.0*aa*cc;
     delta1 = 2.0*bb*bb*bb - 9.0*aa*bb*cc + 27.0*aa*aa*dd;
     temp = delta1*delta1 - 4.0*delta0*delta0*delta0;
     
     if (fabs(delta0) < 1.0e-8 && fabs(delta1) < 1.0e-8) {
         // Then we have a triple root
         x1temp = -bb/(3.0*aa);
         x2temp = x1temp;
         x3temp = x1temp;
     } else {
         if (temp < 0.0) {
             C0 = (delta1 + sqrt(-temp)*1.0i)/2.0;
        } else {
            C0 = (delta1 + sqrt(temp))/2.0;
        }
        C0 = pow(C0, 1.0/3.0);
        
        // Cubic roots
        CC = C0; 
        temp_Val = bb + CC + delta0/CC;
        
        if (!flag_Imaginary_Part) {
            // Then return real part
            x1temp = - real(temp_Val)/(3.0*aa);
            
//             if(fabs(imag(temp_Val)) > 1.0e-8) {
//                 cout << "imag(temp_Val) = " << imag(temp_Val) << "  " << "real(temp_Val) = " << real(temp_Val) << "  " << "delta = " << delta << endl;
//             }
        } else {
            // Then return imaginary part
            x1temp = - imag(temp_Val)/(3.0*aa);
        }
        
        if (x1temp != x1temp) {
            cout << "x1_temp is nan" << "  " << "temp_Val = " << temp_Val << "  " << "aa = " << aa << "  " << "delta0 = " << delta0 << "  " << "CC = " << CC << endl;
            exit(0);
        }
        
        CTEMP = -0.5 - 0.5*sqrt(3.0)*1.0i;
        CC = C0*CTEMP; 
        temp_Val = bb + CC + delta0/CC;
        
        if (!flag_Imaginary_Part) {
            // Then return real part
            x2temp = - real(temp_Val)/(3.0*aa);
            
//             if(fabs(imag(temp_Val)) > 1.0e-8) {
//                 cout << "imag(temp_Val) = " << imag(temp_Val) << "  " << "real(temp_Val) = " << real(temp_Val) << "  " << "delta = " << delta << endl;
//             }
        } else {
            // Then return imaginary part
            x2temp = - imag(temp_Val)/(3.0*aa);
        }
        
        if (x2temp != x2temp) {
            cout << "x2_temp is nan" << endl;
            exit(0);
        }
        
        CTEMP = -0.5 + 0.5*sqrt(3.0)*1.0i;
        CC = C0*CTEMP; 
        temp_Val = bb + CC + delta0/CC;
        
        if (!flag_Imaginary_Part) {
            // Then return real part
            x3temp = - real(temp_Val)/(3.0*aa);
            
//             if(fabs(imag(temp_Val)) > 1.0e-8) {
//                 cout << "imag(temp_Val) = " << imag(temp_Val) << "  " << "real(temp_Val) = " << real(temp_Val) << "  " << "delta = " << delta << endl;
//             }
        } else {
            // Then return imaginary part
            x3temp = - imag(temp_Val)/(3.0*aa);
        }
        
        if (x3temp != x3temp) {
            cout << "x3_temp is nan" << endl;
            exit(0);
        }
     }
     
     // Reorganizing eigenvalues (increasing order)
     if (x1temp<x2temp){
         if (x2temp<x3temp){
             x1 = x1temp; x2 = x2temp; x3 = x3temp;
        } else {
            x1 = x1temp; x2 = x3temp; x3 = x2temp;
        }
    } else if (x2temp<x3temp) {
        if(x1temp<x3temp){
            x1 = x2temp; x2 = x1temp; x3 = x3temp;
        } else {
            x1 = x2temp; x2 = x3temp; x3 = x1temp;
        }
    } else {
        if(x1temp<x2temp) {
            x1 = x3temp; x2 = x1temp; x3 = x2temp;
        } else {
            x1 = x3temp; x2 = x2temp; x3 = x1temp;
        }
    }
    
    // cout << "aa = " << aa << "  " << "bb = " << bb << "   " << "cc = " << cc << "   " << "dd = " << dd << endl;
    
    // cout << "root 1 = " << aa*pow(x1,3) + bb*pow(x1,2) + cc*x1 + dd << endl;
    // cout << "root 2 = " << aa*pow(x2,3) + bb*pow(x2,2) + cc*x2 + dd << endl;
    // cout << "root 3 = " << aa*pow(x3,3) + bb*pow(x3,2) + cc*x3 + dd << endl;
}

void Evaluate_eigenvalues_1D(const long double &r, const long double &N1, const long double &Chi2, const long double &E_dChi2_dE, const long double &dChi2_dN1, long double &lambda_1, long double &lambda_2) {
    long double temp_val, temp_val_temp;
    long double TOLERANCE = 1.0e-10;
    
    temp_val_temp = pow(dChi2_dN1, 2) + 4.0*(Chi2 + E_dChi2_dE);
    if (temp_val_temp < -TOLERANCE) {
        cout << "temp_val is negative. Value is : " << temp_val_temp << endl;
        exit(0);
    } else if (fabs(temp_val_temp) < TOLERANCE) {
        temp_val_temp = 0.0;
    }
    
    temp_val = sqrt(temp_val_temp);
    
    lambda_1 = (c/2.0)*(dChi2_dN1 + temp_val);
    lambda_2 = (c/2.0)*(dChi2_dN1 - temp_val);
    
    if (temp_val_temp < -TOLERANCE) {
        cout << "r = " << r << "    " << "N1 = " << N1 << "\n";
        cout << "Chi2 = " << Chi2 << "\n";
        cout << "E_dChi2_dE = " << E_dChi2_dE << "    " << "dChi2_dN1 = " << dChi2_dN1 << "    " << "Chi2 + E_dChi2_dE = " << Chi2 + E_dChi2_dE << "\n";
        cout << "temp_val_temp = " << temp_val_temp << "  "  << "temp_val = " << temp_val << "\n";
        cout << "lambda_1 = " << lambda_1 << "    " << "lambda_2 = " << lambda_2 << "\n";
        cout << "\n";
        exit(0);
    }
    
//     if (fabs(lambda_1) > 1.0 + 1.0e-2 || fabs(lambda_2) > 1.0 + 1.0e-2) {
//         cout << "r = " << r << "    " << "N1 = " << N1 << "\n";
//         cout << "Chi2 = " << Chi2 << "\n";
//         cout << "E_dChi2_dE = " << E_dChi2_dE << "    " << "dChi2_dN1 = " << dChi2_dN1 << "    " << "Chi2 + E_dChi2_dE = " << Chi2 + E_dChi2_dE << "\n";
//         cout << "temp_val = " << temp_val << "\n";
//         cout << "lambda_1 = " << lambda_1 << "    " << "lambda_2 = " << lambda_2 << "\n";
//         cout << "diff lambda_1 = " << 1.0 - fabs(lambda_1) << "    " << "diff lambda_2 = " << 1.0 - fabs(lambda_2) << "\n";
//         cout << "\n";
//         exit(0);
//     }
}

long double xi(const long double &N1) {
    long double xi_val;
    xi_val = sqrt(4.0-3.0*pow(N1,2));
    return xi_val;
}

long double Chi_2_Gray(const long double &N1) {
    long double Chi_2_Gray_val, xi_val;
    xi_val = xi(N1);
    Chi_2_Gray_val = (3.0 + 4.0*pow(N1, 2))/(5.0 + 2.0*xi_val);
    return Chi_2_Gray_val;
}

void Evaluate_Chi2_derivatives_Gray(const long double &r, const long double &N1, long double &E_dChi2_dE, long double &dChi2_df) {
    dChi2_df = 2.0 * N1 / xi(N1);
    E_dChi2_dE = -N1*dChi2_df;
}

long double dChi2_dN1_Gray(const long double &N1) {
    long double dChi2_df;
    dChi2_df = 2.0 * N1 / xi(N1);
    return dChi2_df;
}

long double d2Chi2_dN1_2_Gray(const long double &N1) {
    long double d2_Chi2_dN1_2;
    d2_Chi2_dN1_2 = 8.0 / pow(xi(N1), 3.0/2.0);
    return d2_Chi2_dN1_2;
}

long double d3Chi2_dN1_3_Gray(const long double &N1) {
    long double d3_Chi2_dN1_3;
    d3_Chi2_dN1_3 = 72.0 * N1 / pow(xi(N1), 5.0/2.0);
    return d3_Chi2_dN1_3;
}

void Eigenvalues_Gray_M1_Closure_2D(long double &lambda_1, long double &lambda_2, long double &lambda_3, const long double &N1, const long double &phi) {
    long double norm_f_2 = pow(N1, 2);
    long double zeta = sqrt(4.0 - 3.0*norm_f_2);
    long double N1_1, N1_2;
    long double temp_val;
    
    N1_1 = N1*cos(phi);
    N1_2 = N1*sin(phi);
    
    temp_val = 2.0*(zeta-1.0)*(zeta+2.0)*(2.0*(zeta-1.0)*(zeta+2.0) + 3.0*pow(N1_2, 2));
    temp_val = sqrt(temp_val)/(sqrt(3.0)*zeta*(zeta+2.0));
    
    lambda_1 = (N1_1/zeta) - temp_val;
    lambda_2 = (N1_1/zeta) + temp_val;
    lambda_3 = (N1_1/norm_f_2)*(2.0 - zeta);
}

void Evaluate_Lambdas_x(long double &lambda_1, long double &lambda_2, long double &lambda_3, const long double &r_E, const long double &N1, const long double &phi, const long double &Chi2, const long double &E_dChi2_dE, const long double &dChi2_df, bool flag_Imaginary_Part) {
    long double aa, bb, cc, dd;
    long double b, b1, b2;
    long double dFdU_12, dFdU_21, dFdU_22, dFdU_23, dFdU_31, dFdU_32, dFdU_33;
    long double w1, w2;
    long double temp_val_1, temp_val_2, temp_val_3;
    
    if (fabs(N1) < TOLER_HYPER){
        if (!flag_Imaginary_Part) {
            lambda_1 = 0.0;
            lambda_2 = -c/sqrt(3.0);
            lambda_3 = c/sqrt(3.0);
        } else {
            lambda_1 = 0.0;
            lambda_2 = 0.0;
            lambda_3 = 0.0;
        }
    } else {
        // cout << "E_dChi2_dE = " << E_dChi2_dE << "    " << "dChi2_df = " << dChi2_df << "\n";
        
        dFdU_12 = c;
        dFdU_21 = c*dI2ij_dI0(1, 1, r_E, N1, phi, Chi2, E_dChi2_dE);
        dFdU_22 = c*dI2ij_dIi(1, 1, 1, r_E, N1, phi, Chi2, dChi2_df);
        dFdU_23 = c*dI2ij_dIi(1, 1, 2, r_E, N1, phi, Chi2, dChi2_df);
        dFdU_31 = c*dI2ij_dI0(1, 2, r_E, N1, phi, Chi2, E_dChi2_dE);
        dFdU_32 = c*dI2ij_dIi(1, 2, 1, r_E, N1, phi, Chi2, dChi2_df);
        dFdU_33 = c*dI2ij_dIi(1, 2, 2, r_E, N1, phi, Chi2, dChi2_df);
        
        // cout << endl;
        // cout << 0.0 << "    " << 1.0 << "    " << 0.0 << "\n";
        // cout << setprecision(12) << dFdU_21 << "    " << dFdU_22 << "    " << dFdU_23 << "\n";
        // cout << setprecision(12) << dFdU_31 << "    " << dFdU_32 << "    " << dFdU_33 << "\n";
        // cout << endl;
        
        aa = -1.0;
        bb = dFdU_22 + dFdU_33;
        cc = dFdU_12*dFdU_21 - dFdU_22*dFdU_33 + dFdU_23*dFdU_32;
        dd = -dFdU_12*dFdU_21*dFdU_33 + dFdU_12*dFdU_23*dFdU_31;
        
        Cubic_solver(aa, bb, cc, dd, lambda_1, lambda_2, lambda_3, flag_Imaginary_Part);
        
        // cout << setprecision(12) << "lambda_1 = " << lambda_1 << "    " << "lambda_2 = " << lambda_2 << "    " << "lambda_3 = " << lambda_3 << "\n";
    }
    
    if (1.0 - fabs(lambda_1) < -1.0e-2 || 1.0 - fabs(lambda_2) < -1.0e-2 || 1.0 - fabs(lambda_3) < -1.0e-2) {
        cout << setprecision(12) << "lambda_1 = " << lambda_1 << "    " << "lambda_2 = " << lambda_2 << "    " << "lambda_3 = " << lambda_3 << "\n";
        cout << "dFdU_21 = " << dFdU_21 << "    " << "dFdU_22 = " << dFdU_22 << "    " << "dFdU_23 = " << dFdU_23 << "\n";
        cout << "dFdU_31 = " << dFdU_31 << "    " << "dFdU_32 = " << dFdU_32 << "    " << "dFdU_33 = " << dFdU_33 << "\n";
        exit(0);
    }
}

void Setup_rc(long double rc_vec[][3], const int &index_U, const long double &r_E, const long double &N1, const long double &phi, const long double &Chi2, const long double &E_dChi2_dE, const long double &dChi2_df) {
    long double lam_val, lam_val_normal, temp_val;
    long double lambda_1, lambda_2, lambda_3;
    long double dFdU_21, dFdU_22, dFdU_23;
    long double dFdU_31, dFdU_32, dFdU_33;
    long double temp_val_check;
    
    Evaluate_Lambdas_x(lambda_1, lambda_2, lambda_3, r_E, N1, phi, Chi2, E_dChi2_dE, dChi2_df);
    
    // cout << "lambda_1 1D = " << lambda_1 <<  "   " << "lambda_2 1D = " << lambda_2 <<  "   " << "lambda_3 = " << lambda_3 << endl;
    
    if (fabs(N1) < TOLER_HYPER) {
        // This is the case where we are in the isotropic regime
        switch (index_U) {
            case 1 :
                rc_vec[0][index_U-1] = 0.0; 
                rc_vec[1][index_U-1] = 0.0;
                rc_vec[2][index_U-1] = 1.0;
                break;
            case 2 :
                rc_vec[0][index_U-1] = -sqrt(3.0);
                rc_vec[1][index_U-1] = 1.0;
                rc_vec[2][index_U-1] = 0.0;
                break;
            case 3 :
                rc_vec[0][index_U-1] = sqrt(3.0);
                rc_vec[1][index_U-1] = 1.0;
                rc_vec[2][index_U-1] = 0.0;
                break;
            default:
                cout << "Incorrect value for index_U = " << index_U << endl;
                exit(0);
                break;
        }
    } else if (fabs(N1*sin(phi)) < TOLER_HYPER) {
        // This is the case where we are in one dimensional radiative transfer
        // In 1D, we have: A_23 = 0
        switch (index_U) {
            case 1 :
                lam_val = lambda_1;
                break;
            case 2 :
                lam_val = lambda_2;
                break;
            case 3 :
                lam_val = lambda_3;
                break;
            default:
                cout << "Incorrect value for index_U = " << index_U << endl;
                exit(0);
                break;
        }
        
        dFdU_21 = dI2ij_dI0(1,1, r_E, N1, phi, Chi2, E_dChi2_dE);
        dFdU_22 = dI2ij_dIi(1,1,1, r_E, N1, phi, Chi2, dChi2_df);
        dFdU_23 = dI2ij_dIi(1,1,2, r_E, N1, phi, Chi2, dChi2_df);
        dFdU_31 = dI2ij_dI0(1,2, r_E, N1, phi, Chi2, E_dChi2_dE);
        dFdU_32 = dI2ij_dIi(1,2,1, r_E, N1, phi, Chi2, dChi2_df);
        dFdU_33 = dI2ij_dIi(1,2,2, r_E, N1, phi, Chi2, dChi2_df);
        
        lam_val_normal = lam_val/c;
        
        temp_val_check = c + lam_val*(dFdU_22 - lam_val)/dFdU_21;
        
        if (fabs(lambda_1 - lambda_2) < 1.0e-8 && fabs(lambda_1 - lambda_3) < 1.0e-8) {
            // In this case we have three degenerate eigenvalues
            switch (index_U) {
                case 1 :
                    rc_vec[0][index_U-1] = 1.0;
                    rc_vec[1][index_U-1] = lam_val_normal;
                    rc_vec[2][index_U-1] = 0.0;
                    break;
                case 2 :
                    rc_vec[0][index_U-1] = lam_val_normal;
                    rc_vec[1][index_U-1] = -1.0;
                    rc_vec[2][index_U-1] = 0.0;
                    break;
                case 3 :
                    rc_vec[0][index_U-1] = 0.0;
                    rc_vec[1][index_U-1] = 0.0;
                    rc_vec[2][index_U-1] = 1.0;
                    break;
                default:
                    cout << "Incorrect value for index_U = " << index_U << endl;
                    exit(0);
                    break;
            }
        } else {
            if (fabs(temp_val_check) < 1.0e-8) {
                rc_vec[0][index_U-1] = 1.0;
                rc_vec[1][index_U-1] = lam_val_normal;
                rc_vec[2][index_U-1] = 0.0;
            } else {
                rc_vec[0][index_U-1] = 0.0;
                rc_vec[1][index_U-1] = 0.0;
                rc_vec[2][index_U-1] = 1.0;
            }
        }
    } else {
        // Then we are in 2D
        switch (index_U) {
            case 1 :
                lam_val = lambda_1;
                break;
            case 2 :
                lam_val = lambda_2;
                break;
            case 3 :
                lam_val = lambda_3;
                break;
            default:
                cout << "Incorrect value for index_U = " << index_U << endl;
                exit(0);
                break;
        }
        
        dFdU_21 = dI2ij_dI0(1,1, r_E, N1, phi, Chi2, E_dChi2_dE);
        dFdU_22 = dI2ij_dIi(1,1,1, r_E, N1, phi, Chi2, dChi2_df);
        dFdU_23 = dI2ij_dIi(1,1,2, r_E, N1, phi, Chi2, dChi2_df);
                    
        lam_val_normal = lam_val/c;
        
        // -(A21 + (lam/c) (A22 - lam))/A23
        temp_val = -(dFdU_21 + lam_val_normal*(dFdU_22 - lam_val));
        temp_val /= dFdU_23;
                    
        dFdU_31 = dI2ij_dI0(1,2, r_E, N1, phi, Chi2, E_dChi2_dE);
        dFdU_32 = dI2ij_dIi(1,2,1, r_E, N1, phi, Chi2, dChi2_df);
        dFdU_33 = dI2ij_dIi(1,2,2, r_E, N1, phi, Chi2, dChi2_df);
        
        temp_val_check = dFdU_23 - (dFdU_33 - lam_val)*(dFdU_21 + lam_val_normal*(dFdU_22 - lam_val))/(dFdU_31 + lam_val_normal*dFdU_32);
        
        if (fabs(temp_val_check) > TOLER_HYPER) {
            cout << "temp_val_check = " << temp_val_check <<  "   " << "lam_val_normal = " << lam_val_normal << endl;
        }
        
        rc_vec[0][index_U-1] = 1.0;
        rc_vec[1][index_U-1] = lam_val_normal;
        rc_vec[2][index_U-1] = temp_val;
    }
}

void Setup_lc(long double lc_vec[][3], const int &index_U, const long double &r_E, const long double &N1, const long double &phi, const long double &Chi2, const long double &E_dChi2_dE, const long double &dChi2_df) {
    long double lambda_1, lambda_2, lambda_3;
    long double dFdU_21, dFdU_22, dFdU_23;
    long double dFdU_31, dFdU_32, dFdU_33;
    long double diff_lam_val, lam_val, temp_val;
    long double deter, temp_val_check_lam1, temp_val_check_lam2, temp_val_check_lam3;
    
    Evaluate_Lambdas_x(lambda_1, lambda_2, lambda_3, r_E, N1, phi, Chi2, E_dChi2_dE, dChi2_df);
    
    if (fabs(N1) < TOLER_HYPER) {
        switch (index_U) {
            case 1 :
                lc_vec[index_U-1][0] = 0.0;
                lc_vec[index_U-1][1] = 0.0;
                lc_vec[index_U-1][2] = 1.0;
                break;
            case 2 :
                lc_vec[index_U-1][0] = -sqrt(3.0)/6.0;
                lc_vec[index_U-1][1] = 1.0/2.0;
                lc_vec[index_U-1][2] = 0.0;
                break;
            case 3 :
                lc_vec[index_U-1][0] = sqrt(3.0)/6.0;
                lc_vec[index_U-1][1] = 1.0/2.0;
                lc_vec[index_U-1][2] = 0.0;
                break;
            default:
                cout << "Incorrect value for index_U = " << index_U << endl;
                exit(0);
                break;
        }
    } else if (fabs(N1*sin(phi)) < TOLER_HYPER) {
        // This is the case where we are in one dimensional radiative transfer
        dFdU_21 = dI2ij_dI0(1,1, r_E, N1, phi, Chi2, E_dChi2_dE);
        dFdU_22 = dI2ij_dIi(1,1,1, r_E, N1, phi, Chi2, dChi2_df);
        dFdU_23 = dI2ij_dIi(1,1,2, r_E, N1, phi, Chi2, dChi2_df);
        // dFdU_31 = dI2ij_dI0(1,2, r_E, N1, phi, Chi2, E_dChi2_dE);
        // dFdU_32 = dI2ij_dIi(1,2,1, r_E, N1, phi, Chi2, dChi2_df);
        // dFdU_33 = dI2ij_dIi(1,2,2, r_E, N1, phi, Chi2, dChi2_df);
        
        temp_val_check_lam1 = c + lambda_1*(dFdU_22 - lambda_1)/dFdU_21;
        temp_val_check_lam2 = c + lambda_2*(dFdU_22 - lambda_2)/dFdU_21;
        temp_val_check_lam3 = c + lambda_3*(dFdU_22 - lambda_3)/dFdU_21;
        
        // cout << "temp_val_check_lam1 = " << temp_val_check_lam1 << "  " << "temp_val_check_lam2 = " << temp_val_check_lam2 << "  " << "temp_val_check_lam3 = " << temp_val_check_lam3 << endl; 
        // cout << "lambda_1 = " << lambda_1 << "  " << "lambda_2 = " << lambda_2 << "  " << "lambda_3 = " << lambda_3 << endl; 
        
        if (fabs(lambda_1 - lambda_2) < 1.0e-8 && fabs(lambda_1 - lambda_3) < 1.0e-8) {
            // In this case we have three degenerate eigenvalues
            temp_val = pow(c,2) + pow(lambda_1, 2);
            switch (index_U) {
                case 1 :
                    lc_vec[index_U-1][0] = pow(c,2)/temp_val;
                    lc_vec[index_U-1][1] = c*lambda_1/temp_val;
                    lc_vec[index_U-1][2] = 0.0;
                    break;
                case 2 :
                    lc_vec[index_U-1][0] = c*lambda_1/temp_val;
                    lc_vec[index_U-1][1] = -pow(c,2)/temp_val;
                    lc_vec[index_U-1][2] = 0.0;
                    break;
                case 3 :
                    lc_vec[index_U-1][0] = 0.0;
                    lc_vec[index_U-1][1] = 0.0;
                    lc_vec[index_U-1][2] = 1.0;
                    break;
                default:
                    cout << "Incorrect value for index_U = " << index_U << endl;
                    exit(0);
                    break;
            }
        } else {
            if (fabs(temp_val_check_lam1) > 1.0e-8) {
                deter = -(lambda_2 - lambda_3);
                switch (index_U) {
                    case 1 :
                        lc_vec[index_U-1][0] = 0.0;
                        lc_vec[index_U-1][1] = 0.0;
                        lc_vec[index_U-1][2] = lambda_3 - lambda_2;
                        break;
                    case 2 :
                        lc_vec[index_U-1][0] = lambda_3;
                        lc_vec[index_U-1][1] = -c;
                        lc_vec[index_U-1][2] = 0.0;
                        break;
                    case 3 :
                        lc_vec[index_U-1][0] = -lambda_2;
                        lc_vec[index_U-1][1] = c;
                        lc_vec[index_U-1][2] = 0.0;
                        break;
                    default:
                        cout << "Incorrect value for index_U = " << index_U << endl;
                        exit(0);
                        break;
                }
            } else if (fabs(temp_val_check_lam2) > 1.0e-8) {
                deter = lambda_1 - lambda_3;
                switch (index_U) {
                    case 1 :
                        lc_vec[index_U-1][0] = -lambda_3;
                        lc_vec[index_U-1][1] = c;
                        lc_vec[index_U-1][2] = 0.0;
                        break;
                    case 2 :
                        lc_vec[index_U-1][0] = 0.0;
                        lc_vec[index_U-1][1] = 0.0;
                        lc_vec[index_U-1][2] = lambda_1 - lambda_3;
                        break;
                    case 3 :
                        lc_vec[index_U-1][0] = lambda_1;
                        lc_vec[index_U-1][1] = -c;
                        lc_vec[index_U-1][2] = 0.0;
                        break;
                    default:
                        cout << "Incorrect value for index_U = " << index_U << endl;
                        exit(0);
                        break;
                }
            } else if (fabs(temp_val_check_lam3) > 1.0e-8) {
                deter = -(lambda_1 - lambda_2);
                switch (index_U) {
                    case 1 :
                        lc_vec[index_U-1][0] = lambda_2;
                        lc_vec[index_U-1][1] = -c;
                        lc_vec[index_U-1][2] = 0.0;
                        break;
                    case 2 :
                        lc_vec[index_U-1][0] = -lambda_1;
                        lc_vec[index_U-1][1] = c;
                        lc_vec[index_U-1][2] = 0.0;
                        break;
                    case 3 :
                        lc_vec[index_U-1][0] = 0.0;
                        lc_vec[index_U-1][1] = 0.0;
                        lc_vec[index_U-1][2] = lambda_2 - lambda_1;
                        break;
                    default:
                        cout << "Incorrect value for index_U = " << index_U << endl;
                        exit(0);
                        break;
                }
            } else {
                cout << "Error in eigenstructure analysis in 1D!!!!" << endl;
                exit(0);
            }
            
            lc_vec[index_U-1][0] /= deter;
            lc_vec[index_U-1][1] /= deter;
            lc_vec[index_U-1][2] /= deter;
        }
        
    } else {
        // Then we are in 2D
        
        dFdU_21 = dI2ij_dI0(1,1, r_E, N1, phi, Chi2, E_dChi2_dE);
        dFdU_22 = dI2ij_dIi(1,1,1, r_E, N1, phi, Chi2, dChi2_df);
        dFdU_23 = dI2ij_dIi(1,1,2, r_E, N1, phi, Chi2, dChi2_df);
        dFdU_31 = dI2ij_dI0(1,2, r_E, N1, phi, Chi2, E_dChi2_dE);
        dFdU_32 = dI2ij_dIi(1,2,1, r_E, N1, phi, Chi2, dChi2_df);
        dFdU_33 = dI2ij_dIi(1,2,2, r_E, N1, phi, Chi2, dChi2_df);
        
        deter = -(lambda_1 - lambda_2)*(lambda_1 - lambda_3)*(lambda_2 - lambda_3);
        
        switch (index_U) {
            case 1 :
                diff_lam_val = lambda_2 - lambda_3;
                lc_vec[index_U-1][0] = -diff_lam_val*(c*dFdU_21 + lambda_2*lambda_3);
                lc_vec[index_U-1][1] = c*diff_lam_val*(lambda_2 + lambda_3 - dFdU_22);
                lc_vec[index_U-1][2] = -c*dFdU_23*diff_lam_val;
                break;
            case 2 :
                diff_lam_val = lambda_1 - lambda_3;
                lc_vec[index_U-1][0] = diff_lam_val*(c*dFdU_21 + lambda_1*lambda_3);
                lc_vec[index_U-1][1] = -c*diff_lam_val*(lambda_1 + lambda_3 - dFdU_22);
                lc_vec[index_U-1][2] = c*dFdU_23*diff_lam_val;
                break;
            case 3 :
                diff_lam_val = lambda_1 - lambda_2;
                lc_vec[index_U-1][0] = -diff_lam_val*(c*dFdU_21 + lambda_1*lambda_2);
                lc_vec[index_U-1][1] = c*diff_lam_val*(lambda_1 + lambda_2 - dFdU_22);
                lc_vec[index_U-1][2] = -c*dFdU_23*diff_lam_val;
                break;
            default:
                cout << "Incorrect value for index_U = " << index_U << endl;
                exit(0);
                break;
        }
        
        lc_vec[index_U-1][0] /= deter;
        lc_vec[index_U-1][1] /= deter;
        lc_vec[index_U-1][2] /= deter;
    }
}

void Eigenstructure_Check(const long double &r_E, const long double &N1, const long double &phi, const long double &Chi2, const long double &E_dChi2_dE, const long double &dChi2_df, const int &Problem_Type) {
    static long double rc_vec[3][3], lc_vec[3][3], lambdas[3][3], dFdU[3][3];
    long double temp_val;
    long double lambda_1, lambda_2, lambda_3;
    long double dFdU_21, dFdU_22, dFdU_23;
    long double dFdU_31, dFdU_32, dFdU_33;
    
    // Compute eigenvalues
    Evaluate_Lambdas_x(lambda_1, lambda_2, lambda_3, r_E, N1, phi, Chi2, E_dChi2_dE, dChi2_df);
    
    long double lambda_1_gray, lambda_2_gray, lambda_3_gray;
    if (Problem_Type == GRAY) {
        Eigenvalues_Gray_M1_Closure_2D(lambda_1_gray, lambda_2_gray, lambda_3_gray, N1, phi);
        
        // cout << "lambda_1 = " << lambda_1 << "  " << "lambda_2 = " << lambda_2 << "  " << "lambda_3 = " << lambda_3 << "  " << endl;
        
        // cout << "lambda_1_gray = " << lambda_1_gray << "  " << "lambda_2_gray = " << lambda_2_gray << "  " << "lambda_3_gray = " << lambda_3_gray << "  " << endl;
        
        if (fabs(lambda_1 - lambda_1_gray) > 1.0e-4 && fabs(lambda_1 - lambda_2_gray) > 1.0e-4 && fabs(lambda_1 - lambda_3_gray) > 1.0e-4) {
            cout << "Eigenvalues for gray M1 closure not correct" << endl;
            cout << "lambda_1 = " << lambda_1 << "  " << "lambda_1_gray = " << lambda_1_gray << "  " << "lambda_2_gray = " << lambda_2_gray << "  " << "lambda_3_gray = " << lambda_3_gray << "  " << endl;
            exit(0);
        }
        
        if (fabs(lambda_2 - lambda_1_gray) > 1.0e-4 && fabs(lambda_2 - lambda_2_gray) > 1.0e-4 && fabs(lambda_2 - lambda_3_gray) > 1.0e-4) {
            cout << "Eigenvalues for gray M1 closure not correct" << endl;
            cout << "lambda_2 = " << lambda_2 << "  " << "lambda_1_gray = " << lambda_1_gray << "  " << "lambda_2_gray = " << lambda_2_gray << "  " << "lambda_3_gray = " << lambda_3_gray << "  " << endl;
            exit(0);
        }
        
        if (fabs(lambda_3 - lambda_1_gray) > 1.0e-4 && fabs(lambda_3 - lambda_2_gray) > 1.0e-4 && fabs(lambda_3 - lambda_3_gray) > 1.0e-4) {
            cout << "Eigenvalues for gray M1 closure not correct" << endl;
            cout << "lambda_3 = " << lambda_3 << "  " << "lambda_1_gray = " << lambda_1_gray << "  " << "lambda_2_gray = " << lambda_2_gray << "  " << "lambda_3_gray = " << lambda_3_gray << "  " << endl;
            exit(0);
        }
    }
    
    // Compute and setup flux Jacobian matrix
    dFdU_21 = dI2ij_dI0(1,1, r_E, N1, phi, Chi2, E_dChi2_dE);
    dFdU_22 = dI2ij_dIi(1,1,1, r_E, N1, phi, Chi2, dChi2_df);
    dFdU_23 = dI2ij_dIi(1,1,2, r_E, N1, phi, Chi2, dChi2_df);
    dFdU_31 = dI2ij_dI0(1,2, r_E, N1, phi, Chi2, E_dChi2_dE);
    dFdU_32 = dI2ij_dIi(1,2,1, r_E, N1, phi, Chi2, dChi2_df);
    dFdU_33 = dI2ij_dIi(1,2,2, r_E, N1, phi, Chi2, dChi2_df);
    
//     cout << endl;
//     cout << 0.0 << "    " << 1.0 << "    " << 0.0 << "\n";
//     cout << setprecision(12) << dFdU_21 << "    " << dFdU_22 << "    " << dFdU_23 << "\n";
//     cout << setprecision(12) << dFdU_31 << "    " << dFdU_32 << "    " << dFdU_33 << "\n";
//     cout << endl;
    
    dFdU[0][0] = 0.0;
    dFdU[0][1] = 1.0;
    dFdU[0][2] = 0.0;
    dFdU[1][0] = dFdU_21;
    dFdU[1][1] = dFdU_22;
    dFdU[1][2] = dFdU_23;
    dFdU[2][0] = dFdU_31;
    dFdU[2][1] = dFdU_32;
    dFdU[2][2] = dFdU_33;
    
    if (fabs(1.0 - N1 /**cos(phi)*/) < TOLER_HYPER) {
        cout << "Flux Jacobian matrix not diagonalizable in the free-streaming limit !!!!" << endl;
    } else {
        // Compute left and right eigenvectors
        for (int i = 1; i <= 3; i++) {
            Setup_lc(lc_vec, i, r_E, N1, phi, Chi2, E_dChi2_dE, dChi2_df);
            Setup_rc(rc_vec, i, r_E, N1, phi, Chi2, E_dChi2_dE, dChi2_df);
        }
        
        // Checking to make sure Rc Lc is the identity matrix
//         long double temp_val_check_eig;
//         for (int i = 0; i < 3; i++) {
//             for (int j = 0; j < 3; j++) {
//                 temp_val_check_eig = 0.0;
//                 for (int k = 0; k < 3; k++) {
//                     temp_val_check_eig += lc_vec[i][k] * rc_vec[k][j];
//                 }
//                 if (i == j) {
//                     if (fabs(1.0 - temp_val_check_eig) > 1.0e-6) {
//                         cout << "Issue with  eigendecomposition" << endl;
//                         cout << "i = " << i << "  " << "j = " << j << "  " << "temp_val_check_eig = " << temp_val_check_eig << endl;
//                         
//                         cout << "lambda_1 = " << lambda_1 << "  " << "lambda_2 = " << lambda_2 << "  " << "lambda_3 = " << lambda_3 << "  " << endl;
//         
//                         cout << endl;
//                         cout << "rc_vec = " << endl;
//                         cout << setprecision(12) << rc_vec[0][0] << "    " << rc_vec[0][1] << "    " << rc_vec[0][2] << "\n";
//                         cout << setprecision(12) << rc_vec[1][0] << "    " << rc_vec[1][1] << "    " << rc_vec[1][2] << "\n";
//                         cout << setprecision(12) << rc_vec[2][0] << "    " << rc_vec[2][1] << "    " << rc_vec[2][2] << "\n";
//                         cout << endl;
//                         
//                         cout << endl;
//                         cout << "lc_vec = " << endl;
//                         cout << setprecision(12) << lc_vec[0][0] << "    " << lc_vec[0][1] << "    " << lc_vec[0][2] << "\n";
//                         cout << setprecision(12) << lc_vec[1][0] << "    " << lc_vec[1][1] << "    " << lc_vec[1][2] << "\n";
//                         cout << setprecision(12) << lc_vec[2][0] << "    " << lc_vec[2][1] << "    " << lc_vec[2][2] << "\n";
//                         cout << endl;
//                         
//                         exit(0);
//                     }
//                 } else {
//                     if (fabs(temp_val_check_eig) > 1.0e-6) {
//                         cout << "Issue with  eigendecomposition" << endl;
//                         cout << "i = " << i << "  " << "j = " << j << "  " << "temp_val_check_eig = " << temp_val_check_eig << endl;
//                         
//                         cout << "lambda_1 = " << lambda_1 << "  " << "lambda_2 = " << lambda_2 << "  " << "lambda_3 = " << lambda_3 << "  " << endl;
//                         
//                         cout << endl;
//                         cout << "rc_vec = " << endl;
//                         cout << setprecision(12) << rc_vec[0][0] << "    " << rc_vec[0][1] << "    " << rc_vec[0][2] << "\n";
//                         cout << setprecision(12) << rc_vec[1][0] << "    " << rc_vec[1][1] << "    " << rc_vec[1][2] << "\n";
//                         cout << setprecision(12) << rc_vec[2][0] << "    " << rc_vec[2][1] << "    " << rc_vec[2][2] << "\n";
//                         cout << endl;
//                         
//                         cout << endl;
//                         cout << "lc_vec = " << endl;
//                         cout << setprecision(12) << lc_vec[0][0] << "    " << lc_vec[0][1] << "    " << lc_vec[0][2] << "\n";
//                         cout << setprecision(12) << lc_vec[1][0] << "    " << lc_vec[1][1] << "    " << lc_vec[1][2] << "\n";
//                         cout << setprecision(12) << lc_vec[2][0] << "    " << lc_vec[2][1] << "    " << lc_vec[2][2] << "\n";
//                         cout << endl;
//                         
//                         exit(0);
//                     }
//                 }
//             }
        }
        
//         for (int i = 0; i < 3; i++) {
//             for (int j = 0; j < 3; j++) {
//                 if (i == j) {
//                     if (i == 0) {
//                         lambdas[i][j] = lambda_1;
//                     } else if (i == 1) {
//                         lambdas[i][j] = lambda_2;
//                     } else {
//                         lambdas[i][j] = lambda_3;
//                     }
//                 } else {
//                     lambdas[i][j] = 0.0;
//                 }
//             }
//         }
//         
//         // Check to make sure that dFdU = [rc] [lambdas] [lc]
//         for (int i = 0; i < 3; i++) {
//             for (int j = 0; j < 3; j++) {
//                 temp_val = 0.0;
//                 for (int k = 0; k < 3; k++) {
//                     for (int l = 0; l < 3; l++) {
//                         // Compute flux Jacobian based on eigen-decomposition
//                         temp_val += rc_vec[i][k] * lambdas[k][l] * lc_vec[l][j];
//                         //                     cout << "rc = " << rc_vec[k][l] << "  " <<  "lc = " << lc_vec[k][l] << "  " <<  "lambdas = " << lambdas[k][l] << "  " <<  "k = " << k << "  " <<  "l = " << l << endl;
//                     }
//                 }
//                 
//                 if (fabs(dFdU[i][j] - temp_val) > TOLER_HYPER) {
//                     cout << "Eigenstructure decomposition not correct: temp_val = " << temp_val << "  " <<  "dFdU[i][j] = " << dFdU[i][j] << "  " <<  "i = " << i << "  " <<  "j = " << j << endl;
//                     cout <<  "N1_1 = " << N1*cos(phi) << "  " <<  "N1_2 = " << N1*sin(phi) << endl;
//                     cout << "Eigenvalues are: lambda_1 = " << lambda_1 << "  " <<  "lambda_2 = " << lambda_2 << "  " <<  "lambda_3 = " << lambda_3 << endl;
//                     
//                     cout << endl;
//                     cout << 0.0 << "    " << 1.0 << "    " << 0.0 << "\n";
//                     cout << setprecision(12) << dFdU_21 << "    " << dFdU_22 << "    " << dFdU_23 << "\n";
//                     cout << setprecision(12) << dFdU_31 << "    " << dFdU_32 << "    " << dFdU_33 << "\n";
//                     cout << endl;
//                     
//                     for (int ii = 0; ii < 3; ii++) {
//                         for (int jj = 0; jj < 3; jj++) {
//                             cout << "rc = " << rc_vec[ii][jj] << "  " <<  "lc = " << lc_vec[ii][jj] << "  " <<  "lambdas = " << lambdas[ii][jj] << "  " <<  "ii = " << ii << "  " <<  "jj = " << jj << endl;
//                         }
//                     }
//                     exit(0);
//                 }
//             }
//         }
//     }
    
//             if (i == j) {
//                 if (fabs(1.0 - temp_val) > TOLER_HYPER) {
//                     cout << "Matrix product lc rc is not the identity: diagonal elements temp_val = " << temp_val << endl;
//                     exit(0);
//                 }
//             } else {
//                 if (fabs(temp_val) > TOLER_HYPER4) {
//                     cout << "Matrix product lc rc is not the identity: off-diagonal elements temp_val = " << temp_val << endl;
//                     exit(0);
//                 }
//             }
}
                
void Evaluate_eigenvalues_2D(const long double &r, const long double &N1, const long double &phi, const long double &Chi2, const long double &E_dChi2_dE, const long double &dChi2_dN1, long double &lambda_1, long double &lambda_2, long double &lambda_3) {
//     long double E, dr_dE;
//     long double E_dChi2_dE, dChi2_dN1;
//     long double fx, fy;
//     long double *dFdU;
//     bool error_flag;
//     int NVARS = 3;
//     long double Chi2;
//     
//     dFdU = new long double[NVARS*NVARS];
//     
//     cout << "r = " << r << "    " << "N1 = " << N1 << "    " << "phi = " << phi << "\n";
//     cout << "Chi2 = " << Chi2 << "\n";
//     cout << "E_dChi2_dE = " << E_dChi2_dE << "    " << "dChi2_dN1 = " << dChi2_dN1 << "\n";
//     
//     if (N1 > 1e-7) {
//         dFdU[0*NVARS+0] = 0.0;
//         dFdU[0*NVARS+1] = c;
//         dFdU[0*NVARS+2] = 0.0;
//         dFdU[1*NVARS+0] = c*dI2ij_dI0(1, 1, Chi2, E_dChi2_dE, r, N1, phi);
//         dFdU[1*NVARS+1] = c*dI2ij_dIi(1, 1, 1, Chi2, dChi2_dN1, r, N1, phi);
//         dFdU[1*NVARS+2] = c*dI2ij_dIi(1, 1, 2, Chi2, dChi2_dN1, r, N1, phi);
//         dFdU[2*NVARS+0] = c*dI2ij_dI0(1, 2, Chi2, E_dChi2_dE, r, N1, phi);
//         dFdU[2*NVARS+1] = c*dI2ij_dIi(1, 2, 1, Chi2, dChi2_dN1, r, N1, phi);
//         dFdU[2*NVARS+2] = c*dI2ij_dIi(1, 2, 2, Chi2, dChi2_dN1, r, N1, phi);
//     } else {
//         dFdU[0*NVARS+0] = 0.0;
//         dFdU[0*NVARS+1] = c;
//         dFdU[0*NVARS+2] = 0.0;
//         dFdU[1*NVARS+0] = c/3.0;
//         dFdU[1*NVARS+1] = 0.0;
//         dFdU[1*NVARS+2] = 0.0;
//         dFdU[2*NVARS+0] = 0.0;
//         dFdU[2*NVARS+1] = 0.0;
//         dFdU[2*NVARS+2] = 0.0;
//     }
//     
//     real_2d_array a;
//     real_1d_array wr, wi;
//     real_2d_array vl, vr;
//     a.setcontent(NVARS, NVARS, dFdU);
//     
//     error_flag = rmatrixevd(a, NVARS, 0, wr, wi, vl, vr);
//     
//     lambda_1 = wr[0];
//     lambda_2 = wr[1];
//     lambda_3 = wr[2];
//     
//     long double dChi2_dN1_temp;
//     if (fabs(wi[0]) != 0.0 || fabs(wi[1]) != 0.0 || fabs(wi[2]) != 0.0) {
//         cout << "r = " << r << "    " << "N1 = " << N1 << "    " << "phi = " << phi << "\n";
//         cout << "Chi2 = " << Chi2 << "\n";
//         cout << "E_dChi2_dE = " << E_dChi2_dE << "    " << "dChi2_dN1 = " << dChi2_dN1 << "    " << "dChi2_dN1_temp = " << dChi2_dN1_temp << "\n";
//         cout << "lambda_1 = " << lambda_1 << "    " << "lambda_2 = " << lambda_2 << "    " << "lambda_3 = " << lambda_3 << "\n";
//     
//         cout << "wi[0] = " << wi[0] << "    " << "wi[1] = " << wi[1] << "    " << "wi[2] = " << wi[2] << "\n";
//         cout << "\n";
//     }
//     
//     delete[] dFdU;
}

void Test_Finite_Difference() {
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
    
    long double dChi2_dN1_actual, d2Chi2_dN1_2_actual, d3Chi2_dN1_3_actual;
    long double f_knot = 1.0;
    long double f_val_temp;

    dChi2_dN1_actual = dChi2_dN1_Gray(f_knot);
    d2Chi2_dN1_2_actual = d2Chi2_dN1_2_Gray(f_knot);
    d3Chi2_dN1_3_actual = d3Chi2_dN1_3_Gray(f_knot);
    
    // backward finite difference approximation to third derivative
    differ_backward ( h, order, prec, c, x );
    
    for (int id_diff = 0; id_diff < n; id_diff++) {
        f_val_temp = f_knot + x[id_diff];
        
        Chi2_vals[id_diff] = Chi_2_Gray(f_val_temp);
    }
    
    diff_3_chi2_dN1_3 = 0.0;
    for (int i = 0; i < n; i++) {
        diff_3_chi2_dN1_3 += c[n-i-1] * Chi2_vals[nmax-i-1];
    }
    
    // backward finite difference approximation to second derivative
    order = 2;
    n = order + prec;
    differ_backward ( h, order, prec, c, x );
    diff_2_chi2_dN1_2 = 0.0;
    for (int i = 0; i < n; i++) {
        diff_2_chi2_dN1_2 += c[n-i-1] * Chi2_vals[nmax-i-1];
    }
    
    // backward finite difference approximation to first derivative
    order = 1;
    n = order + prec;
    differ_backward ( h, order, prec, c, x );
    diff_chi2_dN1 = 0.0;
    for (int i = 0; i < n; i++) {
        diff_chi2_dN1 += c[n-i-1] * Chi2_vals[nmax-i-1];
    }
    
    cout << "N1 = " << f_knot << "   " << "dChi2_dN1 = " << dChi2_dN1_actual << "   " << "d2Chi2_dN1_2 = " << d2Chi2_dN1_2_actual << "   " << "d3Chi2_dN1_3 = " << d3Chi2_dN1_3_actual << endl;
    
    cout << "N1 = " << f_knot << "   " << "diff dChi2_dN1 = " << dChi2_dN1_actual - diff_chi2_dN1 << "   " << "diff d2Chi2_dN1_2 = " << d2Chi2_dN1_2_actual - diff_2_chi2_dN1_2 << "   " << "diff d3Chi2_dN1_3 = " << d3Chi2_dN1_3_actual - diff_3_chi2_dN1_3 << endl;
    
    delete[] c;
    delete[] x;
    delete[] Chi2_vals;
    
    exit(0);
 }
