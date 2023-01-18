#include "../../Non_Gray_Model/Chebyshev_Interpolation/M1_Model_1D.h"

double Lebedev_Quadrature_Orthog_SH(const int &index_SH_1, const int &index_SH_2, const int &Quad_Rule);
double Lebedev_Quadrature_Matrix_SH_Temp(const double *Matrix, const double *x_Fit, const double *y_Fit, const double *z_Fit, const int &index_SH, const int &Quad_Rule, const int &index, const int &N_points);
void Orthogonality_Test(const int &Quad_Rule);

void LU_Solve_Back(const int &n, double *x, const double *A, const double *b);
void LU_Solve_Forward(const int &n, double *x, const double *A, const double *b);
// void Cholesky_Factorization_Hessian(double *A, const int &n);
void LUdecomposition(const double *A, double *L, double *U, const int &n);
void Solve_A_x_b(const double *VanderMonde_Matrix, double *Coeff_Vand_Matrix, const double *VanderMonde_Vector, const int &n, const int &index_Cheby);
double Cartesian_Spherical_harmonics(const double &x, const double &y, const double &z, const int &Index);
double Test_Spherical_Harmonics(const double *Vals, const int &index_Cheby, const int &N_points_Cheby, const double &x, const double &y, const double &z, const int &N_Coeffs_SH);
void Duffy_Transformation(double &zeta_p, double &eta_p, const double &zeta, const double &eta);
double VanderMonde_Matrix_2_vars(const double &var_1, const double &var_2, const int &Order_poly_1, const int &Order_poly_2);
double VanderMonde_Vector_2_vars(const int &Index_Entry, const int &Index_Point);
double Lebedev_Quadrature_Matrix_SH(const double *Matrix, const int &index_SH, const int &Quad_Rule);

int main() {
    char path1_out[256], prefix[256], extension[256];
    double Chi_2_Fit;
    double E_plus_Fit, N1_1_plus_Fit, N2_11_plus_Fit, N2_12_plus_Fit;
    int N_Coeffs_E, N_Coeffs_f, N_Coeffs_Lebed;
    int index, index_rec_Chi_2, index_Ncoeffs, index_Cheby, index_full;
    int index_Cheby_lin, index_Cheby_col;
    double ratio_E;
    record_Chi2 rec_Chi_2;
    record_Ncoeffs rec_Ncoeffs;
    double *VanderMonde_Matrix, *VanderMonde_Vector, *Coeff_Vand_Matrix;
    double norm_f;
    double L2_Norm_Chi_2, L_inf_Norm_Chi_2, max_err_Chi_2_E, max_err_Chi_2_N1_1;
    
    double L2_Norm_E_plus, L2_Norm_N1_1_plus, L2_Norm_N2_11_plus, L2_Norm_N2_12_plus;
    double L_inf_Norm_E_plus, L_inf_Norm_N1_1_plus, L_inf_Norm_N2_11_plus, L_inf_Norm_N2_12_plus;
    
    double *E_NON_GRAY, *N1_1_NON_GRAY, *Chi_2_NON_GRAY;
    double *E_plus_NON_GRAY, *N1_1_plus_NON_GRAY, *N2_11_plus_NON_GRAY, *N2_12_plus_NON_GRAY;
    double *E_NON_GRAY_Cheby, *N1_1_NON_GRAY_Cheby, *Chi_2_NON_GRAY_Cheby;
    double *E_plus_NON_GRAY_Cheby, *N1_1_plus_NON_GRAY_Cheby, *N2_11_plus_NON_GRAY_Cheby, *N2_12_plus_NON_GRAY_Cheby;
    double *Coeff_Matrix_Chi_2, *Coefficient_Matrix_Fit_Chi_2;
    double *Coefficient_Matrix_Fit_E_plus, *Coefficient_Matrix_Fit_N1_1_plus, *Coefficient_Matrix_Fit_N2_11_plus, *Coefficient_Matrix_Fit_N2_12_plus;
    int id_max = 15;
    double Y_l_m;
    double max_err_E = 0, max_err_N1 = 0;
    int N_points_E, N_points_f, N_points_Theta, N_points_Phi;
    int N_points_Non_Gray_Uniform, N_points_Non_Gray_Cheby;
    
    N_points_E = 1;
    N_points_f = 100;
    N_points_Theta = 20;
    N_points_Phi = 20;
    int N_Coeffs_SH = 25;
    
    // Orthogonality_Test(8);
    
    N_points_Non_Gray_Uniform = N_points_E*N_points_f;
    
    E_NON_GRAY = new double [N_points_Non_Gray_Uniform];
    N1_1_NON_GRAY = new double [N_points_Non_Gray_Uniform];
    Chi_2_NON_GRAY = new double [N_points_Non_Gray_Uniform];
    E_plus_NON_GRAY = new double [N_points_Non_Gray_Uniform];
    N1_1_plus_NON_GRAY = new double [N_points_Non_Gray_Uniform];
    N2_11_plus_NON_GRAY = new double [N_points_Non_Gray_Uniform];
    N2_12_plus_NON_GRAY = new double [N_points_Non_Gray_Uniform];
    
    
    char path_Chi_2_out[256], path_E_plus_out[256], path_N1_1_plus_out[256], path_N2_11_plus_out[256], path_N2_12_plus_out[256];
    ofstream out_Chi2, out_E_plus, out_N1_1_plus, out_N2_11_plus, out_N2_12_plus;
    
    fstream in1_out;
    int Lebedev_Order;
    
    strcpy(path1_out, getenv(PATHVAR));
    strcat(path1_out, "/M1_Model/Gray_Model/Chebyshev_Interpolation/Uniform_Nodes_M1");
    sprintf(extension, "_%.6d", 0);
    strcat(extension, ".dat");
    strcat(path1_out, extension);
    
    in1_out.open(path1_out, ios::in|ios::binary);
        
    if (!in1_out) {
        cout << "Uniform_Nodes_M1.dat could not be accessed!" << endl;
    }
    
    if (in1_out.good()) {
        for (int index_e = 0 ; index_e < 1; index_e++){
            for (int i_f = 0; i_f < N_points_f; i_f++) {
                
                index_rec_Chi_2 = index_e*N_points_f + i_f;
                read<record_Chi2>(in1_out, rec_Chi_2, index_rec_Chi_2);
                
                E_NON_GRAY[index_rec_Chi_2] = rec_Chi_2.ratio_I0;
                N1_1_NON_GRAY[index_rec_Chi_2] = rec_Chi_2.N1;
                Chi_2_NON_GRAY[index_rec_Chi_2] = rec_Chi_2.Chi2;
                
                E_plus_NON_GRAY[index] = rec_Chi_2.I0_plus;
                N1_1_plus_NON_GRAY[index] = rec_Chi_2.N1_plus;
                N2_11_plus_NON_GRAY[index] = rec_Chi_2.N2_11_plus;
                N2_12_plus_NON_GRAY[index] = rec_Chi_2.N2_12_plus;
            }
        }
        in1_out.close();
    }
    
    
    strcpy(path1_out, getenv(PATHVAR));
    strcat(path1_out, "/M1_Model/Gray_Model/Chebyshev_Interpolation/Cheby_Nodes_M1");
    sprintf(extension, "_%.6d", 0);
    strcat(extension, ".dat");
    strcat(path1_out, extension);
    
    fstream in2_out;
    in2_out.open(path1_out, ios::in|ios::binary);
    
    if (!in2_out) {
        cout << "Cheby_Nodes_M1.dat could not be accessed!" << endl;
    }
    
    index_rec_Chi_2 = 0;
    index_Ncoeffs = 0;
    
    if (in2_out.good()) {
        for (int id_E = 1; id_E <= 1; id_E++) {
            for (int id_f = 1; id_f <= id_max; id_f++) {
                
                read<record_Ncoeffs, record_Chi2>(in2_out, rec_Ncoeffs, index_Ncoeffs, index_rec_Chi_2);
                index_Ncoeffs = index_Ncoeffs + 1;
                N_Coeffs_E = rec_Ncoeffs.N_Coeffs_E;
                N_Coeffs_f = rec_Ncoeffs.N_Coeffs_f;
                
                cout << "Set of Coefficients .............." << "N_Coeffs_E = " << N_Coeffs_E << "   " << "N_Coeffs_f = " << N_Coeffs_f << endl;
                
                N_points_Non_Gray_Cheby = N_Coeffs_E*N_Coeffs_f;
                            
                E_NON_GRAY_Cheby = new double [N_points_Non_Gray_Cheby];
                N1_1_NON_GRAY_Cheby = new double [N_points_Non_Gray_Cheby];
                Chi_2_NON_GRAY_Cheby = new double [N_points_Non_Gray_Cheby];
                
                E_plus_NON_GRAY_Cheby = new double [N_points_Non_Gray_Cheby];
                N1_1_plus_NON_GRAY_Cheby = new double [N_points_Non_Gray_Cheby];
                N2_11_plus_NON_GRAY_Cheby = new double [N_points_Non_Gray_Cheby];
                N2_12_plus_NON_GRAY_Cheby = new double [N_points_Non_Gray_Cheby];
                
                VanderMonde_Matrix = new double[N_points_Non_Gray_Cheby*N_points_Non_Gray_Cheby];
                
                Coeff_Vand_Matrix = new double[N_points_Non_Gray_Cheby*N_points_Non_Gray_Cheby];
                
                VanderMonde_Vector = new double[N_points_Non_Gray_Cheby];

                for (int i_Cheby_E = 0; i_Cheby_E < N_Coeffs_E; i_Cheby_E++) {
                    for (int i_Cheby_f = 0; i_Cheby_f < N_Coeffs_f; i_Cheby_f++) {
                        
                        index = i_Cheby_E*N_Coeffs_f + i_Cheby_f;
                        
                        read<record_Chi2, record_Ncoeffs>(in2_out, rec_Chi_2, index_rec_Chi_2, index_Ncoeffs);
                        
                        E_NON_GRAY_Cheby[index] = rec_Chi_2.ratio_I0;
                        N1_1_NON_GRAY_Cheby[index] = rec_Chi_2.N1;
                        Chi_2_NON_GRAY_Cheby[index] = rec_Chi_2.Chi2;
                        
                        Chi_2_NON_GRAY_Cheby[index] = log((3.0*Chi_2_NON_GRAY_Cheby[index] - 1.0)/2.0)/log(N1_1_NON_GRAY_Cheby[index]);
                        
                        E_plus_NON_GRAY_Cheby[index] = rec_Chi_2.I0_plus;
                        N1_1_plus_NON_GRAY_Cheby[index] = rec_Chi_2.N1_plus;
                        N2_11_plus_NON_GRAY_Cheby[index] = rec_Chi_2.N2_11_plus;
                        N2_12_plus_NON_GRAY_Cheby[index] = rec_Chi_2.N2_12_plus;
                        index_rec_Chi_2 = index_rec_Chi_2 + 1;
                    } // end for i_Cheby_f
                } // end for i_Cheby_E
                                
                Coefficient_Matrix_Fit_Chi_2 = new double [N_points_Non_Gray_Cheby];
                Coefficient_Matrix_Fit_E_plus = new double [N_points_Non_Gray_Cheby];
                Coefficient_Matrix_Fit_N1_1_plus = new double [N_points_Non_Gray_Cheby];
                Coefficient_Matrix_Fit_N2_11_plus = new double [N_points_Non_Gray_Cheby];
                Coefficient_Matrix_Fit_N2_12_plus = new double [N_points_Non_Gray_Cheby];
                
                // Chebyshev quadrature
                cout  << "*************************** Setting up the Vandermonde System ************************"<< endl;
                cout  << "*************************** Performing Chebyshev Quadrature ************************"<< endl;
                // Chebyshev quadrature
                for (int i_lin_Cheby_E = 0; i_lin_Cheby_E < N_Coeffs_E; i_lin_Cheby_E++) {
                    for (int i_lin_Cheby_f = 0; i_lin_Cheby_f < N_Coeffs_f; i_lin_Cheby_f++) {
                        index_Cheby_lin = i_lin_Cheby_E*N_Coeffs_f + i_lin_Cheby_f;
                        ratio_E = E_NON_GRAY_Cheby[index_Cheby_lin];
                        norm_f = N1_1_NON_GRAY_Cheby[index_Cheby_lin];
                        
                        for (int i_col_Cheby_E = 0; i_col_Cheby_E < N_Coeffs_E; i_col_Cheby_E++) {
                            for (int i_col_Cheby_f = 0; i_col_Cheby_f < N_Coeffs_f; i_col_Cheby_f++) {
                                
                                index_Cheby_col = i_col_Cheby_E*N_Coeffs_f + i_col_Cheby_f;
                                
                                VanderMonde_Matrix[index_Cheby_lin*N_points_Non_Gray_Cheby+index_Cheby_col] = VanderMonde_Matrix_2_vars(2.0*ratio_E-1.0, 2.0*norm_f-1.0, i_col_Cheby_E, i_col_Cheby_f);
                                
                            } // end for i_col_Cheby_f
                        } // end for i_col_Cheby_E
                                     
                    } // end for i_lin_Cheby_f
                } // end for i_lin_Cheby_E
                
                for (int i_Cheby_E = 0; i_Cheby_E < N_Coeffs_E; i_Cheby_E++) {
                    for (int i_Cheby_f = 0; i_Cheby_f < N_Coeffs_f; i_Cheby_f++) {
                        index_Cheby = i_Cheby_E*N_Coeffs_f + i_Cheby_f;
                        
                        for (int i_Coeff_Cheby = 0; i_Coeff_Cheby < N_points_Non_Gray_Cheby; i_Coeff_Cheby++) {
                            VanderMonde_Vector[i_Coeff_Cheby] = VanderMonde_Vector_2_vars(i_Coeff_Cheby, index_Cheby);
                        } // end for i_Coeff_Cheby
                        
                        // cout  << "*************************** solving the Vandermonde System ************************"<< endl;
                        
                        Solve_A_x_b(VanderMonde_Matrix, Coeff_Vand_Matrix, VanderMonde_Vector, N_points_Non_Gray_Cheby, index_Cheby);
                    } // end for i_Cheby_f
                } // end for i_Cheby_E
                
                for (int i_Coeff_Cheby = 0; i_Coeff_Cheby < N_points_Non_Gray_Cheby; i_Coeff_Cheby++) {
                    Coefficient_Matrix_Fit_Chi_2[i_Coeff_Cheby] = 0.0;
                    Coefficient_Matrix_Fit_E_plus[i_Coeff_Cheby] = 0.0;
                    Coefficient_Matrix_Fit_N1_1_plus[i_Coeff_Cheby] = 0.0;
                    Coefficient_Matrix_Fit_N2_11_plus[i_Coeff_Cheby] = 0.0;
                    Coefficient_Matrix_Fit_N2_12_plus[i_Coeff_Cheby] = 0.0;
                    
                    for (int j_Coeff_Cheby = 0; j_Coeff_Cheby < N_points_Non_Gray_Cheby; j_Coeff_Cheby++) {
                        Coefficient_Matrix_Fit_Chi_2[i_Coeff_Cheby] += Coeff_Vand_Matrix[i_Coeff_Cheby*N_points_Non_Gray_Cheby+j_Coeff_Cheby]*Chi_2_NON_GRAY_Cheby[j_Coeff_Cheby];
                    
                        Coefficient_Matrix_Fit_E_plus[i_Coeff_Cheby] += Coeff_Vand_Matrix[i_Coeff_Cheby*N_points_Non_Gray_Cheby+j_Coeff_Cheby]*E_plus_NON_GRAY_Cheby[j_Coeff_Cheby];
                        
                        Coefficient_Matrix_Fit_N1_1_plus[i_Coeff_Cheby] += Coeff_Vand_Matrix[i_Coeff_Cheby*N_points_Non_Gray_Cheby+j_Coeff_Cheby]*N1_1_plus_NON_GRAY_Cheby[j_Coeff_Cheby];
                        
                        Coefficient_Matrix_Fit_N2_11_plus[i_Coeff_Cheby] += Coeff_Vand_Matrix[i_Coeff_Cheby*N_points_Non_Gray_Cheby+j_Coeff_Cheby]*N2_11_plus_NON_GRAY_Cheby[j_Coeff_Cheby];
                        
                        Coefficient_Matrix_Fit_N2_12_plus[i_Coeff_Cheby] += Coeff_Vand_Matrix[i_Coeff_Cheby*N_points_Non_Gray_Cheby+j_Coeff_Cheby]*N2_12_plus_NON_GRAY_Cheby[j_Coeff_Cheby];
                    }
                }
                
                cout  << "*************************** Chebyshev Quadrature Completed ************************"<< endl;
                delete[] Coeff_Vand_Matrix;
                delete[] VanderMonde_Matrix;
                delete[] VanderMonde_Vector;
                            
                cout  << "*************************** Testing fit based on selected node for interpolation ************************"<< endl;
                L2_Norm_Chi_2 = 0.0;
                L_inf_Norm_Chi_2 = 0.0;
                
                for (int i_Cheby_E = 0; i_Cheby_E < N_Coeffs_E; i_Cheby_E++) {
                    for (int i_Cheby_f = 0; i_Cheby_f < N_Coeffs_f; i_Cheby_f++) {
                        index = i_Cheby_E*N_Coeffs_f + i_Cheby_f;
                        
                        Chi_2_Fit = 0.0;
                        for (int i_fit_E = 0; i_fit_E < N_Coeffs_E; i_fit_E++) {
                            for (int i_fit_f = 0; i_fit_f < N_Coeffs_f; i_fit_f++) {
                                index_Cheby = i_fit_E*N_Coeffs_f + i_fit_f;
                                
                                ratio_E = E_NON_GRAY_Cheby[index];
                                norm_f = N1_1_NON_GRAY_Cheby[index];
                                
                                Chi_2_Fit += Coefficient_Matrix_Fit_Chi_2[index_Cheby]*VanderMonde_Matrix_2_vars(2.0*ratio_E-1.0, 2.0*norm_f-1.0, i_fit_E, i_fit_f);
                                
                            }
                        }
                        
                        L2_Norm_Chi_2 += pow(Chi_2_Fit - Chi_2_NON_GRAY_Cheby[index], 2);
                        L_inf_Norm_Chi_2 = max(L_inf_Norm_Chi_2, fabs(Chi_2_Fit - Chi_2_NON_GRAY_Cheby[index]));
                                                        
                        if (L2_Norm_Chi_2 == fabs(Chi_2_Fit - Chi_2_NON_GRAY[index])) {
                            max_err_Chi_2_E = E_NON_GRAY[index];
                            max_err_Chi_2_N1_1 = N1_1_NON_GRAY[index];
                        }
                        
                    }
                }
                
                delete[] E_NON_GRAY_Cheby;
                delete[] N1_1_NON_GRAY_Cheby;
                delete[] Chi_2_NON_GRAY_Cheby;
                            
                delete[] E_plus_NON_GRAY_Cheby;
                delete[] N1_1_plus_NON_GRAY_Cheby;
                delete[] N2_11_plus_NON_GRAY_Cheby;
                delete[] N2_12_plus_NON_GRAY_Cheby;
                            
                cout  << "*************************** Testing fit based on selected node for interpolation Completed ************************"<< endl;
                
                cout << "Errors ................." << "L2_Norm_Chi_2 = " << L2_Norm_Chi_2 << "     "  << "L_inf_Norm_Chi_2 = " << L_inf_Norm_Chi_2 << endl;
                            
                L2_Norm_Chi_2 = 0.0;
                L_inf_Norm_Chi_2 = 0.0;
                
//                 cout << "Set of Coefficients .............." << "N_Coeffs_E = " << N_Coeffs_E << "   " << "N_Coeffs_f = " << N_Coeffs_f << endl;
                
                for (int i = 0; i < N_points_E ; i++) {
                    for (int j = 0; j < N_points_f ; j++) {
                        index = i*N_points_f + j;
                        
                        Chi_2_Fit = 0.0;
                        E_plus_Fit = 0.0;
                        N1_1_plus_Fit = 0.0;
                        N2_11_plus_Fit = 0.0;
                        N2_12_plus_Fit = 0.0;
                        
                        ratio_E = E_NON_GRAY[index];
                        norm_f = N1_1_NON_GRAY[index];
                                
                        for (int i_fit_E = 0; i_fit_E < N_Coeffs_E; i_fit_E++) {
                            for (int i_fit_f = 0; i_fit_f < N_Coeffs_f; i_fit_f++) {
                                index_Cheby = i_fit_E*N_Coeffs_f + i_fit_f;
                                
                                Chi_2_Fit += Coefficient_Matrix_Fit_Chi_2[index_Cheby]*VanderMonde_Matrix_2_vars(2.0*ratio_E-1.0, 2.0*norm_f-1.0, i_fit_E, i_fit_f);
                                
//                                 E_plus_Fit += Coefficient_Matrix_Fit_E_plus[index_Cheby]*VanderMonde_Matrix_2_vars(2.0*ratio_E-1.0, 2.0*norm_f-1.0, i_fit_E, i_fit_f);
//                                 
//                                 N1_1_plus_Fit += Coefficient_Matrix_Fit_N1_1_plus[index_Cheby]*VanderMonde_Matrix_2_vars(2.0*ratio_E-1.0, 2.0*norm_f-1.0, i_fit_E, i_fit_f);
//                                                                          
//                                 N2_11_plus_Fit += Coefficient_Matrix_Fit_N2_11_plus[index_Cheby]*VanderMonde_Matrix_2_vars(2.0*ratio_E-1.0, 2.0*norm_f-1.0, i_fit_E, i_fit_f);
//                                 
//                                 N2_12_plus_Fit += Coefficient_Matrix_Fit_N2_12_plus[index_Cheby]*VanderMonde_Matrix_2_vars(2.0*ratio_E-1.0, 2.0*norm_f-1.0, i_fit_E, i_fit_f);
                            }
                        }   
                        
                        Chi_2_Fit = (1.0 + 2.0*pow(norm_f, Chi_2_Fit))/3.0;
                        
//                         cout << "Chi_2_Fit = " << Chi_2_Fit << "     " << "Chi_2_NON_GRAY[index] = " << Chi_2_NON_GRAY[index] << "   " << "Diff = " << fabs(Chi_2_Fit - Chi_2_NON_GRAY[index]) << endl;
                        
                        L2_Norm_Chi_2 += pow(Chi_2_Fit - Chi_2_NON_GRAY[index], 2);
                        L_inf_Norm_Chi_2 = max(L_inf_Norm_Chi_2, fabs(Chi_2_Fit - Chi_2_NON_GRAY[index]));
                                                     
                        L2_Norm_E_plus += pow(E_plus_Fit - E_plus_NON_GRAY[index], 2);
                        L_inf_Norm_E_plus = max(L_inf_Norm_E_plus, fabs(E_plus_Fit - E_plus_NON_GRAY[index]));
                                                     
                        L2_Norm_N1_1_plus += pow(N1_1_plus_Fit - N1_1_plus_NON_GRAY[index], 2);
                        L_inf_Norm_N1_1_plus = max(L_inf_Norm_N1_1_plus, fabs(N1_1_plus_Fit - N1_1_plus_NON_GRAY[index]));
                                                     
                        L2_Norm_N2_11_plus += pow(N2_11_plus_Fit - N2_11_plus_NON_GRAY[index], 2);
                        L_inf_Norm_N2_11_plus = max(L_inf_Norm_N2_11_plus, fabs(N2_11_plus_Fit - N2_11_plus_NON_GRAY[index]));
                                                     
                        L2_Norm_N2_12_plus += pow(N2_12_plus_Fit - N2_12_plus_NON_GRAY[index], 2);
                        L_inf_Norm_N2_12_plus = max(L_inf_Norm_N2_12_plus, fabs(N2_12_plus_Fit - N2_12_plus_NON_GRAY[index]));
                        
                        if (L_inf_Norm_Chi_2 == fabs(Chi_2_Fit - Chi_2_NON_GRAY[index])) {
                            max_err_Chi_2_E = E_NON_GRAY[index];
                            max_err_Chi_2_N1_1 = N1_1_NON_GRAY[index];
                        }
                    }
                }
                
                cout << "Errors ................." << "L2_Norm_Chi_2 = " << L2_Norm_Chi_2 << "     "  << "L_inf_Norm_Chi_2 = " << L_inf_Norm_Chi_2 << "    " << "max_err_Chi_2_E = " << max_err_Chi_2_E << "    " << "max_err_Chi_2_N1_1 = " << max_err_Chi_2_N1_1 << endl;
                
//                 cout << "Errors ................." << "L2_Norm_E_plus = " << L2_Norm_E_plus << "     "  << "L_inf_Norm_E_plus = " << L_inf_Norm_E_plus << endl;
//                 
//                 cout << "Errors ................." << "L2_Norm_N1_1_plus = " << L2_Norm_N1_1_plus << "     "  << "L_inf_Norm_N1_1_plus = " << L_inf_Norm_N1_1_plus << endl; 
//                 
//                 cout << "Errors ................." << "L2_Norm_N2_11_plus = " << L2_Norm_N2_11_plus << "     "  << "L_inf_Norm_N2_11_plus = " << L_inf_Norm_N2_11_plus << endl;
//                 
//                 cout << "Errors ................." << "L2_Norm_N2_12_plus = " << L2_Norm_N2_12_plus << "     "  << "L_inf_Norm_N2_12_plus = " << L_inf_Norm_N2_12_plus << endl;  
//                 
//                 //                              cout << "Location ................." << "max_err_N3_111_E = " << max_err_N3_111_E << "     "  << "max_err_N1 = " << max_err_N1 << endl;
                cout << endl;
                            
                strcpy(path_Chi_2_out, getenv(PATHVAR));
                strcat(path_Chi_2_out, "/M1_Model/Non_Gray_Model/Chebyshev_Interpolation/Coefficients_Fit_Chi_2_Non_Gray_BE.dat");
                
                strcpy(path_E_plus_out, getenv(PATHVAR));
                strcat(path_E_plus_out, "/M1_Model/Non_Gray_Model/Chebyshev_Interpolation/Coefficients_Fit_E_plus_Non_Gray_BE.dat");
                            
                strcpy(path_N1_1_plus_out, getenv(PATHVAR));
                strcat(path_N1_1_plus_out, "/M1_Model/Non_Gray_Model/Chebyshev_Interpolation/Coefficients_Fit_Fx_plus_Non_Gray_BE.dat");
                            
                strcpy(path_N2_11_plus_out, getenv(PATHVAR));
                strcat(path_N2_11_plus_out, "/M1_Model/Non_Gray_Model/Chebyshev_Interpolation/Coefficients_Fit_Pxx_plus_Non_Gray_BE.dat");
                            
                strcpy(path_N2_12_plus_out, getenv(PATHVAR));
                strcat(path_N2_12_plus_out, "/M1_Model/Non_Gray_Model/Chebyshev_Interpolation/Coefficients_Fit_Pxy_plus_Non_Gray_BE.dat");
                            
                out_Chi2.open(path_Chi_2_out);
                out_E_plus.open(path_E_plus_out);
                out_N1_1_plus.open(path_N1_1_plus_out);
                out_N2_11_plus.open(path_N2_11_plus_out);
                out_N2_12_plus.open(path_N2_12_plus_out);
                
                if (!out_Chi2 || !out_E_plus || !out_N1_1_plus || !out_N2_11_plus || !out_N2_12_plus) {
                    cout << "One of the output files could not be accessed" << endl;
                }
                
                for (int index_Cheby = 0; index_Cheby < N_points_Non_Gray_Cheby; index_Cheby++) {
                    out_Chi2 << Coefficient_Matrix_Fit_Chi_2[index_Cheby] << endl;
                    
                    out_E_plus << Coefficient_Matrix_Fit_E_plus[index_Cheby] << endl;
                    
                    out_N1_1_plus << Coefficient_Matrix_Fit_N1_1_plus[index_Cheby] << endl;
                    
                    out_N2_11_plus << Coefficient_Matrix_Fit_N2_11_plus[index_Cheby] << endl;
                    
                    out_N2_12_plus << Coefficient_Matrix_Fit_N2_12_plus[index_Cheby] << endl;
                }
                
                out_Chi2.close();
                out_E_plus.close();
                out_N1_1_plus.close();
                out_N2_11_plus.close();
                out_N2_12_plus.close();
                
                delete[] Coefficient_Matrix_Fit_Chi_2;
                delete[] Coefficient_Matrix_Fit_E_plus;
                delete[] Coefficient_Matrix_Fit_N1_1_plus;
                delete[] Coefficient_Matrix_Fit_N2_11_plus;
                delete[] Coefficient_Matrix_Fit_N2_12_plus;
                
            } // end for id_f
        } // end for id_E
    } // end if for files
    in2_out.close();
    
    delete[] E_NON_GRAY; E_NON_GRAY = NULL;
    delete[] N1_1_NON_GRAY; N1_1_NON_GRAY= NULL;
    delete[] Chi_2_NON_GRAY; Chi_2_NON_GRAY= NULL;
    delete[] E_plus_NON_GRAY; E_plus_NON_GRAY = NULL;
    delete[] N1_1_plus_NON_GRAY; N1_1_plus_NON_GRAY = NULL;
    delete[] N2_11_plus_NON_GRAY; N2_11_plus_NON_GRAY = NULL;
    delete[] N2_12_plus_NON_GRAY; N2_12_plus_NON_GRAY = NULL;
}

double Cartesian_Spherical_harmonics(const double &x, const double &y, const double &z, const int &Index) {
    double f_SH;
    double y_0_0;
    double y_1_m1, y_1_0, y_1_1;
    double y_2_m2, y_2_m1, y_2_0, y_2_1, y_2_2;
    double y_3_m3, y_3_m2, y_3_m1, y_3_0, y_3_1, y_3_2, y_3_3;
    double y_4_m4, y_4_m3, y_4_m2, y_4_m1, y_4_0, y_4_1, y_4_2, y_4_3, y_4_4;
    
    y_0_0 = sqrt(1.0/(4.0*PI));
    y_1_m1 = sqrt(3.0/(4.0*PI))*y;
    y_1_0 = sqrt(3.0/(4.0*PI))*z;
    y_1_1 = sqrt(3.0/(4.0*PI))*x;
    
    y_2_m2 = sqrt(15.0/(4.0*PI))*x*y;
    y_2_m1 = sqrt(15.0/(4.0*PI))*y*z;
    y_2_0 = sqrt(5.0/(16.0*PI))*(3.0*pow(z,2) - 1.0);
    y_2_1 = sqrt(15.0/(4.0*PI))*x*z;
    y_2_2 = sqrt(15.0/(16.0*PI))*(pow(x,2) - pow(y,2));
    
    y_3_m3 = sqrt(35.0/(32.0*PI))*y*(3.0*pow(x,2) - pow(y,2));
    y_3_m2 = sqrt(105.0/(4.0*PI))*x*y*z;
    y_3_m1 = sqrt(21.0/(32.0*PI))*y*(5.0*pow(z,2) - 1.0);
    y_3_0 = sqrt(7.0/(16.0*PI))*z*(5.0*pow(z,2) - 3.0);
    y_3_1 = sqrt(21.0/(32.0*PI))*x*(5.0*pow(z,2) - 1.0);
    y_3_2 = sqrt(105.0/(16.0*PI))*z*(pow(x,2) - pow(y,2));
    y_3_3 = sqrt(35.0/(32.0*PI))*x*(pow(x,2) - 3.0*pow(y,2));
    
    y_4_m4 = 3.0*sqrt(35.0/(16.0*PI))*x*y*(pow(x,2) - pow(y,2));
    y_4_m3 = 3.0*sqrt(35.0/(32.0*PI))*y*z*(3.0*pow(x,2) - pow(y,2));
    y_4_m2 = 3.0*sqrt(5.0/(16.0*PI))*x*y*(7.0*pow(z,2) - 1.0);
    y_4_m1 = 3.0*sqrt(5.0/(32.0*PI))*y*z*(7.0*pow(z,2) - 3.0);
    y_4_0 = (3.0/16.0)*sqrt(1.0/PI)*(35.0*pow(z,4) - 30.0*pow(z,2) + 3.0);
    y_4_1 = 3.0*sqrt(5.0/(32.0*PI))*x*z*(7.0*pow(z,2) - 3.0);
    y_4_2 = (3.0/8.0)*sqrt(5.0/PI)*(pow(x,2) - pow(y,2))*(7.0*pow(z,2) - 1.0);
    y_4_3 = (3.0/4.0)*sqrt(35.0/(2.0*PI))*x*z*(pow(x,2) - 3.0*pow(y,2));
    y_4_4 = (3.0/16.0)*sqrt(35.0/PI)*(pow(x,2)*(pow(x,2) - 3.0*pow(y,2)) - pow(y,2)*(3.0*pow(x,2) - pow(y,2)));
    
    if (Index == 0) { 
        f_SH = y_0_0;
    } else if (Index == 1) {
        f_SH = y_1_m1;
    } else if (Index == 2) {
        f_SH = y_1_0;
    } else if (Index == 3) {
        f_SH = y_1_1;
    } else if (Index == 4) {
        f_SH = y_2_m2;
    } else if (Index == 5) {
        f_SH = y_2_m1;
    } else if (Index == 6) {
        f_SH = y_2_0;
    } else if (Index == 7) {
        f_SH = y_2_1;
    } else if (Index == 8) {
        f_SH = y_2_2;
    } else if (Index == 9) {
        f_SH = y_3_m3;
    } else if (Index == 10) {
        f_SH = y_3_m2;
    } else if (Index == 11) {
        f_SH = y_3_m1;
    } else if (Index == 12) {
        f_SH = y_3_0;
    } else if (Index == 13) {
        f_SH = y_3_1;
    } else if (Index == 14) {
        f_SH = y_3_2;
    } else if (Index == 15) {
        f_SH = y_3_3;
    } else if (Index == 16) {
        f_SH = y_4_m4;
    } else if (Index == 17) {
        f_SH = y_4_m3;
    } else if (Index == 18) {
        f_SH = y_4_m2;
    } else if (Index == 19) {
        f_SH = y_4_m1;
    } else if (Index == 20) {
        f_SH = y_4_0;
    } else if (Index == 21) {
        f_SH = y_4_1;
    } else if (Index == 22) {
        f_SH = y_4_2;
    } else if (Index == 23) {
        f_SH = y_4_3;
    } else if (Index == 24) {
        f_SH = y_4_4;
    }
    
    return f_SH;
}

double VanderMonde_Matrix_2_vars(const double &var_1, const double &var_2, const int &Order_poly_1, const int &Order_poly_2) {
    double Van_Entry;
    double poly_1, poly_2;
    poly_1 = Chebyshev_Polynomial_Basis(var_1, Order_poly_1);
    poly_2 = Chebyshev_Polynomial_Basis(var_2, Order_poly_2);
    
    Van_Entry = poly_1*poly_2;
    return Van_Entry;
}

double VanderMonde_Vector_2_vars(const int &Index_Entry, const int &Index_Point) {
    double Van_Vector;
    
    if (Index_Entry == Index_Point) {
        Van_Vector = 1.0;
    } else {
        Van_Vector = 0.0;
    }
    return Van_Vector;
}

void Solve_A_x_b(const double *VanderMonde_Matrix, double *Coeff_Vand_Matrix, const double *VanderMonde_Vector, const int &n, const int &index_Cheby) {
    // First perform Cholesky decomposition of A
    double *Temp_Mat_L, *Temp_Mat_U, *Temp_Vec;
    Temp_Mat_L = new double[n*n];
    Temp_Mat_U = new double[n*n];
    Temp_Vec = new double[n];
    
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


void LUdecomposition(const double *A, double *L, double *U, const int &n) {
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

void LU_Solve_Back(const int &n, double *x, const double *A, const double *b)
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

void LU_Solve_Forward(const int &n, double *x, const double *A, const double *b)
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

// double Lebedev_Quadrature_Matrix_SH_Temp(const double *Matrix, const double *x_Fit, const double *y_Fit, const double *z_Fit, const int &index_SH, const int &Quad_Rule, const int &index, const int &N_points) {
//     int order;
//     double *x, *y, *z, *w;
//     double integral_approx, Y_l_m;
//     double norm_f;
//     order = order_table ( Quad_Rule );
//     
//     x = new double[order];
//     y = new double[order];
//     z = new double[order];
//     w = new double[order];
//     
//     ld_by_order ( order, x, y, z, w );
//     
//     integral_approx = 0.0;
//     
//     for ( int m = 0; m < order; m++ ) {
//         Y_l_m = Cartesian_Spherical_harmonics(x[m], y[m], z[m], index_SH);
//         integral_approx = integral_approx + w[m] * Matrix[m] * Y_l_m;
//     }
//     
//     integral_approx = 4.0*PI*integral_approx;
//     
//     delete[] x;
//     delete[] y;
//     delete[] z;
//     delete[] w;
//     
//     return integral_approx;
// }
// 
// double Lebedev_Quadrature_Matrix_SH(const double *Matrix, const int &index_SH, const int &Quad_Rule) {
//     int order;
//     double *x, *y, *z, *w;
//     double integral_approx, Y_l_m;
//     
//     order = order_table ( Quad_Rule );
//     
//     x = new double[order];
//     y = new double[order];
//     z = new double[order];
//     w = new double[order];
//     
//     ld_by_order ( order, x, y, z, w );
//     
//     integral_approx = 0.0;
//     
//     for ( int m = 0; m < order; m++ ) {
//         Y_l_m = Cartesian_Spherical_harmonics(x[m], y[m], z[m], index_SH);
//         integral_approx = integral_approx + w[m] * Matrix[m] * Y_l_m;
//     }
//     
//     integral_approx = 4.0*PI*integral_approx;
//     
//     delete[] x;
//     delete[] y;
//     delete[] z;
//     delete[] w;
//     
//     return integral_approx;
// }
// 
// double Lebedev_Quadrature_Orthog_SH(const int &index_SH_1, const int &index_SH_2, const int &Quad_Rule) {
//     int order;
//     double *x, *y, *z, *w;
//     double integral_approx, Y_l_m_1, Y_l_m_2;
//     
//     order = order_table ( Quad_Rule );
//     
//     x = new double[order];
//     y = new double[order];
//     z = new double[order];
//     w = new double[order];
//     
//     ld_by_order ( order, x, y, z, w );
//     
//     integral_approx = 0.0;
//     
//     for ( int m = 0; m < order; m++ ) {
//         Y_l_m_1 = Cartesian_Spherical_harmonics(x[m], y[m], z[m], index_SH_1);
//         Y_l_m_2 = Cartesian_Spherical_harmonics(x[m], y[m], z[m], index_SH_2);
//         integral_approx = integral_approx + w[m] * Y_l_m_1 * Y_l_m_2;
//     }
//     
//     integral_approx = 4.0*PI*integral_approx;
//     
//     delete[] x;
//     delete[] y;
//     delete[] z;
//     delete[] w;
//     
//     return integral_approx;
// }

double Test_Spherical_Harmonics(const double *Vals, const int &index_Cheby, const int &N_points_Cheby, const double &x, const double &y, const double &z, const int &N_Coeffs_SH) {
    double f_SH;
    f_SH = 0.0;
    
    for (int i_SH = 0; i_SH < N_Coeffs_SH; i_SH++) {
        f_SH += Vals[i_SH*N_points_Cheby + index_Cheby]*Cartesian_Spherical_harmonics(x,y,z,i_SH);
    }
    return f_SH;
}

// void Orthogonality_Test(const int &Quad_Rule) {
//     double delta_ij;
//     for (int i = 0; i < 25; i++) {
//         for (int j = 0; j < 25; j++) {
//             delta_ij = Lebedev_Quadrature_Orthog_SH(i, j, 8);
//             cout << "i = " << i << "    "  << "j = " << j << "   " << "delta_ij = " << delta_ij << endl;
//         }
//     }
// }
