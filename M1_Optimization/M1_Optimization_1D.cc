#ifndef _M1_Model_1D_H_INCLUDED
#include "../M1_Optimization/M1_Optimization.h"
#endif // _M1_Model_1D_H_INCLUDED

// g++ -I/usr/local/include M1_NLOPT_TEST.cc -L/usr/local/lib -lnlopt -lm -o tuutorial

long double I0_val_max = 5.0*1.0e6;
    
long double L_vals[11] = {1.0e-6, 1.0e-1, 0.5, 1.0, 5.0, 1.0e1, 5.0e1, 1.0e2, 5.0e2, 1.0e3, 1.0e4};
long double E_vals[11] = {1.0e-6, 1.0e-1, 0.5, 1.0, 5.0, 1.0e1, 5.0e1, 1.0e2, 5.0e2, 1.0e3, 1.0e4};

long double f_L_Chi2_vals[10] = {0.0, 0.9, 0.94, 0.96, 0.97, 0.98, 0.99, 0.995, 0.998, 0.9995};
long double f_L_In_Plus_vals[10] = {0.0, 0.7, 0.9, 0.94, 0.96, 0.97, 0.98, 0.99, 0.994, 0.998};

long double r_l[4] = {0, 1.0e-8, 1.0e-6, 1.0e-4}; //declare extern in here and declare whithout extern in .cc file to make it a global variable
int Lebed_Rule_Set[6] = {10, 20, 32, 44, 50, 65};
int N_quad_points_Circle_Set[6] = {20, 30, 64, 127, 129, 255};
int M1_State_Param::Regime = GRAY;

const int M1_State_Param::max_refinement_level;

int M1_State_Param::Max_Ent_Solution_Type = CLOSING_FLUX;

int DISPLAY_ID = 0;
int PRIMARY_ID = 0;

long double finite_diff_h = 2.0*1.0e-2;
long double tol_grad = 1.0e-6;

// ********************************************************************************
// This routine calls the optimization algorithm used for solving the maximum-
// entropy problem for any given set of angular moments up to first-order, for 
// gray or non-gray radiation
// ********************************************************************************
void NLOPT_Optimization_Algorithm_NON_GRAY_M1(record_Chi2 *rec_Chi2_local, record_Partial_Moments *rec_PM_array, const int &index_count, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, const int &id_proc, M1_State_Param &M1_State, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, long double *Sk_final) {
    long double x[M1_State.NVARS];
    long double tol_x[M1_State.NVARS];
    int max_iters_nlopt, max_index_reg;
    
    max_index_reg = 2;
    max_iters_nlopt = 500;
    
    long double minf, min_grad; /* the minimum objective value, upon return */
    nlopt_opt opt = NULL;
    
    if (id_proc == DISPLAY_ID) {
        if (M1_State.display) {
            cout << "i_E = " << Var_index->E << "     "  << "i_f = " << Var_index->f << endl;
        }
    }
    
    Set_Moments(&M1_State, Var_index, num_points, &Mobius_Scale_Params, struc_map_N1);
    
    if (M1_State.flag_Realizability) {
        if (M1_State.display) {
            cout << ".......................Non Realizable Moments......................" << endl;
            cout << "i_E = " << Var_index->E << "     "  << "i_f = " << Var_index->f << endl;
            cout  << "I0 = " << M1_State.I0 << "   " << "N1 = " << M1_State.N1 << endl; 
        }
        goto Continue;
    } else if (M1_State.Boundary_Point == FREE_STREAMING_LIMIT) {
        if (M1_State.display) {
            if (M1_State.id_proc == DISPLAY_ID) {
                cout << "....................... Free-streaming limit ......................" << endl;
            }
        }
        goto Continue;
    } else if (M1_State.flag_Taylor_Series_Expansion) {
        if (M1_State.display) {
            if (M1_State.id_proc == DISPLAY_ID) {
                cout << "....................... Taylor Series near free-streaming limit ......................" << endl;
            }
        }
        goto Continue;
    }
    
    for (int index_reg = 0; index_reg < max_index_reg; index_reg++) {
        if (index_reg > 0) {
            if (M1_State.display) {
                if (id_proc == DISPLAY_ID) {
                    printf(".................Moment Regularization with r_l = %Le..................\n", r_l[index_reg]);
                }
            }
        }
        
        for ( int id_refinement = 0; id_refinement <= M1_State.max_refinement_level; id_refinement++) {
            M1_State.order_quad = 20;
            // Allocate and setup quadrature scheme
            M1_State.Allocate_Quad_1D();
            // Setup the bounds of integrations for the current level of refinement
            M1_State.set_bounds_quad_refin(id_refinement);
            
            Set_Moments(&M1_State, Var_index, num_points, &Mobius_Scale_Params, struc_map_N1);
            Regularize_Moment(r_l[index_reg], &M1_State);
            
            if (id_refinement != M1_State.refinement_level) {
                cout << "id_refinement = " << id_refinement << "  " << "refinement_level = " << M1_State.refinement_level << endl;
                exit(0);
            }
            
            Set_Initial_Guess_And_Tol(M1_State);
            
            for (int i = 0; i < M1_State.NVARS; i++) {
                x[i] = M1_State.x[i];
                tol_x[i] = M1_State.tol_x[i];
            }
            
            opt = nlopt_create(NLOPT_LD_SLSQP, M1_State.NVARS); /* algorithm and dimensionality */
            nlopt_set_min_objective_orthog(opt, myfunc, Q_func_Gram_Schmidt, &M1_State);
            
            nlopt_set_maxeval(opt, max_iters_nlopt);
            
            nlopt_set_xtol_abs(opt, tol_x);
            nlopt_set_ftol_abs(opt, tol_grad);
            
            my_constraint_data data[2];
            // data[0].a = 1.0;
            // data[0].b = 1.0;
            
            // data[1].a = 1.0;
            // data[1].b = -1.0;
            
            // The constraints only need to be enforced at the quadrature points for the numerical integration
            nlopt_add_inequality_constraint_orthog(opt, myconstraint, &M1_State, 0.0);
            // nlopt_add_inequality_constraint_orthog(opt, myconstraint, &data[0], 0.0);
            // nlopt_add_inequality_constraint_orthog(opt, myconstraint, &data[1], 0.0);
            
            if (nlopt_optimize_orthog(opt, x, M1_State.I0, &minf, &min_grad, Sk_final) < 0) {
                if (id_proc == DISPLAY_ID) {
                    if (M1_State.display) {
                        printf("....................nlopt failed!.....................\n");
                    }
                }
                
                if (index_reg == max_index_reg-1 && id_refinement == M1_State.max_refinement_level) {
                    printf("....................nlopt failed: Exiting!.....................\n");
                    cout << "Failure with set of moments "  << "I0 = " << M1_State.I0 << "   " << "N1 = " << M1_State.N1 << "    " << "id = " << id_proc << "     " << "grad_norm = " << min_grad << endl;
                    exit(0);
                    goto Exiting;
                }
            } else if (min_grad > tol_grad || isnan(min_grad) || isinf(min_grad)) {
                if (M1_State.display) {
                    if (id_proc == DISPLAY_ID) {
                        printf("***********Tolerance on gradient not satisfied g(%Lf,%Lf) = %0.10Le***********\n ", x[0], x[1], min_grad);
                        cout << "x[0] + x[1] = " << x[0] + x[1] << endl;
                    }
                }
                if (index_reg == max_index_reg-1 && id_refinement == M1_State.max_refinement_level) {
                    if (id_proc == DISPLAY_ID) {
                        printf("....................nlopt failed: Exiting!.....................\n");
                    }
                    exit(0);
                    goto Exiting;
                }
            } else {
                if (M1_State.display) {
                    if (id_proc == DISPLAY_ID) {
                        printf("***********found minimum at f(%Lf,%Lf) = %0.10Le***********\n ", x[0], x[1], minf);
                        printf("***********Final gradient is min_grad = %0.10Le *********** \n", min_grad);
                        cout << "x[0] + x[1] = " << x[0] + x[1] << endl;
                    }
                }
                goto Continue;
            }
            
            if (opt != NULL) {
                nlopt_destroy(opt);
                opt = NULL;
            }
        }
    }
    
    Continue:;
    
    // Now compute Eddginton factor or partial moments based on the 
    // regime of radiation encountered
    switch (M1_State.Max_Ent_Solution_Type) {
        case CLOSING_FLUX:
            Compute_Eddington_Factor(rec_Chi2_local, index_count, Var_index, num_points, M1_State, x, Mobius_Scale_Params, struc_map_N1, Sk_final);
            if (!M1_State.flag_Taylor_Series_Expansion) {
                // We only need to compute derivatives separately in the case where we are not
                // performing a taylor series expansion
                Compute_Eddington_Factor_Derivatives(rec_Chi2_local, index_count, Var_index, num_points, M1_State, x, Mobius_Scale_Params, struc_map_N1, Sk_final);
            }
            break;
        case PARTIAL_MOMENTS:
            cout << "Partial Moments must be computed in the 3D case" << endl;
            exit(0);
            break;
        default:
            cout << "Maximum entropy solution type not specified" << endl;
            exit(0);
            break;
    }
    
    // Deallocate now
    if (opt != NULL) {
        nlopt_destroy(opt);
        opt = NULL;
    }
                
    M1_State.Deallocate_Quad_1D();
    
    Exiting:;
}

// ********************************************************************************
// This routines sets up the parameters required to solve the dual maximum-entropy
// problem for any given set of angular moments up to first-order, for gray or 
// non-gray radiation in one dimensional problems
// ********************************************************************************
void OPTIM_NON_GRAY_M1_Fixed_N1(record_Chi2 *rec_Chi2_array, record_Partial_Moments *rec_PM_array, const int &Max_Ent_Solution_Type, const M1_Var_Num_Points *num_points, const int &Problem_Type, const int &Dimension, const int &Node_Distribution_E, const int &Node_Distribution_f, const int &id_proc, const int &num_proc, const int &N_pts_Mob_Scale, const long double *Coefficients_Fit_Mob_Scale, const int &index_N1, const Algebraic_Mapping_N1 *struc_map_N1, const bool &display) {
    M1_Var_Num_Points Var_index;
    int index_count = 0;
    int NVARS = 2;
    record_Lag_Mult rec_Lag_mult;
    M1_State_Param M1_State(NVARS);
    M1_State.Problem_Type = Problem_Type;
    M1_State.Dimension = Dimension;
    M1_State.Node_Dist_E = Node_Distribution_E;
    M1_State.Node_Dist_f = Node_Distribution_f;
    M1_State.display = display;
    M1_State.Max_Ent_Solution_Type = Max_Ent_Solution_Type;
    
    
    Mobius_Scale_Parameters Mobius_Scale_Params(N_pts_Mob_Scale, Coefficients_Fit_Mob_Scale);
    
    long double x[NVARS];
    long double tol_x[NVARS];
    long double *Sk_final;
    
    Sk_final = new long double[NVARS*NVARS];
   
    for (int i = 0; i < NVARS; i++) {
        for (int j = 0; j < NVARS; j++) {
            if (i == j) {
                Sk_final[i*NVARS + j] = 1.0;
            } else {
                Sk_final[i*NVARS + j] = 0.0;
            }
        }
    }
    
    long double minf, min_grad; /* the minimum objective value, upon return */
    nlopt_opt opt = NULL;
    
    for (int i_E = 0 ; i_E < num_points->E; i_E++){
        Var_index.E = i_E;
        Var_index.f = index_N1;
        
        NLOPT_Optimization_Algorithm_NON_GRAY_M1(rec_Chi2_array, rec_PM_array, index_count, &Var_index, num_points, id_proc, M1_State, Mobius_Scale_Params, struc_map_N1, Sk_final);
        
        index_count++;
    }
    
    delete[] Sk_final;
}

// ********************************************************************************
// This routines sets up the parameters required to solve the dual maximum-entropy
// problem for any given set of angular moments up to first-order, for gray or 
// non-gray radiation in one dimensional problems
// ********************************************************************************
void OPTIM_NON_GRAY_M1_Array(record_Chi2 *rec_Chi2_global, record_Partial_Moments *rec_PM_global, const int &Max_Ent_Solution_Type, const M1_Var_Num_Points *num_points, const int &Problem_Type, const int &Dimension, const int &Node_Distribution_E, const int &Node_Distribution_f, const int &id_proc, const int &num_proc, const int &N_pts_Mob_Scale, const long double *Coefficients_Fit_Mob_Scale, const Algebraic_Mapping_N1 *struc_map_N1, const bool &display, const bool flag_use_mpi) {
    M1_Var_Num_Points Var_index;
    record_Chi2 *rec_Chi2_local, *rec_Chi2_global_temp;
    record_Partial_Moments *rec_PM_local, *rec_PM_global_temp;
    int NVARS = 2;
    
    M1_State_Param M1_State(NVARS);
    M1_State.Problem_Type = Problem_Type;
    M1_State.Dimension = Dimension;
    M1_State.Node_Dist_E = Node_Distribution_E;
    M1_State.Node_Dist_f = Node_Distribution_f;
    M1_State.Max_Ent_Solution_Type = Max_Ent_Solution_Type;
    M1_State.display = display;
    M1_State.id_proc = id_proc;
    
     // Create the rec_Chi2 datatype for MPI
    MPI_Datatype rec_Chi2_type, rec_PM_type;
    int lengths[1];
    
    switch (M1_State.Max_Ent_Solution_Type) {
        case CLOSING_FLUX:
            lengths[0] = { 13 };
            break;
        case PARTIAL_MOMENTS:
            lengths[0] = { 7 };
            cout << "Fix this " << endl;
            exit(0);
            break;
        default:
            cout << "Maximum entropy solution type not specified" << endl;
            exit(0);
            break;
    } // end switch
    
    const MPI_Aint displacements[1] = { 0 };
    MPI_Datatype types[1] = { MPI_LONG_DOUBLE };
    
    switch (M1_State.Max_Ent_Solution_Type) {
        case CLOSING_FLUX:
            MPI_Type_create_struct(1, lengths, displacements, types, &rec_Chi2_type);
            MPI_Type_commit(&rec_Chi2_type);
            break;
        case PARTIAL_MOMENTS:
            MPI_Type_create_struct(1, lengths, displacements, types, &rec_PM_type);
            MPI_Type_commit(&rec_PM_type);
            break;
        default:
            cout << "Maximum entropy solution type not specified" << endl;
            exit(0);
            break;
    } // end switch
    
    Mobius_Scale_Parameters Mobius_Scale_Params(N_pts_Mob_Scale, Coefficients_Fit_Mob_Scale);
    
    Mobius_Scale_Params.Length_Scale_Dist_Type = num_points->Length_Scale_Dist_Type;
    Mobius_Scale_Params.Least_Squares_L_Chi2_Mode = num_points->Least_Squares_L_Chi2_Mode;
    
    Finite_Diff_Parameters Finite_Diff_Chi2;
    
    long double *Sk_final;
    
    Sk_final = new long double[NVARS*NVARS];
   
    for (int i = 0; i < NVARS; i++) {
        for (int j = 0; j < NVARS; j++) {
            if (i == j) {
                Sk_final[i*NVARS + j] = 1.0;
            } else {
                Sk_final[i*NVARS + j] = 0.0;
            }
        }
    }
    
    int num_proc_per_var_E, num_proc_per_var_f;
    int id_proc_E, id_proc_f;
    
    if (Problem_Type == NON_GRAY) {
        num_proc_per_var_E = floor(sqrt(num_proc));
        num_proc_per_var_f = num_proc_per_var_E;
        id_proc_E = floor(id_proc/num_proc_per_var_f);
        id_proc_f = id_proc - id_proc_E*num_proc_per_var_f;
    } else if (Problem_Type == GRAY) {
        num_proc_per_var_E = 1;
        num_proc_per_var_f = num_proc;
        id_proc_E = 0;
        id_proc_f = id_proc;   
    }
    
    if (num_points->E % num_proc_per_var_E != 0) {
        if (id_proc == DISPLAY_ID) {
            cout << "The number num_points->E = " << num_points->E << "is not divisible by num_proc_per_var_E = " << num_proc_per_var_E << "the quotient is " << num_points->E/num_proc_per_var_E << endl;
        }
        goto Exiting;
    }
    
    if (num_points->f % num_proc_per_var_f != 0) {
        if (id_proc == DISPLAY_ID) {
            cout << "The number num_points->f = " << num_points->f << "is not divisible by num_proc_per_var_f = " << num_proc_per_var_f << "the quotient is " << num_points->f/num_proc_per_var_f << endl;
        }
        goto Exiting;
    }
    
    int incr_E = ceil(double(num_points->E)/double(num_proc_per_var_E));
    int incr_f = ceil(double(num_points->f)/double(num_proc_per_var_f));
    
    int size_rec;
    size_rec = incr_E*incr_f;
    
    switch (M1_State.Max_Ent_Solution_Type) {
        case CLOSING_FLUX:
            rec_Chi2_local = new record_Chi2[size_rec];
            break;
        case PARTIAL_MOMENTS:
            rec_PM_local = new record_Partial_Moments[size_rec];
            break;
        default:
            cout << "Maximum entropy solution type not specified" << endl;
            exit(0);
            break;
    } // end switch
        
    
    int size_rec_total;
    size_rec_total = size_rec*num_proc;
    
    if (id_proc == PRIMARY_ID) {
        switch (M1_State.Max_Ent_Solution_Type) {
            case CLOSING_FLUX:
                rec_Chi2_global_temp = new record_Chi2[size_rec_total];
                break;
            case PARTIAL_MOMENTS:
                rec_PM_global_temp = new record_Partial_Moments[size_rec_total];
                break;
            default:
                cout << "Maximum entropy solution type not specified" << endl;
                exit(0);
                break;
        } // end switch
    }
    
    int index_count = 0;
    
    for (int i_E = id_proc_E*incr_E; i_E < (id_proc_E+1)*incr_E; i_E++){
//         Finite_Diff_Chi2.flag_finite_diff_setup = 0;
        for (int i_f = id_proc_f*incr_f; i_f < (id_proc_f+1)*incr_f; i_f++) {
            
            // **********************************************************************************
            // Checks for MPI memory errors
            // **********************************************************************************
            if (index_count >= size_rec) {
                cout << "index_count exceeds size_rec ====> Out of bound memory likely to happen !!!!!!!!!!!!!!!!" << endl;
                cout << "id_proc = " << id_proc << "    " << "index_count = " << index_count << "  " << "size_rec = " << size_rec << endl;
                exit(0);
            }
            
            if (i_E < num_points->E && i_f < num_points->f) {
                Var_index.E = i_E;
                Var_index.f = i_f;
                
                NLOPT_Optimization_Algorithm_NON_GRAY_M1(rec_Chi2_local, rec_PM_local, index_count, &Var_index, num_points, id_proc, M1_State, Mobius_Scale_Params, struc_map_N1, Sk_final);
            }
            index_count++;
        }
    }
    
    // **********************************************************************************
    // Checks for MPI memory errors
    // **********************************************************************************
    if (size_rec_total != num_proc*size_rec) {
        cout << "size_rec_total and maxium size_rec do not match !!!!!!!!!!!!!!!!" << endl;
        cout << "size_rec_total = " << size_rec_total << "  " << "max size_rec = " << size_rec*num_proc << endl;
        exit(0);
    }
    
    int index, increment;
    if (flag_use_mpi) {
        MPI_Barrier(MPI_COMM_WORLD);
        
        switch (Max_Ent_Solution_Type) {
            case CLOSING_FLUX:
                MPI_Gather(rec_Chi2_local, size_rec, rec_Chi2_type, rec_Chi2_global_temp, size_rec, rec_Chi2_type, PRIMARY_ID, MPI_COMM_WORLD);
                break;
            case PARTIAL_MOMENTS:
                MPI_Gather(rec_PM_local, size_rec, rec_PM_type, rec_PM_global_temp, size_rec, rec_PM_type, PRIMARY_ID, MPI_COMM_WORLD);
                break;
            default:
                cout << "Maximum entropy solution type not specified" << endl;
                exit(0);
                break;
        }
    
        MPI_Barrier(MPI_COMM_WORLD);
        
        increment = 0;
        if (id_proc == PRIMARY_ID) {
            for (int i_proc_E = 0; i_proc_E < num_proc_per_var_E; i_proc_E++) {
                for (int id_E = 0; id_E < incr_E; id_E++) {
                    for (int i_proc_f = 0; i_proc_f < num_proc_per_var_f; i_proc_f++) {
                        for (int id_f = 0; id_f < incr_f; id_f++) {
                            index = i_proc_E*num_proc_per_var_f + i_proc_f;
                            index = index*size_rec;
                            index = index + id_E*incr_f + id_f;
                            
                            if ((i_proc_E*incr_E + id_E < num_points->E) && 
                                (i_proc_f*incr_f + id_f < num_points->f)) {   
                                
                                switch (Max_Ent_Solution_Type) {
                                    case CLOSING_FLUX:
                                        rec_Chi2_global[increment] = rec_Chi2_global_temp[index];
                                        break;
                                    case PARTIAL_MOMENTS:
                                        rec_PM_global[increment] = rec_PM_global_temp[index];
                                        break;
                                    default:
                                        cout << "Maximum entropy solution type not specified" << endl;
                                        exit(0);
                                        break;
                                }
                                increment++;
                            }
                        }
                    }
                }
            }
        }
        
        if (increment > num_points->E*num_points->f) {
            cout << "Value of increment exceeds maximum allowable value" << endl;
            cout << "increment = " << increment << " " << "num_points->E*num_points->f = " << num_points->E*num_points->f << endl;
        }
    } else {
        for (int id_E = 0; id_E < num_points->E; id_E++) {
            for (int id_f = 0; id_f < num_points->f; id_f++) {
                index = id_E*num_points->f + id_f;
                
                switch (Max_Ent_Solution_Type) {
                    case CLOSING_FLUX:
                        rec_Chi2_global[index] = rec_Chi2_local[index];
                        break;
                    case PARTIAL_MOMENTS:
                        rec_PM_global[index] = rec_PM_local[index];
                        break;
                    default:
                        cout << "Maximum entropy solution type not specified" << endl;
                        exit(0);
                        break;
                } // end switch
            }
        }
    }
    
    Exiting:;
    
    delete[] Sk_final;
    
    switch (Max_Ent_Solution_Type) {
        case CLOSING_FLUX:
            delete[] rec_Chi2_local;
            
            if (id_proc == PRIMARY_ID) {
                delete[] rec_Chi2_global_temp;
            }
            break;
        case PARTIAL_MOMENTS:
            delete[] rec_PM_local;
            
            if (id_proc == PRIMARY_ID) {
                delete[] rec_PM_global_temp;
            }
            break;
        default:
            cout << "Maximum entropy solution type not specified" << endl;
            exit(0);
            break;
    }
}

void OPTIM_NON_GRAY_M1(const M1_Var_Num_Points *num_points, const int &Max_Ent_Solution_Type, const int &Problem_Type, const int &Dimension, const int &Node_Distribution_E, const int &Node_Distribution_f, const int &id_proc, const int &num_proc, const int &N_pts_Mob_Scale, const long double *Coefficients_Fit_Mob_Scale, fstream &in_Chi2_out, const Algebraic_Mapping_N1 *struc_map_N1, const bool &display, const bool flag_use_mpi) {
    M1_Var_Num_Points Var_index;
    record_Chi2 *rec_Chi2_global;
    record_Partial_Moments *rec_PM_global;
    
    if (id_proc == PRIMARY_ID) {
        switch (Max_Ent_Solution_Type) {
            case CLOSING_FLUX:
                rec_Chi2_global = new record_Chi2[num_points->E*num_points->f];
                break;
            case PARTIAL_MOMENTS:
                rec_PM_global = new record_Partial_Moments[num_points->E*num_points->f];
                break;
            default:
                cout << "Maximum entropy solution type not specified" << endl;
                exit(0);
                break;
        }
    }
    
    OPTIM_NON_GRAY_M1_Array(rec_Chi2_global, rec_PM_global, Max_Ent_Solution_Type, num_points, Problem_Type, Dimension, Node_Distribution_E, Node_Distribution_f, id_proc, num_proc, N_pts_Mob_Scale, Coefficients_Fit_Mob_Scale, struc_map_N1, display, flag_use_mpi);
    
    int index;
    if (id_proc == PRIMARY_ID) {
        for (int id_E = 0; id_E < num_points->E; id_E++) {
            for (int id_f = 0; id_f < num_points->f; id_f++) {
                index = id_E*num_points->f + id_f;
                
                switch (Max_Ent_Solution_Type) {
                    case CLOSING_FLUX:
                        write<record_Chi2>(in_Chi2_out, rec_Chi2_global[index]);
                        break;
                    case PARTIAL_MOMENTS:
                        write<record_Partial_Moments>(in_Chi2_out, rec_PM_global[index]);
                        break;
                    default:
                        cout << "Maximum entropy solution type not specified" << endl;
                        exit(0);
                        break;
                }
            }
        }
        if (display) {
            cout << "***************************** Writing to file Completed *****************************" << endl;
        }
    }
    
    Exiting:;
    
    if (id_proc == PRIMARY_ID) {
        switch (Max_Ent_Solution_Type) {
            case CLOSING_FLUX:
                delete[] rec_Chi2_global;
                break;
            case PARTIAL_MOMENTS:
                delete[] rec_PM_global;
                break;
            default:
                cout << "Maximum entropy solution type not specified" << endl;
                exit(0);
                break;
        }
    }
}

int Check_Realizability(M1_State_Param *M1_State) {
    int flag_Realizability = 0;
    long double norm_f;
    M1_State->Boundary_Point = 0;
    M1_State->flag_Taylor_Series_Expansion = false;
    
    norm_f = pow(M1_State->N1_1, 2) + pow(M1_State->N1_2, 2) + pow(M1_State->N1_3, 2);
    norm_f = sqrt(norm_f);
    
    if (norm_f > 1.0 || M1_State->I0 < 0.0) {
        flag_Realizability = 1;
    }
    
    if (fabs(fabs(norm_f) - 1.0) < 1.0e-6) {
        M1_State->Boundary_Point = FREE_STREAMING_LIMIT;
    } else if (norm_f > 1.0 - finite_diff_h) {
//         if (M1_State->flag_finite_difference_N1) {
//             cout << "Inconsistency in finite differencing for N1 !!!!!!!!!!!!!" << endl;
//             cout << "norm_f = " << norm_f << "  " << "N1_1 = " << M1_State->N1_1 << "  " << "N1_2 = " << M1_State->N1_2 << "  " << "N1_3 = " << M1_State->N1_3 << endl;
//             exit(0);
//         } else {
//             M1_State->flag_Taylor_Series_Expansion = true;
//         }
    }
    
    return flag_Realizability;
}

long double sum_Lag_Mom(const long double *x, const int &NVARS, const long double* MOM_VEC) {
     long double sum;
     sum = 0.0;
     for (int i = 0 ; i < NVARS; i++){
         sum = sum + x[i]*MOM_VEC[i];
     }
     
     return sum;
}

long double sum_Lagrange(const long double *x, const long double *poly_Sk, const int &NVARS) {   
    long double sum = 0;
    for (int i = 0; i < NVARS; i++) {
        sum += x[i]*poly_Sk[i];
    }
    return sum;
}

void generate_polynomials_1D(long double* poly, const int &Id_angle, const M1_State_Param &M1_State) {
    poly[0] = 1.0;
    poly[1] = M1_State.x_quad[Id_angle];
//     poly[1] = sin(M1_State.x_quad[Id_angle] * PI/2.0); // sine mapping to decrease peaks near boundaries
}

long double Map_Jacob(const int &Id_angle, const M1_State_Param &M1_State) {
    long double Map_jac;
    Map_jac = (PI/2.0)*cos(M1_State.x_quad[Id_angle] * PI/2.0);
//     Map_jac = 1.0;
    return Map_jac;
}

void Regularize_Moment(const long double r_l, M1_State_Param *M1_State) {
    if (M1_State->display) {
        if (M1_State->id_proc == DISPLAY_ID) {
            printf("r_l = %Le \n", r_l);
        }
    }
    
    for (int i = 1; i < M1_State->NVARS; i++) {
        M1_State->MOM_VEC[i] = (1.0 - r_l)*M1_State->MOM_VEC[i];
    }
    M1_State->N1 = (1.0 - r_l)*M1_State->N1;
}

// long double generate_monomials_basis(const long double &Omega1, const long double &Omega2, const long double &Omega3, const int &index) {
//     long double poly;
//     //     cout << "Omega1 = " << Omega1 << "     " << "Omega2 = " << Omega2 << "     " << "Omega3 = " << Omega3 << endl;
//     
//     switch (index) {
//         case 0:
//             poly = 1.0;
//             break;
//         case 1:
//             poly = Omega1;
//             break;
//         case 2:
//             poly = Omega2;
//             break;
//         case 3:
//             poly = Omega3;
//             break;
//     };
//     
//     return poly;
// }

int obj_function(long double *F_obj, const int &NFUN, const int &Id_angle, void *fdata) {
     long double coeff_1;
     long double sum_Lag_temp;
     M1_State_Param *M1_State = (M1_State_Param *) fdata;
     
     coeff_1 = 1.0;
     
     long double *poly, *poly_Sk;
     poly = new long double[M1_State->NVARS];
     poly_Sk = new long double[M1_State->NVARS];
     
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
     
     for (int i = 0; i < M1_State->NVARS; i++) {
         poly_Sk[i] = 0.0;
         for (int j = 0; j < M1_State->NVARS; j++) {
             poly_Sk[i] += M1_State->Sk[i*M1_State->NVARS+j]*poly[j];
        }
    }

     switch (M1_State->Regime) {
         case GRAY:
             F_obj[0] = (-1.0/3.0)/pow(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS),3);
             break;
         case BOSE_EINSTEIN:
             F_obj[0] = -coeff_1*log(1.0 - exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS)));
             break;
         case HYPERBOLIC_LIMIT:
             F_obj[0] = coeff_1*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             break;
         case LOGARITHMIC_LIMIT:
             F_obj[0] = -coeff_1*log(-sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             break;
     }
     
     delete[] poly;
     delete[] poly_Sk;
     return 0;
 }
 
 int gradient(long double *grad, const int &NFUN, const int &Id_angle, void *fdata) {
     long double coeff_1;
     long double sum_Lag_temp;
     M1_State_Param *M1_State = (M1_State_Param *) fdata;
     
     coeff_1 = 1.0;
     
     long double *poly, *poly_Sk;
     poly = new long double[M1_State->NVARS];
     poly_Sk = new long double[M1_State->NVARS];
     
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
     
     for (int i = 0; i < M1_State->NVARS; i++) {
         poly_Sk[i] = 0.0;
         for (int j = 0; j < M1_State->NVARS; j++) {
             poly_Sk[i] += M1_State->Sk[i*M1_State->NVARS+j]*poly[j];
        }
    }
     
     switch (M1_State->Regime) {
        case GRAY:
            for (int i = 0;  i < M1_State->NVARS; i++) {
                grad[i] = 1.0*poly_Sk[i]/(pow(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 4));
            }
            break;
         case BOSE_EINSTEIN:
             for (int i = 0;  i < M1_State->NVARS; i++) {
                 grad[i] = coeff_1*poly_Sk[i]*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS))/(1.0 - exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS)));
             }
             break;
         case HYPERBOLIC_LIMIT:
             for (int i = 0;  i < M1_State->NVARS; i++) {
                 grad[i] = coeff_1*poly_Sk[i]*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             }
             break;
         case LOGARITHMIC_LIMIT:
             for (int i = 0;  i < M1_State->NVARS; i++) {
                 grad[i] = coeff_1*poly_Sk[i]/(-sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             }
             break;
     }
     
     delete[] poly;
     delete[] poly_Sk;
     return 0;
 }

void new_basis_Moments(const int &n, long double *MOMENTS_NEW_BASIS, const long double *MOMENTS, const long double *Sk_cur) {
    for (int i = 0; i < n; i++) {
        MOMENTS_NEW_BASIS[i] = 0.0;
        for (int j = 0; j < n; j++) {
            MOMENTS_NEW_BASIS[i] += Sk_cur[i*n+j]*MOMENTS[j];
        }
    }
} // added by jojo

long double myfunc(unsigned n, const long double *x, const long double *Sk, long double *grad, void *my_func_data) {
     M1_State_Param *M1_State = (M1_State_Param *) my_func_data;
    long double *obj_func_val, *grad_vals, *MOM_VEC_new_basis;
    long double objective_function;
    int NF_obj, NF_grad;
    obj_func_val = new long double[1];
    grad_vals = new long double[M1_State->NVARS];
    MOM_VEC_new_basis = new long double[M1_State->NVARS];
    
    NF_obj = 1;
    NF_grad = M1_State->NVARS;
    
    M1_State->set_x(x);
    M1_State->set_Sk(Sk);
    
    new_basis_Moments(M1_State->NVARS, MOM_VEC_new_basis, M1_State->MOM_VEC, Sk);
    
//     cout << "M1_State->MOM_VEC[0] = " << M1_State->MOM_VEC[0] << "    " << "MOM_VEC_new_basis[0] = " << MOM_VEC_new_basis[0] << endl;
    
    if (grad) {
        switch (M1_State->Dimension) {
         case ONE_DIMENSIONAL:
             One_Dimensional_Quadrature(grad_vals, M1_State->NVARS, gradient, M1_State);
             break;
         case THREE_DIMENSIONAL:
             Quadrature_3D_Func(grad_vals, M1_State->NVARS, gradient, M1_State);
             break;
         default:
             cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
             exit(0);
             break;
        }
        
        for (int i = 0;  i < M1_State->NVARS; i++) {
            grad[i] = grad_vals[i] - MOM_VEC_new_basis[i];
        }
    }	
    
    switch (M1_State->Dimension) {
        case ONE_DIMENSIONAL:
            One_Dimensional_Quadrature(obj_func_val, NF_obj, obj_function, M1_State);
            break;
        case THREE_DIMENSIONAL:
            Quadrature_3D_Func(obj_func_val, NF_obj, obj_function, M1_State);
            break;
        default:
            cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
            exit(0);
            break;
    }
    
    objective_function = obj_func_val[0] - sum_Lag_Mom(x,M1_State->NVARS,MOM_VEC_new_basis);
    
    delete[] obj_func_val;
    delete[] grad_vals;
    delete[] MOM_VEC_new_basis;
    
    return objective_function;
}
 
// int Q_vals_Gram_Schmidt(long double *F_Hess, const int &NFUN, const int &Id_angle, void *fdata) {
//      long double coeff_1;
//      int index_p, index_a;
//      long double sum_Lag_temp;
//      M1_State_Param *M1_State = (M1_State_Param *) fdata;
//      index_p = M1_State->index_p;
//      index_a = M1_State->index_a;
//      
//      coeff_1 = 1.0;
//     
//      long double *poly, *poly_Sk;
//      long double poly_Sk_cur, poly_ak;
//      poly = new long double[M1_State->NVARS];
//      
//      switch (M1_State->Dimension) {
//          case ONE_DIMENSIONAL:
//              generate_polynomials_1D(poly, Id_angle, *M1_State);
//              break;
//          case THREE_DIMENSIONAL:
//              generate_polynomials_3D(poly, Id_angle, *M1_State);
//              break;
//          default:
//              cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
//              exit(0);
//              break;
//      }
//      
//      poly_Sk = new long double[M1_State->NVARS];
//      for (int i = 0; i < M1_State->NVARS; i++) {
//          poly_Sk[i] = 0.0;
//          for (int j = 0; j < M1_State->NVARS; j++) {
//              poly_Sk[i] += M1_State->Sk[i*M1_State->NVARS+j]*poly[j];
//         }
//      }
//      
//      poly_Sk_cur = 0.0;
//      poly_ak = 0.0;
//      for (int i = 0; i < M1_State->NVARS; i++) {
//          poly_Sk_cur += M1_State->Sk_cur[index_p*M1_State->NVARS+i]*poly[i];
//          poly_ak += M1_State->ak[index_a*M1_State->NVARS+i]*poly[i];
//      }
//      
//      switch (M1_State->Regime) {
//          case GRAY:
//              F_Hess[0] = -4.0*poly_Sk_cur*poly_ak/(pow(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 5));
//              break;
//          case BOSE_EINSTEIN:
//              F_Hess[0] = coeff_1*poly_Sk_cur*poly_ak*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS))/pow(1.0 - exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS)), 2.0);
//             break;
//          case HYPERBOLIC_LIMIT:
//              F_Hess[0] = coeff_1*poly_Sk_cur*poly_ak*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
//              break;
//          case LOGARITHMIC_LIMIT:
//              F_Hess[0] = coeff_1*poly_Sk_cur*poly_ak/pow(-sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 2.0);
//             break;
//      }
//      
// //      cout << "Map_Jacob(Id_angle, *M1_State) = " << Map_Jacob(Id_angle, *M1_State) << endl;
//      
// //      F_Hess[0] *= Map_Jacob(Id_angle, *M1_State);
//      
//      delete[] poly;
//      delete[] poly_Sk;
//      return 0;
//  }
 
 long double dIn_dInm1(const int &i, const int &j, const long double &I0, const long double &N1_1, const long double &N1_2, const long double &N1_3) {
     long double temp_val;
     if (i == j) {
         temp_val = 1.0;
     } else {
         temp_val = 0.0;
     }
     return temp_val;
 }
 
 void Q_func_Gram_Schmidt(unsigned n, const long double *x, const long double *Sk, long double *Q_data, void *my_func_data, void *f_data_MN_State) {
     
     M1_State_Param *M1_State_temp = (M1_State_Param *) f_data_MN_State; 
     NLOPT_State *NLP_State = (NLOPT_State *) my_func_data;  
     M1_State_temp->set_x(x);
     M1_State_temp->set_Sk(Sk);
     M1_State_temp->set_Sk_cur(NLP_State->Sk_cur);
     M1_State_temp->set_ak(NLP_State->ak);
     M1_State_temp->index_p = NLP_State->index_p;
     M1_State_temp->index_a = NLP_State->index_a;
     
     switch (M1_State_temp->Dimension) {
         case ONE_DIMENSIONAL:
             One_Dimensional_Quadrature(Q_data, 1, H_ONE_Matrix_Entries, M1_State_temp);
             break;
         case THREE_DIMENSIONAL:
             Quadrature_3D_Func(Q_data, 1, H_ONE_Matrix_Entries, M1_State_temp);
             break;
         default:
             cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
             exit(0);
             break;
     }
     
//      if (NLP_State->index_p == NLP_State->index_a) {
//          Q_data[0] = 1.0;
//     } else {
//         Q_data[0] = 0.0;
//     }
     
 }

void Set_Moments(M1_State_Param *M1_State, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, Mobius_Scale_Parameters *Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1) {
    long double f_test;
    long double r_N1;
    
    if (M1_State->flag_finite_difference_N1) {
        if (M1_State->finite_difference_type_N1 == TAYLOR_SERIES_EXPANSION || 
            M1_State->finite_difference_type_N1 == DERIVATIVE_BOUNDARY) {
            f_test = M1_State->f_test_knot + M1_State->x_finite_diff_N1;
        } else {
            cout << "Finite difference type for N1 not specified!!!!!!!!!" << endl;
            exit(0);
        }
    } else {
        if (M1_State->Node_Dist_f == UNIFORM_DISTRIBUTION) {
            f_test = Uniform_Distribution(Var_index->f, num_points->f - 1, 0.0, 1.0);
        } else if (M1_State->Node_Dist_f == UNIFORM_ALGEBRAIC_MAPPING_N1) {
            r_N1 = Uniform_Distribution(Var_index->f, num_points->f - 1, -1.0, 1.0);
            // f_test = Inverse_mapping_L_Chi2(struc_map_N1->f_L_Chi2, r_N1);
        } else if (M1_State->Node_Dist_f == CHEBYSHEV_ALGEBRAIC_MAPPING_N1) {
            r_N1 = zeros_shifted( Var_index->f, num_points->f, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
            // f_test = Inverse_mapping_L_Chi2(struc_map_N1->f_L_Chi2, r_N1);
            // cout << "r_N1 = " << r_N1 << "  "  << "f_test = " << f_test << endl;
        } else {
            if (num_points->f == 1) {
                f_test = 0.0;  
            } else {
                // if we want to construct a polynomial of degree n, we need to consider
                // the root of the (n+1)^th order polynomial
                //             f_test = zeros_shifted( Var_index->f, num_points->f, 0.0, 1.0, M1_State->Node_Dist_f);
                f_test = zeros_shifted( Var_index->f + num_points->f - 1, 2*(num_points->f - 1) + 1, -1.0, 1.0, M1_State->Node_Dist_f);
            }
        }
    }
    
    if (fabs(f_test) < 1.0e-8) {
        f_test = 0.005;
    }
    
    if (M1_State->flag_finite_difference_r_I0) {
        if (M1_State->finite_difference_type_r_I0 == DERIVATIVE_BOUNDARY_RADIATIVE_ENERGY) {
            M1_State->ratio_I0 = M1_State->ratio_I0_knot + M1_State->x_finite_diff_r_I0;
        } else {
            cout << "Finite difference type for exponential map of I0 not specified!!!!!!!!!" << endl;
            exit(0);
        }
    } else {
        if (M1_State->Node_Dist_E == DISCRETE_SET_ENERGY) {
            if (Var_index->E == 0) {
                M1_State->ratio_I0 = -1.0;
            } else if (Var_index->E == num_points->E - 1) {
                M1_State->ratio_I0 = 1.0;
            } else {
                M1_State->ratio_I0 = 0.0;
            }
            M1_State->MOM_VEC[0] = E_vals[Var_index->E];
        } else {
            if (M1_State->Node_Dist_E == UNIFORM_DISTRIBUTION) {
//                 if (Var_index->E < 10) {
//                     M1_State->ratio_I0 = Uniform_Distribution(Var_index->E, 10 - 1, -1.0, -0.9);
//                 } else if (Var_index->E < num_points->E - 10) {
//                     M1_State->ratio_I0 = Uniform_Distribution(Var_index->E - 10, num_points->E - 20 - 1, -0.9, 0.9);
//                 } else {
//                     M1_State->ratio_I0 = Uniform_Distribution(Var_index->E - (num_points->E - 10), 10 - 1, 0.9, 1.0);
//                 }
                M1_State->ratio_I0 = Uniform_Distribution(Var_index->E, num_points->E - 1, -1.0, 1.0);
            } else {
                M1_State->ratio_I0 = zeros_shifted(Var_index->E, num_points->E, -1.0, 1.0, M1_State->Node_Dist_E);
//                 M1_State->ratio_I0 = zeros_shifted( Var_index->E + num_points->E - 1, 2*(num_points->E - 1) + 1, -1.0, 1.0, M1_State->Node_Dist_E);
//                 M1_State->ratio_I0 = 2.0*M1_State->ratio_I0 - 1.0;
            }
        }
    }
    
    // Now compute the energy density based on the type of mapping and nodal distribution adopted
    if (M1_State->Node_Dist_E != DISCRETE_SET_ENERGY) {
        if (M1_State->flag_finite_difference_N1) {
            // In the case where we are performing finite differencing, the other
            // variables must remain constant
            M1_State->MOM_VEC[0] = Inverse_Mobius_Transformation(M1_State->ratio_I0, Mobius_Scale_Params->Evaluate_Length_Scale(f_test));
            // M1_State->MOM_VEC[0] = Inverse_Mobius_Transformation(M1_State->ratio_I0, Mobius_Scale_Params->Evaluate_Length_Scale(M1_State->f_test_knot));
        } else {
            if (M1_State->Node_Dist_f == UNIFORM_ALGEBRAIC_MAPPING_N1 || M1_State->Node_Dist_f == CHEBYSHEV_ALGEBRAIC_MAPPING_N1) {
                // M1_State->MOM_VEC[0] = Inverse_Mobius_Transformation(M1_State->ratio_I0, Mobius_Scale_Params->Evaluate_Length_Scale_r_N1(r_N1));
            } else if (struc_map_N1->flag_Algebraic_Map_N1) {
                // r_N1 = Mapping_L_Chi2(struc_map_N1->f_L_Chi2, f_test);
                // M1_State->MOM_VEC[0] = Inverse_Mobius_Transformation(M1_State->ratio_I0, Mobius_Scale_Params->Evaluate_Length_Scale_r_N1(r_N1));
            } else {
                M1_State->MOM_VEC[0] = Inverse_Mobius_Transformation(M1_State->ratio_I0, Mobius_Scale_Params->Evaluate_Length_Scale(f_test));
            }
        }
    }
    
    if (num_points->E == 1 && M1_State->Problem_Type == NON_GRAY) {
        switch (num_points->Maximum_Entropy_Solution_Regime) {
            case HYPERBOLIC_LIMIT:
                M1_State->ratio_I0 = -1.0;
                break;
            case LOGARITHMIC_LIMIT:
                M1_State->ratio_I0 = 1.0;
                break;
            default:
                cout << "Invalid value for Maximum_Entropy_Solution_Regime" << endl;
                exit(0);
                break;
        }
    }
    
    Set_Regime(M1_State);
    
    M1_State->I0 = M1_State->MOM_VEC[0];
    M1_State->N1 = f_test;
    
    M1_State->N1_1 = M1_State->N1;
    M1_State->N1_2 = 0.0;
    M1_State->N1_3 = 0.0;
    
    M1_State->flag_Realizability = Check_Realizability(M1_State);
    
    M1_State->MOM_VEC[1] = f_test*M1_State->MOM_VEC[0];
    
    if (M1_State->display) {
        if (M1_State->id_proc == DISPLAY_ID) {
            cout << "refinement_level = " << M1_State->refinement_level << "  " << "ratio_I0 = " << M1_State->ratio_I0 << "    "  << "I0 = " << M1_State->I0 << "   " << "N1 = " << M1_State->N1 << endl; 
        }
    }
}

void Set_Initial_Guess_And_Tol(M1_State_Param &M1_State) {
    switch (M1_State.Regime){
        case GRAY:
            M1_State.x[0] = -pow(4.0*PI/M1_State.MOM_VEC[0], 1.0/4.0);
            break;
        case HYPERBOLIC_LIMIT:
            M1_State.x[0] = -log(4.0*PI/M1_State.MOM_VEC[0]);
            break;
        case BOSE_EINSTEIN:
            M1_State.x[0] = -log(1.0 + 4.0*PI/M1_State.MOM_VEC[0]);
            break;
        case LOGARITHMIC_LIMIT:
            M1_State.x[0] = -4.0*PI/M1_State.MOM_VEC[0];
            break;
        default:
            cout << "Invalid Regime type !!!!!!!!!!" << endl;
            exit(0);
            break;
    };
    
    M1_State.x[0] = roundn(M1_State.x[0], 6);
    
    if (M1_State.x[0] == 0.0) {
        M1_State.x[0] = -1.0;
    }
    
    for (int i = 1; i < M1_State.NVARS; i++) {
        M1_State.x[i] = 0.0;  /* some initial guess */
    }
    
    for (int i = 0; i < M1_State.NVARS; i++) {
        M1_State.tol_x[i] = 1e-32;
    }
}

void Set_Regime(M1_State_Param *M1_State) {
    if (M1_State->Problem_Type == NON_GRAY) {
        if (M1_State->MOM_VEC[0] < I0_val_max) {
            if (M1_State->ratio_I0 < -1.0 + 1.0e-6) {
                M1_State->Regime = HYPERBOLIC_LIMIT;
                if (M1_State->display) {
                    if (M1_State->id_proc == DISPLAY_ID) {
                        cout << "Hyperbolic Limit" << endl;
                    }
                }
                M1_State->MOM_VEC[0] = 1.0e0;
                M1_State->ratio_I0 = -1.0;
            } else if (M1_State->ratio_I0 < 1.0 - 1.0e-6) {
                M1_State->Regime = BOSE_EINSTEIN;
                if (M1_State->display) {
                    if (M1_State->id_proc == DISPLAY_ID) {
                        cout << "Bose Einstein statistics" << endl;
                    }
                }
            } else {
                // M1_State->Regime = LOGARITHMIC_LIMIT; 
                M1_State->Regime = BOSE_EINSTEIN;
                if (M1_State->display) {
                    if (M1_State->id_proc == DISPLAY_ID) {
                        cout << "Logarithmic Limit" << endl;
                    }
                }
                M1_State->ratio_I0 = 1.0;
                // M1_State->MOM_VEC[0] = 1.0e0;
                M1_State->MOM_VEC[0] = I0_val_max;
            }
        } else {
            // M1_State->Regime = LOGARITHMIC_LIMIT; 
            M1_State->Regime = BOSE_EINSTEIN;
            if (M1_State->display) {
                if (M1_State->id_proc == DISPLAY_ID) {
                    cout << "Logarithmic Limit" << endl;
                }
            }
            M1_State->ratio_I0 = 1.0;
            // M1_State->MOM_VEC[0] = 1.0e0;
            M1_State->MOM_VEC[0] = I0_val_max;
        }
    } else if (M1_State->Problem_Type == GRAY) {
        M1_State->Regime = GRAY;
        M1_State->MOM_VEC[0] = 1.0;
    } else {
        if (M1_State->display) {
            cout << "Problem Type not specified" << endl;
        }
    }
    
    // M1_State->Regime = BOSE_EINSTEIN;
    // M1_State->MOM_VEC[0] = 1.0e1;
}

// ******************************************************************************
// This routine computes entries of the Hessian matrix of second-derivatives
// of the objective function (or first-derivatives of the angular moments)
// ===> H^{(1)}_{ij} = {\partial I^{(i)}} {\partial \lam_j}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
int H_ONE_Matrix_Entries(long double *H_ONE, const int &NFUN, const int &Id_angle, void *fdata) {
     long double coeff_1;
     int index_p, index_a;
     M1_State_Param *M1_State = (M1_State_Param *) fdata;
     index_p = M1_State->index_p;
     index_a = M1_State->index_a;
     coeff_1 = 1.0;
    
     long double *poly, *poly_Sk;
     long double poly_Sk_cur, poly_ak;
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
     
     poly_Sk_cur = 0.0;
     poly_ak = 0.0;
     for (int i = 0; i < M1_State->NVARS; i++) {
         poly_Sk_cur += M1_State->Sk_cur[index_p*M1_State->NVARS+i]*poly[i];
         poly_ak += M1_State->ak[index_a*M1_State->NVARS+i]*poly[i];
     }
     
     long double summation;
     
     switch (M1_State->Regime) {
         case GRAY:
             if (sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS) > 0.0) {
                 cout << "sum_Lagrange = " << sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS) << endl; 
             }
             H_ONE[0] = -4.0*poly_Sk_cur*poly_ak/(pow(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 5));
             break;
         case BOSE_EINSTEIN:
             // use this to avoid overflow for the exponential function
             summation = sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS);
             H_ONE[0] = coeff_1*poly_Sk_cur*poly_ak*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS))/pow((1.0 - exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS))), 2.0);
             break;
         case HYPERBOLIC_LIMIT:
             H_ONE[0] = coeff_1*poly_Sk_cur*poly_ak*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             break;
         case LOGARITHMIC_LIMIT:
             H_ONE[0] = coeff_1*poly_Sk_cur*poly_ak/pow(-sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 2.0);
            break;
     }
     
     delete[] poly;
     delete[] poly_Sk;
     return 0;
 }

// ******************************************************************************
// This routine computes entries of the second-derivatives of the angular moments
// with respect to the Lagrange multipliers
// ===> H^{(2)}_{ijk} = {\partial^{2} I^{(i)}} {\partial \lam_j\partial \lam_k}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
int H_TWO_Matrix_Entries(long double *H_TWO, const int &NFUN, const int &Id_angle, void *fdata) {
    long double coeff_1;
     int index_Lag1, index_Lag2, index_Moment;
     long double *poly, *poly_Sk;
     long double poly_Sk_cur, poly_Lag_1, poly_Lag_2;
     M1_State_Param *M1_State = (M1_State_Param *) fdata;
     coeff_1 = 1.0;
     long double exp_x;
     index_Moment = M1_State->index_Moment;
     index_Lag1 = M1_State->index_Lag_i;
     index_Lag2 = M1_State->index_Lag_j;
    
     poly = new long double[M1_State->NVARS];
     poly_Sk = new long double[M1_State->NVARS];
     
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
     
     for (int i = 0; i < M1_State->NVARS; i++) {
         poly_Sk[i] = 0.0;
         for (int j = 0; j < M1_State->NVARS; j++) {
             poly_Sk[i] += M1_State->Sk[i*M1_State->NVARS+j]*poly[j];
        }
     }
     
     poly_Sk_cur = 0.0;
     poly_Lag_1 = 0.0;
     poly_Lag_2 = 0.0;
     for (int i = 0; i < M1_State->NVARS; i++) {
         poly_Sk_cur += M1_State->Sk_cur[index_Moment*M1_State->NVARS+i]*poly[i];
         poly_Lag_1 += M1_State->ak[index_Lag1*M1_State->NVARS+i]*poly[i];
         poly_Lag_2 += M1_State->ak[index_Lag2*M1_State->NVARS+i]*poly[i];
     }
     
     switch (M1_State->Regime) {
         case GRAY:
             H_TWO[0] = 20.0*poly_Sk_cur*poly_Lag_1*poly_Lag_2/(pow(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 6));
             break;
         case BOSE_EINSTEIN:
             exp_x = exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             H_TWO[0] = coeff_1*poly_Sk_cur*poly_Lag_1*poly_Lag_2*exp_x*(1.0 + exp_x)/pow((1.0 - exp_x), 3.0);
             break;
         case HYPERBOLIC_LIMIT:
             H_TWO[0] = coeff_1*poly_Sk_cur*poly_Lag_1*poly_Lag_2*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             break;
         case LOGARITHMIC_LIMIT:
             H_TWO[0] = 2.0*coeff_1*poly_Sk_cur*poly_Lag_1*poly_Lag_2/pow(-sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 3.0);
            break;
     }
     
     delete[] poly;
     delete[] poly_Sk;
     return 0;
}
 
// ******************************************************************************
// This routine computes entries of the third-derivatives of the angular moments
// with respect to the Lagrange multipliers
// ===> H^{(3)}_{ijkl} = {\partial^{3} I^{(i)}} {\partial \lam_j \partial \lam_k \partial \lam_l}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
int H_THREE_Matrix_Entries(long double *H_THREE, const int &NFUN, const int &Id_angle, void *fdata) {
     long double coeff_1, exp_x;
     int index_Lag1, index_Lag2, index_Lag3, index_Moment;
     M1_State_Param *M1_State = (M1_State_Param *) fdata;
     index_Moment = M1_State->index_Moment;
     index_Lag1 = M1_State->index_Lag_i;
     index_Lag2 = M1_State->index_Lag_j;
     index_Lag3 = M1_State->index_Lag_k;
     coeff_1 = 1.0;
    
     long double *poly, *poly_Sk;
     long double poly_Sk_cur, poly_Lag_1, poly_Lag_2, poly_Lag_3;
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
     
     poly_Sk_cur = 0.0;
     poly_Lag_1 = 0.0;
     poly_Lag_2 = 0.0;
     poly_Lag_3 = 0.0;
     for (int i = 0; i < M1_State->NVARS; i++) {
         poly_Sk_cur += M1_State->Sk_cur[index_Moment*M1_State->NVARS+i]*poly[i];
         poly_Lag_1 += M1_State->ak[index_Lag1*M1_State->NVARS+i]*poly[i];
         poly_Lag_2 += M1_State->ak[index_Lag2*M1_State->NVARS+i]*poly[i];
         poly_Lag_3 += M1_State->ak[index_Lag3*M1_State->NVARS+i]*poly[i];
     }
     
     switch (M1_State->Regime) {
         case GRAY:
             H_THREE[0] = -120.0*poly_Lag_1*poly_Lag_2*poly_Lag_3*poly_Sk_cur/(pow(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 7));
             break;
         case BOSE_EINSTEIN:
             exp_x = exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             H_THREE[0] = coeff_1*poly_Lag_1*poly_Lag_2*poly_Lag_3*poly_Sk_cur*exp_x*(1.0 + 4.0*exp_x + exp_x*exp_x)/pow((1.0 - exp_x), 4.0);
             break;
         case HYPERBOLIC_LIMIT:
             H_THREE[0] = coeff_1*poly_Lag_1*poly_Lag_2*poly_Lag_3*poly_Sk_cur*exp(sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS));
             break;
         case LOGARITHMIC_LIMIT:
             H_THREE[0] = 6.0*coeff_1*poly_Lag_1*poly_Lag_2*poly_Lag_3*poly_Sk_cur/pow(-sum_Lagrange(M1_State->x,poly_Sk,M1_State->NVARS), 4.0);
            break;
     }
     
     delete[] poly;
     delete[] poly_Sk;
     return 0;
 } 
 
// ******************************************************************************
// This routine computes the Hessian matrix of first-derivatives of the 
// angular moments with respect to the Lagrange multipliers
// ===> H^{(1)}_{ij} = {\partial I^{(i)}} {\partial \lam_j}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
 void Compute_H_ONE_Matrix(long double *H_ONE, const long double *x, const long double *Sk, M1_State_Param &M1_State) {
     long double H_ONE_ij;
     M1_State.set_x(x);
     M1_State.set_Sk(Sk);
     M1_State.set_Sk_cur(Sk);
     M1_State.set_ak(Sk);
     
     for (int i = 0; i < M1_State.NVARS; i++) {
         for (int j = 0; j < M1_State.NVARS; j++) {
             M1_State.index_p = i;
             M1_State.index_a = j;
             
             switch (M1_State.Dimension) {
                 case ONE_DIMENSIONAL:
                     One_Dimensional_Quadrature(&H_ONE_ij, 1, H_ONE_Matrix_Entries, &M1_State);
                     break;
                 case THREE_DIMENSIONAL:
                     Quadrature_3D_Func(&H_ONE_ij, 1, H_ONE_Matrix_Entries, &M1_State);  
                     break;
                 default:
                     cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
                     exit(0);
                     break;
            }      
            
            H_ONE[i*M1_State.NVARS + j] = H_ONE_ij;
         }
     }
 }
 
// ******************************************************************************
// This routine computes the matrix of second-derivatives of the angular moments
// with respect to the Lagrange multipliers
// ===> H^{(2)}_{ijk} = {\partial^{2} I^{(i)}} {\partial \lam_j\partial \lam_k}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
 void Compute_H_TWO_Matrix(long double *H_TWO, const long double *x, const long double *Sk, M1_State_Param &M1_State) {
     long double H_TWO_ijk;
     M1_State.set_x(x);
     M1_State.set_Sk(Sk);
     M1_State.set_Sk_cur(Sk);
     M1_State.set_ak(Sk);
     
     for (int i = 0; i < M1_State.NVARS; i++) {
         for (int j = 0; j < M1_State.NVARS; j++) {
             for (int k = 0; k < M1_State.NVARS; k++) {
                 M1_State.index_Moment = i;
                 M1_State.index_Lag_i = j;
                 M1_State.index_Lag_j = k;
                 
                 switch (M1_State.Dimension) {
                     case ONE_DIMENSIONAL:
                         One_Dimensional_Quadrature(&H_TWO_ijk, 1, H_TWO_Matrix_Entries, &M1_State);
                         break;
                     case THREE_DIMENSIONAL:
                         Quadrature_3D_Func(&H_TWO_ijk, 1, H_TWO_Matrix_Entries, &M1_State);  
                         break;
                     default:
                         cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
                         exit(0);
                         break;
                }  
                H_TWO[(i*M1_State.NVARS + j)*M1_State.NVARS + k] = H_TWO_ijk;
             }
         }
     } 
 }
 
 void transpose_matrix(const int &nl, const int &nc, long double *A_transpose, const long double *A) {
     // This routine aims at storing elements of the transpose of 
     // matrix A in the matrix A_transpose
     
     if (nl != nc) {
         cout << "Cannot take transpose of matrix as it is not square" << endl; 
         exit(0);
     }
     
     for (int i = 0; i < nl; i++) {
         for (int j = 0; j < nc; j++) {
             A_transpose[i*nc + j] = A[j*nc + i];
         }
     }
 }
 
 void Backward_Substitution(const int &nl, const int &nc, long double *A, long double *Q, const long double *R) {
     // This routine aims at solving the system of linear equations: R Q = A for Q
     // where R is an upper triangular matrix
     // This can be achieved using backward substitution
     
     // nl: number of lines
     // nc: number of columns
     
     /* Backward substitution for discovering values of unknowns */
     for(int k = 0; k < nc; k++) {                   
         for(int i = nl-1; i >= 0; i--) {                     
             Q[i*nc + k] = A[i*nc + k];
             for(int j = i+1; j < nl;j++) {
                 if(i != j) {
                     Q[i*nc + k] = Q[i*nc + k] - R[i*nl + j]*Q[j*nc + k];
                }          
            }
            Q[i*nc + k] = Q[i*nc + k]/R[i*nl + i];
        }
     }
 }
 
// ******************************************************************************
// This routine computes the matrix of first-derivatives of the Lagrange 
// multipliers with respect to the angular moments (which is also the inverse of
// the Hessian matrix)
// ===> A^{(1)}_{ij} = {\partial \lam_i} {\partial I^{(j)}}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
 void Compute_A_ONE_Matrix(long double *A_ONE, const long double *x, const long double *Sk, M1_State_Param &M1_State) {
     long double *Q, *R, *H_ONE;
     H_ONE = new long double[M1_State.NVARS*M1_State.NVARS];
     Q = new long double[M1_State.NVARS*M1_State.NVARS];
     R = new long double[M1_State.NVARS*M1_State.NVARS];
     
     // Compute Hessian Matrix H_ONE
     Compute_H_ONE_Matrix(H_ONE, x, Sk, M1_State);
     
     
//      cout << "************************** Original !!!!!!!! ***********************" << endl;
// 
//      for (int i = 0; i < M1_State.NVARS; i++) {
//          for (int j = 0; j < M1_State.NVARS; j++) {
//              cout << "i = " << i << "   " << "j = " << j << "   " << "H_ONE = " << H_ONE[i*M1_State.NVARS+j] << endl;
//         }
//     }
     
     // ******************************************************************************
     // Perform QR decomposition on H_ONE using Modified Gramd Schmidt Factorization
     // ==> H_ONE = Q R
     // where R is upper triangular
     // ******************************************************************************
     Modified_Gram_Schmidt_Factorization(M1_State.NVARS, H_ONE, Q, R);
     
//      for (int i = 0; i < M1_State.NVARS; i++) {
//          for (int j = 0; j < M1_State.NVARS; j++) {
//              cout << "i = " << i << "   " << "j = " << j << "   " << "Q = " << Q[i*M1_State.NVARS+j] << "   " << "R = " << R[i*M1_State.NVARS+j] << endl;
//         }
//     }
     
     // ******************************************************************************
     // The inverse of H_ONE satisfies: 
     // ==> H_ONE inv(H_ONE) = \delta_{ij}
     // Using the QR decomposition of H_ONE, we then have
     // Q R inv(H_ONE) = \delta_{ij}
     // This system is solved in two steps here:
     // First : Q A = \delta_{ij} where A = R inv(H_ONE)
     // And then : R inv(H_ONE) = A
     // ******************************************************************************
     
     
     // ******************************************************************************
     // Here we solve Q A = \delta_{ij} ==> A = inv(Q) \delta_{ij} = Q^{T} \delta_{ij}
     // In fact inv(Q) = Q^{T} since Q is an orthogonal matrix
     // ******************************************************************************
     transpose_matrix(M1_State.NVARS, M1_State.NVARS, H_ONE, Q); // now A ==> H_ONE
     
     // ******************************************************************************
     // Now we solve : R inv(H_ONE) = A using backward substitution since
     // R is an upper triangular matrix
     // ******************************************************************************
     Backward_Substitution(M1_State.NVARS, M1_State.NVARS, H_ONE, A_ONE, R);
     
//      for (int i = 0; i < M1_State.NVARS; i++) {
//          for (int j = 0; j < M1_State.NVARS; j++) {
//              cout << "i = " << i << "   " << "j = " << j << "   " << "A_ONE = " << A_ONE[i*M1_State.NVARS+j] << endl;
//         }
//     }
     
     delete[] H_ONE;
     delete[] Q;
     delete[] R;
 }
 
// ******************************************************************************
// This routine computes the matrix of second-derivatives of the Lagrange 
// multipliers with respect to the angular moments
// ===> A^{(2)}_{ijk} = {\partial^{2} \lam_i} {\partial I^{(j)} \partial I^{(k)}}
// where \lam_j are the Lagrange multipliers
// ******************************************************************************
 void Compute_A_TWO_Matrix(long double *A_TWO, long double *A_ONE, const int &j, const int &l, M1_State_Param &M1_State) {
     long double H_TWO_imk;
     long double *temp_mat_M_ijl;
     temp_mat_M_ijl = new long double[M1_State.NVARS];
     
     // A^{(1)}_{ij} = Inv_Hessian_{ij}
     
     for (int i = 0; i < M1_State.NVARS; i++) {
         M1_State.index_Moment = i;
         temp_mat_M_ijl[i] = 0.0;
         for (int m = 0; m < M1_State.NVARS; m++) {
             for (int k = 0; k < M1_State.NVARS; k++) {
                 M1_State.index_Lag_i = m;
                 M1_State.index_Lag_j = k;
                 // Compute H^{(2)}_{imk}
                 
                 switch (M1_State.Dimension) {
                     case ONE_DIMENSIONAL:
                         One_Dimensional_Quadrature(&H_TWO_imk, 1, H_TWO_Matrix_Entries, &M1_State);
                         break;
                     case THREE_DIMENSIONAL:
                         Quadrature_3D_Func(&H_TWO_imk, 1, H_TWO_Matrix_Entries, &M1_State);  
                         break;
                     default:
                         cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
                         exit(0);
                         break;
                } 
                
                // Compute M_{ijl} =  H^{(1)}_{ip} A^{(2)}_{pjl} = - H^{(2)}_{imk} A^{(1)}_{mj} A^{(1)}_{kl}
                temp_mat_M_ijl[i] += -H_TWO_imk*A_ONE[m*M1_State.NVARS+j]*A_ONE[k*M1_State.NVARS+l];
             }
         }
     }
     
     for (int k = 0; k < M1_State.NVARS; k++) {
         A_TWO[k] = 0.0;
         // Now compute A^{(2)}_{kjl} = inv(H^{(1)})_{ki} M_{ijl}
         for (int i = 0; i < M1_State.NVARS; i++) {
             A_TWO[k] += temp_mat_M_ijl[i]*A_ONE[k*M1_State.NVARS+i];
         }
     }
     
     delete[] temp_mat_M_ijl;
 }
  
void Compute_A_TWO_Matrix(long double *A_TWO, long double *A_ONE, M1_State_Param &M1_State) {
     long double H_TWO_ipq;
     long double *temp_mat_M_ijl;
     int index_temp;
     temp_mat_M_ijl = new long double[M1_State.NVARS*M1_State.NVARS*M1_State.NVARS];
     
     // A^{(1)}_{ij} = Inv_Hessian_{ij}
     
     for (int i = 0; i < M1_State.NVARS; i++) {
         M1_State.index_Moment = i;
         for (int j = 0; j < M1_State.NVARS; j++) {
             for (int k = 0; k < M1_State.NVARS; k++) {
                 index_temp = (i*M1_State.NVARS + j)*M1_State.NVARS + k;
                 temp_mat_M_ijl[index_temp] = 0.0;
                 for (int p = 0; p < M1_State.NVARS; p++) {
                     for (int q = 0; q < M1_State.NVARS; q++) {
                         M1_State.index_Lag_i = p;
                         M1_State.index_Lag_j = q;
                         // Compute H^{(2)}_{ipq}
                         switch (M1_State.Dimension) {
                             case ONE_DIMENSIONAL:
                                 One_Dimensional_Quadrature(&H_TWO_ipq, 1, H_TWO_Matrix_Entries, &M1_State);
                                 break;
                             case THREE_DIMENSIONAL:
                                 Quadrature_3D_Func(&H_TWO_ipq, 1, H_TWO_Matrix_Entries, &M1_State);  
                                 break;
                             default:
                                 cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
                                 exit(0);
                                 break;
                        } 
                         // Compute M_{ijk} =  H^{(1)}_{ip} A^{(2)}_{pjk} = - H^{(2)}_{ipq} A^{(1)}_{pj} A^{(1)}_{qk}
                         temp_mat_M_ijl[index_temp] += -H_TWO_ipq*A_ONE[p*M1_State.NVARS+j]*A_ONE[q*M1_State.NVARS+k];
                    }
                 }
             }
         }
     }
     
     int index_temp_2;
     for (int i = 0; i < M1_State.NVARS; i++) {
         for (int j = 0; j < M1_State.NVARS; j++) {
             for (int k = 0; k < M1_State.NVARS; k++) {
                 index_temp = (i*M1_State.NVARS + j)*M1_State.NVARS + k;
                 A_TWO[index_temp] = 0.0;
                 // Now compute A^{(2)}_{ijk} = inv(H^{(1)})_{ip} M_{pjk}
                 for (int p = 0; p < M1_State.NVARS; p++) {
                     index_temp_2 = (p*M1_State.NVARS + j)*M1_State.NVARS + k;
                     A_TWO[index_temp] += A_ONE[i*M1_State.NVARS+p]*temp_mat_M_ijl[index_temp_2];
                }
            }
         }
     }
     
     delete[] temp_mat_M_ijl;
 }
 
// *************************************************************************************************
// This routine computes the matrix of third-derivatives of the Lagrange 
// multipliers with respect to the angular moments
// ===> A^{(2)}_{ijkl} = {\partial^{2} \lam_i} {\partial I^{(j)} \partial I^{(k)} \partial I^{(l)}}
// where \lam_j are the Lagrange multipliers
// *************************************************************************************************
 void Compute_A_THREE_Matrix(long double *A_THREE, const long double *A_ONE, const long double *A_TWO, const int &j, const int &k, const int &l, const int &p, M1_State_Param &M1_State) {
     long double d2_In_lower_dlambda_2, d3_In_lower_dlambda_3;
     long double *temp_mat_V_ijkl;
     temp_mat_V_ijkl = new long double[M1_State.NVARS];
     int index_temp_1, index_temp_2, index_temp_3;
     
     // A^{(1)}_{ij} = Inv_Hessian_{ij}
     
     for (int i = 0; i < M1_State.NVARS; i++) {
         M1_State.index_Moment = i;
         temp_mat_V_ijkl[i] = 0.0;
         for (int m = 0; m < M1_State.NVARS; m++) {
             for (int p = 0; p < M1_State.NVARS; p++) {
                 for (int q = 0; q < M1_State.NVARS; q++) {
                     M1_State.index_Lag_i = m;
                     M1_State.index_Lag_j = p;
                     M1_State.index_Lag_k = q; 
                     
                     // Compute H^{(3)}_{impq}
                     switch (M1_State.Dimension) {
                         case ONE_DIMENSIONAL:
                             One_Dimensional_Quadrature(&d3_In_lower_dlambda_3, 1, H_THREE_Matrix_Entries, &M1_State);
                             break;
                         case THREE_DIMENSIONAL:
                             Quadrature_3D_Func(&d3_In_lower_dlambda_3, 1, H_THREE_Matrix_Entries, &M1_State);  
                             break;
                         default:
                             cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
                             exit(0);
                             break;
                    } 
                     // Compute M_{ijkl} = - H^{(3)}_{impq} A^{(1)}_{mj} A^{(1)}_{pk} A^{(1)}_{ql}
                     
                     temp_mat_V_ijkl[i] += -d3_In_lower_dlambda_3*A_ONE[m*M1_State.NVARS+j]*A_ONE[p*M1_State.NVARS+k]*A_ONE[q*M1_State.NVARS+l];
                 }
             }
         }
     }
     
     // Compute the matrix A^{(2)}: a third-order tensor
     // Compute_A_TWO_Matrix(A_TWO, A_ONE, M1_State);
         
     for (int i = 0; i < M1_State.NVARS; i++) {
         M1_State.index_Moment = i;
         
         for (int p = 0; p < M1_State.NVARS; p++) {
             for (int q = 0; q < M1_State.NVARS; q++) {
                 M1_State.index_Lag_i = p;
                 M1_State.index_Lag_j = q;
                 
                 // Compute H^{(2)}_{imp}
                 switch (M1_State.Dimension) {
                     case ONE_DIMENSIONAL:
                         One_Dimensional_Quadrature(&d2_In_lower_dlambda_2, 1, H_TWO_Matrix_Entries, &M1_State);
                         break;
                     case THREE_DIMENSIONAL:
                         Quadrature_3D_Func(&d2_In_lower_dlambda_2, 1, H_TWO_Matrix_Entries, &M1_State);  
                         break;
                     default:
                         cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
                         exit(0);
                         break;
                }
                      
                 // Compute V_{ijkl} = M_{ijkl} - H^{(2)}_{ipq} (A^{(2)}_{qkl} A^{(1)}_{pj} + A^{(2)}_{pjl} A^{(1)}_{qk} + A^{(2)}_{pjk} A^{(1)}_{ql})
                 // where V_{ijkl} = H^{(1)}_{ip} A^{(3)}_{pjkl}
                 
                 index_temp_1 = (q*M1_State.NVARS + k)*M1_State.NVARS + l;
                 index_temp_2 = (p*M1_State.NVARS + j)*M1_State.NVARS + l;
                 index_temp_3 = (p*M1_State.NVARS + j)*M1_State.NVARS + k;
                 
                 temp_mat_V_ijkl[i] += -d2_In_lower_dlambda_2*(A_TWO[index_temp_1]*A_ONE[p*M1_State.NVARS+j] + A_TWO[index_temp_2]*A_ONE[q*M1_State.NVARS+k] + A_TWO[index_temp_3]*A_ONE[q*M1_State.NVARS+l]);
             }
         }
     }
     
     for (int i = 0; i < M1_State.NVARS; i++) {
         A_THREE[i] = 0.0;
         // Now compute A^{(3)}_{ijkl} = inv(H^{(1)})_{ip} V_{pjkl}
         for (int p = 0; p < M1_State.NVARS; p++) {
             A_THREE[i] += A_ONE[i*M1_State.NVARS+p]*temp_mat_V_ijkl[p];
         }
     }
     
     delete[] temp_mat_V_ijkl;
 }

// ******************************************************************************
// This routine computes the matrix of second-derivatives of the Eddington factor
// or the partial angular moments with respect to the lower-order angular moments
// ===> {\partial^{2} I^{(2)}} {\partial I^{(j)} \partial I^{(k)}}
// ******************************************************************************
 long double Calculate_Angular_Moments_Second_Derivatives(long double *d_Chi2_In_plus_dlambda, long double *A_ONE, const int &j, const int &k, M1_State_Param &M1_State) {
     long double d2_I2_In_plus_d_lambda_2;
     int index_temp;
     long double *A_TWO;
     long double d2_Chi2_In_plus_dIn2 = 0.0;
     A_TWO = new long double[M1_State.NVARS];
     
     Compute_A_TWO_Matrix(A_TWO, A_ONE, j, k, M1_State);
     
     // A^{(1)}_{ij} ==> {\partial \lam_i} {\partial I^{(j)}}
     // A^{(2)}_{ijk} ==> {\partial^{2} \lam_i} {\partial I^{(j)} \partial I^{(k)}}
     
     for (int p = 0; p < M1_State.NVARS; p++) {
         for (int q = 0; q < M1_State.NVARS; q++) {
             M1_State.index_Lag_i = p;
             M1_State.index_Lag_j = q;
             
             // Compute {\partial^{2} I^{(2)}} {\partial \lam_{q} \partial \lam_{p}}
             switch (M1_State.Dimension) {
                 case ONE_DIMENSIONAL:
                     One_Dimensional_Quadrature(&d2_I2_In_plus_d_lambda_2, 1, Angular_Moments_Second_Derivatives, &M1_State);
                     break;
                 case THREE_DIMENSIONAL:
                     Partial_Quadrature_3D_Func(&d2_I2_In_plus_d_lambda_2, 1, Angular_Moments_Second_Derivatives, &M1_State);
                     break;
                 default:
                     cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
                     exit(0);
                     break;
            }
             
             // Now compute {\partial^{2} I^{(2)}} {\partial I^{(k)} \partial I^{(j)}} = 
             // A^{(1)}_{pj} A^{(1)}_{qk} {\partial^{2} I^{(2)}} {\partial \lam_{q} \partial \lam_{p}}
             d2_Chi2_In_plus_dIn2 += d2_I2_In_plus_d_lambda_2*A_ONE[p*M1_State.NVARS+j]*A_ONE[q*M1_State.NVARS+k];
        }
    }
    
    for (int p = 0; p < M1_State.NVARS; p++) {
        // Now compute {\partial^{2} I^{(2)}} {\partial I^{(k)} \partial I^{(j)}} = 
        // A^{(2)}_{pjk} {\partial I^{(2)}} {\partial \lam_{p}}
        d2_Chi2_In_plus_dIn2 += d_Chi2_In_plus_dlambda[p]*A_TWO[p];
    }
    
    // So far we have compute d2_I2_In_plus_dIn2, from which we can deduce d2_Chi2_In_plus_dIn2 as follows
    d2_Chi2_In_plus_dIn2 *= M1_State.I0;
    
    delete[] A_TWO;
    
    return d2_Chi2_In_plus_dIn2;
 }
 
// *********************************************************************************
// This routine computes the matrix of third-derivatives of the Eddington factor 
// or the partial angular moments with respect to the lower-order angular moments
// ===> {\partial^{3} I^{(2)}} {\partial I^{(j)} \partial I^{(k)} \partial I^{(l)}}
// *********************************************************************************
 long double Calculate_Angular_Moments_Third_Derivatives(long double *dChi2dlambda, long double *A_ONE, const int &i, const int &j, const int &k, M1_State_Param &M1_State) {
     long double d2_Chi2_In_plus_d_lambda_2, d3_Chi2_In_plus_d_lambda_3;
     long double *A_TWO, *A_THREE;
     long double d3N3dlIn = 0.0;
     int index_temp;
     A_TWO = new long double[M1_State.NVARS];
     A_THREE = new long double[M1_State.NVARS];
     
     cout << "Seg fault occuring here sometimes !!!!!!!!!!!!!!!!!!!!!!" << endl;
     
     Compute_A_TWO_Matrix(A_TWO, A_ONE, j, k, M1_State);
     
//      Compute_A_THREE_Matrix(A_THREE, A_ONE, j, k, M1_State);
     
     for (int p = 0; p < M1_State.NVARS; p++) {
         for (int q = 0; q < M1_State.NVARS; q++) {
             for (int r = 0; r < M1_State.NVARS; r++) {
                 M1_State.index_Lag_i = p;
                 M1_State.index_Lag_j = q;
                 M1_State.index_Lag_k = r;
                 
                 switch (M1_State.Dimension) {
                     case ONE_DIMENSIONAL:
                         One_Dimensional_Quadrature(&d3_Chi2_In_plus_d_lambda_3, 1, Angular_Moments_Third_Derivatives, &M1_State);
                         break;
                     case THREE_DIMENSIONAL:
                         Partial_Quadrature_3D_Func(&d3_Chi2_In_plus_d_lambda_3, 1, Angular_Moments_Third_Derivatives, &M1_State);
                         break;
                     default:
                         cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
                         exit(0);
                         break;
                }
                 
                 d3N3dlIn += d3_Chi2_In_plus_d_lambda_3*A_ONE[p*M1_State.NVARS+i]*A_ONE[q*M1_State.NVARS+j]*A_ONE[r*M1_State.NVARS+k];
             }
        }
    }
    
    for (int p = 0; p < M1_State.NVARS; p++) {
         for (int q = 0; q < M1_State.NVARS; q++) {
             M1_State.index_Lag_i = p;
             M1_State.index_Lag_j = q;
             
             // Compute {\partial^{2} I^{(3)}} {\partial \lam_{q} \partial \lam_{p}}
             switch (M1_State.Dimension) {
                 case ONE_DIMENSIONAL:
                     One_Dimensional_Quadrature(&d2_Chi2_In_plus_d_lambda_2, 1, Angular_Moments_Second_Derivatives, &M1_State);
                     break;
                 case THREE_DIMENSIONAL:
                     Partial_Quadrature_3D_Func(&d2_Chi2_In_plus_d_lambda_2, 1, Angular_Moments_Second_Derivatives, &M1_State);
                     break;
                 default:
                     cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
                     exit(0);
                     break;
            }
             
             // Now compute {\partial^{2} I^{(3)}} {\partial I^{(k)} \partial I^{(j)}} = 
             // A^{(1)}_{pj} A^{(1)}_{qk} {\partial^{2} I^{(3)}} {\partial \lam_{q} \partial \lam_{p}}
             d3N3dlIn += d2_Chi2_In_plus_d_lambda_2*A_ONE[p*M1_State.NVARS+j]*A_ONE[q*M1_State.NVARS+k];
        }
    }
    
    for (int p = 0; p < M1_State.NVARS; p++) {
        index_temp = ((p*M1_State.NVARS + i)*M1_State.NVARS + j)*M1_State.NVARS + k;
        d3N3dlIn += dChi2dlambda[p]*A_THREE[index_temp];
    }
    
    d3N3dlIn *= pow(M1_State.I0, 2);
    
    delete[] A_TWO;
    delete[] A_THREE;
    
    return d3N3dlIn;
 }
 
 //  void Calculate_Hessian_Matrix(long double *Q_data, const long double *x, const long double *Sk, M1_State_Param &M1_State) {
//      long double Q_val_ij;
//      M1_State.set_x(x);
//      M1_State.set_Sk(Sk);
//      M1_State.set_Sk_cur(Sk);
//      M1_State.set_ak(Sk);
//      
//      for (int i = 0; i < M1_State.NVARS; i++) {
//          for (int j = 0; j < M1_State.NVARS; j++) {
//              M1_State.index_p = i;
//              M1_State.index_a = j;
//              
//              switch (M1_State.Dimension) {
//                  case ONE_DIMENSIONAL:
//                      One_Dimensional_Quadrature(&Q_val_ij, 1, H_ONE_Matrix_Entries, &M1_State);
//                      break;
//                  case THREE_DIMENSIONAL:
//                      Quadrature_3D_Func(&Q_val_ij, 1, H_ONE_Matrix_Entries, &M1_State);
//                      break;
//                  default:
//                      cout << " ********************* Dimension Not Specified !!!!! ********************* " << endl;
//                      exit(0);
//                      break;
//             }
//              Q_data[i*M1_State.NVARS + j] = Q_val_ij;
//          }
//      }
//  }
 
//  void Calculate_Inverse_Hessian_Matrix(long double *Q_data, const long double *x, const long double *Sk, M1_State_Param &M1_State) {
//      long double *Q, *R;
//      Q = new long double[M1_State.NVARS*M1_State.NVARS];
//      R = new long double[M1_State.NVARS*M1_State.NVARS];
//      
//      Calculate_Hessian_Matrix(Q_data, x, Sk, M1_State);
//      
//      
// //      cout << "************************** Original !!!!!!!! ***********************" << endl;
// // 
// //      for (int i = 0; i < M1_State.NVARS; i++) {
// //          for (int j = 0; j < M1_State.NVARS; j++) {
// //              cout << "i = " << i << "   " << "j = " << j << "   " << "Q_data = " << Q_data[i*M1_State.NVARS+j] << endl;
// //         }
// //     }
//      
//      Modified_Gram_Schmidt_Factorization(M1_State.NVARS, Q_data, Q, R);
//      
// //      for (int i = 0; i < M1_State.NVARS; i++) {
// //          for (int j = 0; j < M1_State.NVARS; j++) {
// //              cout << "i = " << i << "   " << "j = " << j << "   " << "Q = " << Q[i*M1_State.NVARS+j] << "   " << "R = " << R[i*M1_State.NVARS+j] << endl;
// //         }
// //     }
//      
//      for (int i = 0; i < M1_State.NVARS; i++) {
//          for (int j = 0; j < M1_State.NVARS; j++) {
//              Q_data[i*M1_State.NVARS+j] = Q[j*M1_State.NVARS+i];
//         }
//     }
//      
//      /* Backward substitution for discovering values of unknowns */
//      for(int k = 0; k < M1_State.NVARS; k++) {                   
//          for(int i = M1_State.NVARS-1; i >= 0; i--) {                     
//              Q[i*M1_State.NVARS + k] = Q_data[i*M1_State.NVARS + k];
//              for(int j = i+1; j < M1_State.NVARS;j++) {
//                  if(i != j) {
//                      Q[i*M1_State.NVARS + k] = Q[i*M1_State.NVARS + k] - R[i*M1_State.NVARS + j]*Q[j*M1_State.NVARS + k];
//                 }          
//             }
//             Q[i*M1_State.NVARS + k] = Q[i*M1_State.NVARS + k]/R[i*M1_State.NVARS + i];
//         }
//      }
//      
// //      cout << "************************** Inverse !!!!!!!! ***********************" << endl;
// 
//      for (int i = 0; i < M1_State.NVARS; i++) {
//          for (int j = 0; j < M1_State.NVARS; j++) {
//              Q_data[i*M1_State.NVARS+j] = Q[i*M1_State.NVARS+j];
// //              cout << "i = " << i << "   " << "j = " << j << "   " << "Q_data = " << Q_data[i*M1_State.NVARS+j] << endl;
//         }
//     }
//      
//      delete[] Q;
//      delete[] R;
//  }
 
//  long double myconstraint(unsigned n, const long double *x, const long double *Sk, long double *grad, void *data) {
//     long double temp_val = 0.0;
//     long double constr_val = 0.0;
//     my_constraint_data *d = (my_constraint_data *) data;
//     
//     if (grad) {
//         for (int i = 0; i < n; i++) {
//             grad[i] = 0.0;
//         }
//     }
//             
//     switch (M1_State_Param::Regime) {
//         case GRAY:
//         case BOSE_EINSTEIN:
//         case LOGARITHMIC_LIMIT:
//             
//             if (grad) {
//                 for (int i = 0; i < n; i++) {
//                     for (int j = 0; j < n; j++) {
//                         grad[i] += Sk[i*n+j]*(*d)[j];
//                     }
//                 }
//             }
//             
//             constr_val = 0.0;
//             for (int i = 0; i < n; i++) {
//                 for (int j = 0; j < n; j++) {
//                     constr_val += Sk[i*n+j]*(*d)[j]*x[i];
//                 }
//             }
//             
//             temp_val = constr_val + 1.0e-32;
//             break;
//         case HYPERBOLIC_LIMIT:
//             // Positivity of the distribution does not impose any constraint
//             // on the Lagrange multipliers for the Boltzmann distribution
//             if (grad) {
//                 for (int i = 0; i < n; i++) {
//                     grad[i] = 0.0;
//                 }
//             }
//             temp_val = 0.0;
//             break;
//     }
//     
//     // temp_val = 1.0e2*temp_val;
//             
//     return temp_val;
//  }


long double myconstraint(unsigned n, const long double *x, const long double *Sk, long double *grad, void *data_constraints) {
    long double temp_val = 0.0;
    long double constr_val = 0.0, constr_val_max = -1.0e12;
    my_constraint_data data_max, data;
    M1_State_Param *M1_State_temp = (M1_State_Param *) data_constraints;
    
    int Id_Angle;
    long double *poly;
    poly = new long double[M1_State_temp->NVARS];
    
    long double phi_start, phi_end;
    long double mu_start, mu_end;
    switch (M1_State_Param::Regime) {
        case GRAY:
        case BOSE_EINSTEIN:
        case LOGARITHMIC_LIMIT:
            for (int i_mu_quad = 0; i_mu_quad < M1_State_temp->N_intervals_mu_quad; i_mu_quad++) {
                mu_start = M1_State_temp->bounds_quad_mu_refin[i_mu_quad];
                mu_end = M1_State_temp->bounds_quad_mu_refin[i_mu_quad+1];
                
                M1_State_temp->compute_x_quad(mu_start, mu_end);
                
                for (int Id_Angle = 0; Id_Angle < M1_State_temp->order_quad; Id_Angle++) {
                    generate_polynomials_1D(poly, Id_Angle, *M1_State_temp);
                    data.a = poly[0];
                    data.b = poly[1];  
                    
                    constr_val = 0.0;
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                            constr_val += Sk[i*n+j]*data[j]*x[i];
                        }
                    }
                    
                    constr_val_max = max(constr_val_max, constr_val);
                    if (constr_val_max == constr_val) {
                        data_max = data;
                    }
                }
            }
            
            if (grad) {
                // Since the location of constr_val_max may change, the gradient should be set
                // to zero to avoid convergence issues??
                for (int i = 0; i < n; i++) {
                    grad[i] = 0.0;
//                     for (int j = 0; j < n; j++) {
//                         grad[i] += Sk[i*n+j]*data_max[j];
//                     }
//                     // Check what is wrong with SLSQP sign convention for gradients of inequality constraints ??
//                     grad[i] = -grad[i];
//                     grad[i] *= 1.0e6;
                }
            }
            
            temp_val = 1.0e6*(constr_val_max + 1.0e-16);
            break;
        case HYPERBOLIC_LIMIT:
            // Positivity of the distribution does not impose any constraint
            // on the Lagrange multipliers for the Boltzmann distribution
            if (grad) {
                for (int i = 0; i < n; i++) {
                    grad[i] = 0.0;
                }
            }
            temp_val = 0.0;
            break;
    }
    
    delete[] poly;
            
    return temp_val;
 }

void One_Dimensional_Quadrature(long double *Vals, const int &NFUN, max_ent_obj_grad_hess func, void *fdata) {
    int order_quad;
    long double Temp_vals[NFUN];
    M1_State_Param *M1_State = (M1_State_Param *) fdata;
    order_quad = M1_State->order_quad;
    
    for (int i = 0; i < NFUN; i++) {
        Vals[i] = 0.0;
    }
    
    long double mu_start, mu_end;
    for (int i_refin = 0; i_refin < M1_State->N_intervals_mu_quad; i_refin++) {
        mu_start = M1_State->bounds_quad_mu_refin[i_refin];
        mu_end = M1_State->bounds_quad_mu_refin[i_refin+1];
        
        // cout << "mu_start = " << mu_start << "  " << "mu_end = " << mu_end << endl;
        
        M1_State->compute_x_quad(mu_start, mu_end);
        
        for ( int m = 0; m < order_quad; m++ ) {
            func(Temp_vals, NFUN, m, fdata);
            for (int i = 0; i < NFUN; i++) {
                Vals[i] = Vals[i] + M1_State->w_quad_1D[m] * Temp_vals[i];
            }
        }
    }
}

void One_Dimensional_Quadrature_Partial_Moments(long double *Vals, const int &NFUN, max_ent_obj_grad_hess func, void *fdata) {
    int order_quad;
    long double Temp_vals[NFUN];
    M1_State_Param *M1_State = (M1_State_Param *) fdata;
    order_quad = M1_State->order_quad;
    
    for (int i = 0; i < NFUN; i++) {
        Vals[i] = 0.0;
    }
    
    long double mu_start, mu_end;
    for (int i_refin = 0; i_refin < M1_State->N_intervals_mu_quad; i_refin++) {
        mu_start = M1_State->bounds_quad_mu_refin[i_refin];
        mu_end = M1_State->bounds_quad_mu_refin[i_refin+1];
        
        M1_State->compute_x_quad(mu_start, mu_end);
        
        for ( int m = 0; m < order_quad; m++ ) {
            if (M1_State->x_quad[m] >= 0.0) {
                func(Temp_vals, NFUN, m, fdata);
                for (int i = 0; i < NFUN; i++) {
                    Vals[i] = Vals[i] + M1_State->w_quad_1D[m] * Temp_vals[i];
                }
            }
        }
    }
}
 
 long double Inverse_Mobius_Transformation(const long double &ratio, const long double &Length_Scale) {
     long double Val;
     Val = -Length_Scale*log((1.0 - ratio)/2.0);
//      Val = Length_Scale*(ratio + 1.0)/(1.0 - ratio);
     
//      if (Val > I0_val_max) {
//          Val = I0_val_max;
//      }
     return Val;
 }
 
 long double Mobius_Transformation(const long double &Moment, const long double &Length_Scale) {
     long double ratio;
     ratio = 1.0 - 2.0*exp(-Moment/Length_Scale);
     
//      if (Moment > I0_val_max) {
//          ratio = 1.0 - 2.0*exp(-I0_val_max/Length_Scale);
//      }
//      ratio = (Moment - Length_Scale)/(Moment + Length_Scale);
     return ratio;
 }
 
 long double d_ratio_E_d_I0(long double &Moment, const long double &Length_Scale) {
     long double ratio;
     ratio = (2.0/Length_Scale)*exp(-Moment/Length_Scale);
//      ratio = 2.0*Length_Scale/pow(Moment + Length_Scale, 2);
     
//      if (Moment > I0_val_max) {
//          ratio = 0.0;
//      }
     
     return ratio;
 }
 
 long double d_ratio_E_d_L_Chi2(long double &Moment, const long double &Length_Scale) {
     long double ratio;
     ratio = -2.0*(Moment/(Length_Scale*Length_Scale))*exp(-Moment/Length_Scale);
//      ratio = -2.0*Moment/pow(Moment + Length_Scale, 2);
     return ratio;
 }

long double Rotation_Matrix_X_axis(const int &i, const int &j, const long double &cos_angle, const long double &sin_angle) {
    long double mat = 0.0;
    switch (i) {
        case 0:
            switch (j) {
                case 0:
                    mat = 1.0;
                    break;
            };
            break;
        case 1:
            switch (j) {
                case 1:
                    mat = cos_angle;
                    break;
                case 2:
                    mat = -sin_angle;
                    break;
            };
            break;
        case 2:
            switch (j) {
                case 1:
                    mat = sin_angle;
                    break;
                case 2:
                    mat = cos_angle;
                    break;
            };
            break;
    };
    
    return mat;
}

long double Rotation_Matrix_Y_axis(const int &i, const int &j, const long double &cos_angle, const long double &sin_angle) {
    long double mat = 0.0;
    switch (i) {
        case 0:
            switch (j) {
                case 0:
                    mat = cos_angle;
                    break;
                case 2:
                    mat = sin_angle;
                    break;
            };
            break;
        case 1:
            switch (j) {
                case 1:
                    mat = 1.0;
                    break;
            };
            break;
        case 2:
            switch (j) {
                case 0:
                    mat = -sin_angle;
                    break;
                case 2:
                    mat = cos_angle;
                    break;
            };
            break;
    };
    
    return mat;
} 

long double Rotation_Matrix_Z_axis(const int &i, const int &j, const long double &cos_angle, const long double &sin_angle) {
    long double mat = 0.0;
    switch (i) {
        case 0:
            switch (j) {
                case 0:
                    mat = cos_angle;
                    break;
                case 1:
                    mat = -sin_angle;
                    break;
            };
            break;
        case 1:
            switch (j) {
                case 0:
                    mat = sin_angle;
                    break;
                case 1:
                    mat = cos_angle;
                    break;
            };
            break;
        case 2:
            switch (j) {
                case 2:
                    mat = 1.0;
                    break;
            };
            break;
    };
    
    return mat;
}

void Rotate_Frame(long double *x, const int &NVARS, const long double &x_Lebed, const long double &y_Lebed, const long double &z_Lebed) {
//     long double *x_rot_interm;
//     long double theta, phi;
//     long double cos_angle, sin_angle;
//     
//     x_rot_interm = new long double[NVARS];
//     
//     //     The Lebedev quadrature is defined such that 
//     //     xi = cos ( thetai ) * sin ( phii )
//     //     yi = sin ( thetai ) * sin ( phii )
//     //     zi = cos ( phii )
//     
//     xyz_to_tp_mine(x_Lebed, y_Lebed, z_Lebed, &theta, &phi);
//     
//     cos_angle = cos(theta);
//     sin_angle = sin(theta);
//     
//     // Rotating Vector
//     for (int i = 1; i < NVARS; i++) {
//         x_rot_interm[i] = 0.0;
//         for (int j = 1; j < NVARS; j++) {
//             x_rot_interm[i] += x[j]*Rotation_Matrix_Z_axis(i-1, j-1, cos_angle, sin_angle);
//         }
//     }
//     
//     cos_angle = cos(phi - PI/2.0);
//     sin_angle = sin(phi - PI/2.0);
//     
//     // Rotating Vector
//     for (int i = 1; i < NVARS; i++) {
//         x[i] = 0.0;
//         for (int j = 1; j < NVARS; j++) {
//             x[i] += x_rot_interm[j]*Rotation_Matrix_Y_axis(i-1, j-1, cos_angle, sin_angle);
//         }
//     }
//     
//     delete[] x_rot_interm;
}

void Rotate_Frame_Vector(long double *x, const int &NVARS, const long double &x_Lebed, const long double &y_Lebed, const long double &z_Lebed) {
    long double *x_rot_interm;
    long double theta, phi;
    long double cos_angle, sin_angle;
    
    x_rot_interm = new long double[NVARS];
    
    //     The Lebedev quadrature is defined such that 
    //     xi = cos ( thetai )
    //     yi = sin ( thetai ) * cos ( phii )
    //     zi = sin ( thetai ) * sin ( phii )
    
    xyz_to_tp_mine(x_Lebed, y_Lebed, z_Lebed, &theta, &phi);
    
    cos_angle = cos(phi);
    sin_angle = sin(phi);
    
    // Rotating Vector
    for (int i = 0; i < NVARS; i++) {
        x_rot_interm[i] = 0.0;
        for (int j = 0; j < NVARS; j++) {
            x_rot_interm[i] += x[j]*Rotation_Matrix_X_axis(i, j, cos_angle, sin_angle);
        }
    }
    
    cos_angle = cos(theta);
    sin_angle = sin(theta);
    
    // Rotating Vector
    for (int i = 0; i < NVARS; i++) {
        x[i] = 0.0;
        for (int j = 0; j < NVARS; j++) {
            x[i] += x_rot_interm[j]*Rotation_Matrix_Z_axis(i, j, cos_angle, sin_angle);
        }
    }
    
    delete[] x_rot_interm;
}

void Rotate_Frame_Vector_back(long double *x, const int &NVARS, const long double &x_Lebed, const long double &y_Lebed, const long double &z_Lebed) {
    long double *x_rot_interm;
    long double theta, phi;
    long double cos_angle, sin_angle;
    
    x_rot_interm = new long double[NVARS];
    
    
    for (int i = 0; i < NVARS; i++) {
        x_rot_interm[i] = x[i];
    }
    
    //     The Lebedev quadrature is defined such that 
    //     xi = cos ( thetai )
    //     yi = sin ( thetai ) * cos ( phii )
    //     zi = sin ( thetai ) * sin ( phii )
    
    xyz_to_tp_mine(x_Lebed, y_Lebed, z_Lebed, &theta, &phi);
    
//     cout << "theta = " << theta << "  " << "phi = " << phi << endl; 
    
    theta = -theta;
    phi = -phi;
    
    cos_angle = cos(theta);
    sin_angle = sin(theta);
    
    // Rotating Vector
    for (int i = 0; i < NVARS; i++) {
        x[i] = 0.0;
        for (int j = 0; j < NVARS; j++) {
            x[i] += x_rot_interm[j]*Rotation_Matrix_Z_axis(i, j, cos_angle, sin_angle);
        }
    }
    
    cos_angle = cos(phi);
    sin_angle = sin(phi);
    
    // Rotating Vector
    for (int i = 0; i < NVARS; i++) {
        x_rot_interm[i] = 0.0;
        for (int j = 0; j < NVARS; j++) {
            x_rot_interm[i] += x[j]*Rotation_Matrix_X_axis(i, j, cos_angle, sin_angle);
        }
    }
    
    for (int i = 0; i < NVARS; i++) {
        x[i] = x_rot_interm[i];
    }
    
    delete[] x_rot_interm;
}

long double Heaviside(const long double &x) {
    long double val;
    long double TOLERANCE = 1.0e-6;
    if (x > TOLERANCE) {
        val = 1.0;
    } else if (x < -TOLERANCE) {
        val = 0.0;
    } else {
        val = 0.5;
    }
    
    return val;
}

void Modified_Gram_Schmidt_Factorization(const int &n, const long double *A, long double *Q, long double *R) {
    // Modified Gram Schmidt algorithm wih reorthogonalization
    long double Q_val_ak_pm, Q_val_pm2, *Sk_prev;
    long double norm_a;
    long double proj_Q_A, proj_Q_Q;
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Q[i*n+j] = 0.0;
        }
    }

    for (int k = 0; k < n; k++){ // k represents the number of polynomials in the non-orthogonal basis
        for (int i = 0; i < n; i++) {
            Q[i*n + k] = A[i*n + k];
        }  // end for loop i
        for (int l = 0; l < 2; l++) { // loop for reorthogonalization
            for(int m = 0; m < k; m++) {
                // evaluate norm`
                proj_Q_A = 0.0;
                proj_Q_Q = 0.0;
                for (int i = 0; i < n; i++) {
                    proj_Q_A += Q[i*n + k]*Q[i*n + m];
                    proj_Q_Q += Q[i*n + m]*Q[i*n + m];
                }  // end for loop i
                for (int i = 0; i < n; i++) {
                    // orthogonalize k^{th} column of A with respect to columns up to (k-1) of Q
                    Q[i*n + k] = Q[i*n + k] - (proj_Q_A/proj_Q_Q)*Q[i*n + m];
                }  // end for loop i
            } // end for loop m
        } // end for loop l
        
        proj_Q_Q = 0.0;
        for (int i = 0; i < n; i++) {
            proj_Q_Q += Q[i*n + k]*Q[i*n + k];
        }  // end for loop i
        for (int i = 0; i < n; i++) {
            Q[i*n+k] = Q[i*n+k]/sqrt(proj_Q_Q); // Orthonormalize
//             cout << "i = " << i << "   " << "k = " << k << "   " << "Q = " << Q[i*n+k] << "   " << "proj_Q_Q = " << proj_Q_Q << endl;
        }
    } // end for loop k
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            R[i*n+j] = 0.0;
            for (int k = 0; k < n; k++) {
                R[i*n+j] += Q[k*n+i]*A[k*n+j];
            }
//             cout << "i = " << i << "   " << "j = " << j << "   " << "R = " << R[i*n+j] << endl;
        }
    }
}

void setup_chebyshev_nodes_data(const int &N_points_cheby, long double *x_cheby) {
    for (int index = 0; index < N_points_cheby; index++) {
        x_cheby[index] = zeros_shifted( index + N_points_cheby - 1, 2*(N_points_cheby - 1) + 1, -1.0, 1.0, CHEBYSHEV_SECOND_KIND_DISTRIBUTION);
    }
}

long double Cartesian_to_spherical_Coordinates_Conversion_Jacobian(const int &i, const int &j, const long double &r, const long double &cos_theta, const long double &sin_theta, const long double &cos_phi, const long double &sin_phi) {
    long double mat = 0.0;
    switch (i) {
        case 0:
            switch (j) {
                case 0:
                    mat = cos_theta;
                    break;
                case 1:
                    mat = -r*sin_theta;
                    break;
            };
            break;
        case 1:
            switch (j) {
                case 0:
                    mat = sin_theta*cos_phi;
                    break;
                case 1:
                    mat = r*cos_theta*cos_phi;
                    break;
                case 2:
                    mat = -r*sin_theta*sin_phi;
                    break;
            };
            break;
        case 2:
            switch (j) {
                case 0:
                    mat = sin_theta*sin_phi;
                    break;
                case 1:
                    mat = r*cos_theta*sin_phi;
                    break;
                case 2:
                    mat = r*sin_theta*cos_phi;
                    break;
            };
            break;
    };
    
    return mat;
}

long double Cartesian_to_spherical_Coordinates_Conversion_Jacobian_Second_Derivatives(const int &i, const int &j, const long double &r, const long double &cos_theta, const long double &sin_theta, const long double &cos_phi, const long double &sin_phi) {
    long double mat = 0.0;
    switch (i) {
        case 0:
            switch (j) {
                case 0:
                    mat = cos_theta;
                    break;
                case 1:
                    mat = -r*cos_theta;
                    break;
            };
            break;
        case 1:
            switch (j) {
                case 0:
                    mat = sin_theta*cos_phi;
                    break;
                case 1:
                    mat = -r*sin_theta*cos_phi;
                    break;
                case 2:
                    mat = -r*sin_theta*sin_phi;
                    break;
            };
            break;
        case 2:
            switch (j) {
                case 0:
                    mat = sin_theta*sin_phi;
                    break;
                case 1:
                    mat = -r*sin_theta*sin_phi;
                    break;
                case 2:
                    mat = r*sin_theta*cos_phi;
                    break;
            };
            break;
    };
    
    return mat;
}

long double Cartesian_to_spherical_Coordinates_Jacobian(const long double &norm_f, const long double &x_val, const long double &y_val, const long double &z_val, const long double &dIn_Plus_dN1_1, const long double &dIn_Plus_dN1_2, const long double &dIn_Plus_dN1_3, const int &VAR_NUM) {
    long double cos_theta, sin_theta, cos_phi, sin_phi;
    long double dIn_Plus_dmu;
    long double dN1_1_dtheta, dN1_2_dtheta, dN1_3_dtheta;
    
    cos_theta = x_val;
    sin_theta = sqrt(1.0 - pow(cos_theta,2));
    
    cos_phi = y_val/sin_theta;
    sin_phi = z_val/sin_theta;
    
    if (sin_theta < 1.0e-8) {
        cos_phi = 1.0;
        sin_phi = 0.0;
    }
    
    dN1_1_dtheta = Cartesian_to_spherical_Coordinates_Conversion_Jacobian(0, VAR_NUM, norm_f, cos_theta, sin_theta, cos_phi, sin_phi);
    dN1_2_dtheta = Cartesian_to_spherical_Coordinates_Conversion_Jacobian(1, VAR_NUM, norm_f, cos_theta, sin_theta, cos_phi, sin_phi);
    dN1_3_dtheta = Cartesian_to_spherical_Coordinates_Conversion_Jacobian(2, VAR_NUM, norm_f, cos_theta, sin_theta, cos_phi, sin_phi);
    
    dIn_Plus_dmu = dIn_Plus_dN1_1*dN1_1_dtheta + dIn_Plus_dN1_2*dN1_2_dtheta + dIn_Plus_dN1_3*dN1_3_dtheta;
    
    if (VAR_NUM == 1) {
        dIn_Plus_dmu /= -sin_theta;
    }
    
    return dIn_Plus_dmu;
}

long double Cartesian_to_spherical_Coordinates_Jacobian_Second_Derivatives(const long double &norm_f, const long double &x_val, const long double &y_val, const long double &z_val, const long double &dIn_Plus_dN1_1, const long double &dIn_Plus_dN1_2, const long double &dIn_Plus_dN1_3, const int &VAR_NUM) {
    long double cos_theta, sin_theta, cos_phi, sin_phi;
    long double dIn_Plus_dmu;
    long double dN1_1_dtheta, dN1_2_dtheta, dN1_3_dtheta;
    
    cos_theta = x_val;
    sin_theta = sqrt(1.0 - pow(cos_theta,2));
    
    cos_phi = y_val/sin_theta;
    sin_phi = z_val/sin_theta;
    
    if (sin_theta < 1.0e-8) {
        cos_phi = 1.0;
        sin_phi = 0.0;
    }
    
    dN1_1_dtheta = Cartesian_to_spherical_Coordinates_Conversion_Jacobian_Second_Derivatives(0, VAR_NUM, norm_f, cos_theta, sin_theta, cos_phi, sin_phi);
    dN1_2_dtheta = Cartesian_to_spherical_Coordinates_Conversion_Jacobian_Second_Derivatives(1, VAR_NUM, norm_f, cos_theta, sin_theta, cos_phi, sin_phi);
    dN1_3_dtheta = Cartesian_to_spherical_Coordinates_Conversion_Jacobian_Second_Derivatives(2, VAR_NUM, norm_f, cos_theta, sin_theta, cos_phi, sin_phi);
    
    dIn_Plus_dmu = dIn_Plus_dN1_1*dN1_1_dtheta + dIn_Plus_dN1_2*dN1_2_dtheta + dIn_Plus_dN1_3*dN1_3_dtheta;
    
    if (VAR_NUM == 1) {
        dIn_Plus_dmu /= -sin_theta;
    }
    
    return dIn_Plus_dmu;
}

void xyz_to_tp_mine ( long double x, long double y, long double z, long double *t, long double *p )

//****************************************************************************80
//
//  Purpose:
//
//    XYZ_TO_TP converts (X,Y,Z) to (Theta,Phi) coordinates on the unit sphere.
//
//  Modified:
//
//    09 September 2010
//
//  Author:
//
//    Dmitri Laikin
//
//  Parameters:
//
//    Input, long double X, Y, Z, the Cartesian coordinates of a point
//    on the unit sphere.
//
//    Output, long double T, P, the Theta and Phi coordinates of
//    the point.
//
{
  long double ang_x;
  long double fact;
  long double pi = 3.14159265358979323846;

  *t = acos ( x );

  fact = sqrt ( y * y + z * z );

  if ( 0 < fact ) 
  {
    ang_x = acos ( y / fact );
  }
  else
  {
    ang_x = acos ( y );
  }

  if ( z < 0 ) 
  {
    ang_x = - ang_x;
  }
  *p = ang_x;

  return;
}

long double roundval( const long double &val ) {
    if ( val < 0.0 ) {
        return -floor(fabs(val) + 0.5);
    }
    return floor(val + 0.5);
}

long double roundn( const long double &val, const int &n ) {
    long double value, temp_val;
    long double prec;
    prec = pow(10.0, n);
    temp_val = val * prec;
    
    value = roundval( temp_val ) / prec;
    
    return value;
}
