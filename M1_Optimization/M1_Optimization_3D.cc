#ifndef _M1_Model_1D_H_INCLUDED
#include "../M1_Optimization/M1_Optimization.h"
#endif // _M1_Model_1D_H_INCLUDED

const int M1_State_Param::NFUNS_PM_3D;

// ********************************************************************************
// This routines call the optimization algorithm used for solving the maximum-
// entropy problem for any given set of angular moments up to first-order, for 
// gray or non-gray radiation in three dimensions
// ********************************************************************************
void NLOPT_Optimization_Algorithm_NON_GRAY_M1_3D(record_Chi2 *rec_Chi2_local, record_Partial_Moments *rec_Partial_Moments, const int &index_count, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, const int &id_proc, M1_State_Param &M1_State, Mobius_Scale_Parameters &Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1, long double *Sk_final, const long double *x_SH, const long double *y_SH, const long double *z_SH) {
    long double x[M1_State.NVARS];
    long double tol_x[M1_State.NVARS];
    long double tol_grad = 1.0e-8; //1.0e-10;
    int max_iters_nlopt;
    
    max_iters_nlopt = 2000;
    
    long double minf, min_grad; /* the minimum objective value, upon return */
    nlopt_opt opt = NULL;
    
    Set_Moments_3D(&M1_State, Var_index, num_points, x_SH, y_SH, z_SH, &Mobius_Scale_Params, struc_map_N1);
    
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
    
    for (int index_reg = 0; index_reg < 2; index_reg++) {
        if (index_reg > 0) {
            if (M1_State.id_proc == DISPLAY_ID) {
                if (M1_State.display) {
                    printf(".................Moment Regularization with r_l = %g..................\n", r_l[index_reg]);
                }
            }
        }
        
        for ( int id_refinement = 0; id_refinement <= M1_State.max_refinement_level; id_refinement++) {
            Set_Moments_3D(&M1_State, Var_index, num_points, x_SH, y_SH, z_SH, &Mobius_Scale_Params, struc_map_N1);
            
            Regularize_Moment(r_l[index_reg], &M1_State);
            
            M1_State.order_quad = 20;
            
            // Allocate and setup quadrature scheme
            M1_State.Allocate_Quad_3D();
            M1_State.compute_Omegas();
            // Setup the bounds of integrations for the current level of refinement
            M1_State.set_bounds_quad_refin(id_refinement);
            
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
            
            my_constraint_data data[M1_State.N_dirs(M1_State.order_quad, 2*M1_State.order_quad)];
            add_constraints(opt, data, M1_State);
            
            if (nlopt_optimize_orthog(opt, x, M1_State.I0, &minf, &min_grad, Sk_final) < 0) {
                if (id_proc == DISPLAY_ID) {
                    if (M1_State.display) {
                        printf("....................nlopt failed!.....................\n");
                    }
                }
                
                if (index_reg == 1 &&  id_refinement == M1_State.max_refinement_level) {
                    printf("....................nlopt failed: Exiting!.....................\n");
                    cout << "Failure with set of moments "  << "I0 = " << M1_State.I0 << "   " << "N1 = " << M1_State.N1 << "    " << "id = " << id_proc << "     " << "grad_norm = " << min_grad << endl;
                    
                    exit(0);
                    goto Exiting;
                }
            } else if (min_grad > tol_grad) {
                if (id_proc == DISPLAY_ID) {
                    if (M1_State.display) {
                        printf("***********Tolerance on gradient not satisfied g(");
                        for (int index_vars = 0; index_vars < M1_State.NVARS; index_vars++) {
                            if (index_vars < M1_State.NVARS - 1) {
                                printf("%g,", x[index_vars]);
                            } else {
                                printf("%g", x[index_vars]);
                            }
                        }
                        printf(") = %0.10g................\n", min_grad);
                    }
                }
                
                if (index_reg == 1 && id_refinement == M1_State.max_refinement_level) {
                    if (id_proc == DISPLAY_ID) {
                        if (M1_State.display) {
                            printf("....................nlopt failed: Exiting!.....................\n");
                        }
                    }
                    exit(0);
                    goto Exiting;
                }
            } else {
                if (id_proc == DISPLAY_ID) {
                    if (M1_State.display) {
                        printf("***********found minimum at f(");
                        for (int index_vars = 0; index_vars < M1_State.NVARS; index_vars++) {
                            if (index_vars < M1_State.NVARS - 1) {
                                printf("%g,", x[index_vars]);
                            } else {
                                printf("%g", x[index_vars]);
                            }
                        }
                        printf(") = %0.10g***********\n ", minf);
                        printf("***********Final gradient is min_grad = %0.10g *********** \n\n", min_grad);
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
            cout << "Eddington factor computed in 1D instead of 3D, exiting!!!" << endl;
            exit(0);
            break;
        case PARTIAL_MOMENTS:
            Compute_Partial_Moments(rec_Partial_Moments, index_count, Var_index, num_points, M1_State, x, Mobius_Scale_Params, struc_map_N1, x_SH, y_SH, z_SH, Sk_final);
            if (!M1_State.flag_Taylor_Series_Expansion) {
                // We only need to compute derivatives separately in the case where we are not
                // performing a taylor series expansion
                Compute_Partial_Moments_Derivatives(rec_Partial_Moments, index_count, Var_index, num_points, M1_State, x, Mobius_Scale_Params, struc_map_N1, x_SH, y_SH, z_SH, Sk_final);
            }
            break;
        default:
            cout << "Maximum entropy solution type not specified" << endl;
            exit(0);
            break;
    }
    
    if (opt != NULL) {
        nlopt_destroy(opt);
        opt = NULL;
    }
    
    Exiting:;
    M1_State.Deallocate_Quad_3D();
}

// ********************************************************************************
// This routines sets up the parameters required to solve the dual maximum-entropy
// problem for any given set of angular moments up to first-order, for gray or 
// non-gray radiation in one dimensional problems
// ********************************************************************************
void OPTIM_NON_GRAY_M1_3D_Array(record_Chi2 *rec_Chi2_global, record_Partial_Moments *rec_PM_global, const M1_Var_Num_Points *num_points, M1_State_Param &M1_State, const int &N_pts_Mob_Scale, const long double *Coefficients_Mobius_Scale_Fit, const Algebraic_Mapping_N1 *struc_map_N1) {
    int num_proc_per_var_E, num_proc_per_var_f;
    int id_proc_E, id_proc_f;
    int size_rec, size_rec_total;;
    int incr_E, incr_f;
    
    MPI_Datatype rec_Chi2_type, rec_PM_type;
    int lengths[1];
    const MPI_Aint displacements[1] = { 0 };
    MPI_Datatype types[1] = { MPI_LONG_DOUBLE };
    
    M1_Var_Num_Points Var_index;
    int N_points_Leg = num_points->N_pts_Leg;
    
    long double *x_SH = NULL;
    long double *y_SH = NULL;
    long double *z_SH = NULL;
    
    long double minf, min_grad; /* the minimum objective value, upon return */
    nlopt_opt opt = NULL;
    
    int N_pts_Leg_Total;
    
    long double *Sk_final;
    
    int index_count;
    
    // Setup number of variable for M1_State
    M1_State.Set_Num_Vars();
    
    // Define data structure for the length of the exponential mapping of the radiative energy
    // density in the case of non-gray radiation
    Mobius_Scale_Parameters Mobius_Scale_Params(N_pts_Mob_Scale, Coefficients_Mobius_Scale_Fit);
    
    // Create the rec_Chi2 or rec_PM datatype for MPI
    switch (M1_State.Max_Ent_Solution_Type) {
        case CLOSING_FLUX:
            lengths[0] = { 13 };
            break;
        case PARTIAL_MOMENTS:
            lengths[0] = { 75 };
            break;
        default:
            cout << "Maximum entropy solution type not specified" << endl;
            exit(0);
            break;
    } // end switch
    
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
    
    // Setup the number processors per variable
    if (M1_State.Problem_Type == NON_GRAY) {
        num_proc_per_var_E = floor(sqrt(M1_State.num_proc));
        num_proc_per_var_f = num_proc_per_var_E;
        id_proc_E = floor(M1_State.id_proc/num_proc_per_var_f);
        id_proc_f = M1_State.id_proc - id_proc_E*num_proc_per_var_f;
    } else if (M1_State.Problem_Type == GRAY) {
        num_proc_per_var_E = 1;
        num_proc_per_var_f = M1_State.num_proc;
        id_proc_E = 0;
        id_proc_f = M1_State.id_proc;   
    }
    
    if (num_points->E % num_proc_per_var_E != 0) {
        if (M1_State.id_proc == DISPLAY_ID) {
            cout << "The number num_points->E = " << num_points->E << "is not divisible by num_proc_per_var_E = " << num_proc_per_var_E << "the quotient is " << num_points->E/num_proc_per_var_E << endl;
        }
        goto Exiting;
    }
    
    if (num_points->f % num_proc_per_var_f != 0) {
        if (M1_State.id_proc == DISPLAY_ID) {
            cout << "The number num_points->f = " << num_points->f << "is not divisible by num_proc_per_var_f = " << num_proc_per_var_f << "the quotient is " << num_points->f/num_proc_per_var_f << endl;
        }
        goto Exiting;
    }
    
    record_Partial_Moments *rec_PM_local, *rec_PM_global_temp;
    // Perform allocation of data structure for Chi2 or PM for each processor
    incr_E = ceil(double(num_points->E)/double(num_proc_per_var_E));
    incr_f = ceil(double(num_points->f)/double(num_proc_per_var_f));
    size_rec = incr_E*incr_f*num_points->phi*num_points->theta;
    switch (M1_State.Max_Ent_Solution_Type) {
        case CLOSING_FLUX:
            cout << "Eddington factor must be computed in the 1D case" << endl;
            break;
        case PARTIAL_MOMENTS:
            rec_PM_local = new record_Partial_Moments[size_rec];
            break;
        default:
            cout << "Maximum entropy solution type not specified" << endl;
            exit(0);
            break;
    } // end switch
    
    // Perform allocation of data structure for Chi2 or PM on primary processor for the purpose
    // of gathering data from all processors
    size_rec_total = size_rec*M1_State.num_proc;
    
    if (M1_State.id_proc == PRIMARY_ID) {
        switch (M1_State.Max_Ent_Solution_Type) {
            case CLOSING_FLUX:
                cout << "Eddington factor must be computed in the 1D case" << endl;
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
    
    if (N_points_Leg == 1) {
        N_pts_Leg_Total = 1;
        
        x_SH = new long double[N_pts_Leg_Total];
        y_SH = new long double[N_pts_Leg_Total];
        z_SH = new long double[N_pts_Leg_Total];
        x_SH[0] = 0.0;
        y_SH[0] = 0.0;
        z_SH[0] = 0.0;
    } else if (num_points->Symmetry_Type == SYMMETRY_N1_1) {
        N_pts_Leg_Total = N_points_Leg;
        x_SH = new long double[N_pts_Leg_Total];
        y_SH = new long double[N_pts_Leg_Total];
        z_SH = new long double[N_pts_Leg_Total];
        setup_chebyshev_nodes_data(N_pts_Leg_Total, x_SH);
        for (int i = 0; i < N_pts_Leg_Total; i++) {
            y_SH[i] = sqrt(1.0 - pow(x_SH[i], 2));
            z_SH[i] = 0.0;
        }
    } else {
        N_pts_Leg_Total = 2*N_points_Leg*N_points_Leg;
        x_SH = new long double[N_pts_Leg_Total];
        y_SH = new long double[N_pts_Leg_Total];
        z_SH = new long double[N_pts_Leg_Total];
        setup_spherical_harmonics_data ( N_points_Leg, x_SH, y_SH, z_SH, NULL );
    }
    
    Sk_final = new long double[M1_State.NVARS*M1_State.NVARS];
   
    for (int i = 0; i < M1_State.NVARS; i++) {
        for (int j = 0; j < M1_State.NVARS; j++) {
            if (i == j) {
                Sk_final[i*M1_State.NVARS + j] = 1.0;
            } else {
                Sk_final[i*M1_State.NVARS + j] = 0.0;
            }
        }
    }
    
    index_count = 0;
    for (int i_E = id_proc_E*incr_E; i_E < (id_proc_E+1)*incr_E; i_E++) {
        for (int i_f = id_proc_f*incr_f; i_f < (id_proc_f+1)*incr_f; i_f++) {
            for (int i_Phi = 0 ; i_Phi < num_points->phi; i_Phi++) {
                for (int i_Theta = 0 ; i_Theta < num_points->theta; i_Theta++) {
                    Var_index.E = i_E;
                    Var_index.f = i_f;
                    Var_index.phi = i_Phi;
                    Var_index.theta = i_Theta;
                    
                    // **********************************************************************************
                    // Checks for MPI memory errors
                    // **********************************************************************************
                    if (index_count >= size_rec) {
                        cout << "index_count exceeds size_rec ====> Out of bound memory likely to happen !!!!!!!!!!!!!!!!" << endl;
                        cout << "id_proc = " << M1_State.id_proc << "    " << "index_count = " << index_count << "  " << "size_rec = " << size_rec << endl;
                        exit(0);
                    }
                    
                    if (M1_State.display) {
                        if (M1_State.id_proc == DISPLAY_ID) {
                            cout << "i_E = " << i_E << "     "  << "i_f = " << i_f << "     "  << "i_Phi = " << i_Phi << "     "  << "i_Theta = " << i_Theta << endl;
                        }
                    }
                    
                    NLOPT_Optimization_Algorithm_NON_GRAY_M1_3D(NULL, rec_PM_local, index_count, &Var_index, num_points, M1_State.id_proc, M1_State, Mobius_Scale_Params, struc_map_N1, Sk_final, x_SH, y_SH, z_SH);   
                    
                    index_count++;
                }
            }
        }
    }
    
    // **********************************************************************************
    // Checks for MPI memory errors
    // **********************************************************************************
    if (size_rec_total != M1_State.num_proc*size_rec) {
        cout << "size_rec_total and maxium size_rec do not match !!!!!!!!!!!!!!!!" << endl;
        cout << "size_rec_total = " << size_rec_total << "  " << "max size_rec = " << size_rec*M1_State.num_proc << endl;
        exit(0);
    }
    
    int index, increment;
    if (M1_State.flag_use_mpi) {
        MPI_Barrier(MPI_COMM_WORLD);
        
        switch (M1_State.Max_Ent_Solution_Type) {
            case CLOSING_FLUX:
                cout << "Eddington factor must be computed in the 1D case" << endl;
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
        if (M1_State.id_proc == PRIMARY_ID) {
            for (int i_proc_E = 0; i_proc_E < num_proc_per_var_E; i_proc_E++) {
                for (int id_E = 0; id_E < incr_E; id_E++) {
                    for (int i_proc_f = 0; i_proc_f < num_proc_per_var_f; i_proc_f++) {
                        for (int id_f = 0; id_f < incr_f; id_f++) {
                            for (int id_Phi = 0 ; id_Phi < num_points->phi; id_Phi++) {
                                for (int id_Theta = 0 ; id_Theta < num_points->theta; id_Theta++) {
                                    index = i_proc_E*num_proc_per_var_f + i_proc_f;
                                    index = index*size_rec;
                                    index = index + id_E*incr_f + id_f;
                                    index = (index * num_points->phi + id_Phi) * num_points->theta + id_Theta;
                                    
                                    if ((i_proc_E*incr_E + id_E < num_points->E) && (i_proc_f*incr_f + id_f < num_points->f)) {   
                                        
                                        switch (M1_State.Max_Ent_Solution_Type) {
                                            case CLOSING_FLUX:
                                                cout << "Eddington factor must be computed in the 1D case" << endl;
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
            }
        }
        
        if (increment > num_points->E*num_points->f*num_points->phi*num_points->theta) {
            cout << "Value of increment exceeds maximum allowable value" << endl;
            cout << "increment = " << increment << " " << "num_points->E*num_points->f = " << num_points->E*num_points->f << endl;
        }
    } else {
        for (int id_E = 0; id_E < num_points->E; id_E++) {
            for (int id_f = 0; id_f < num_points->f; id_f++) {
                for (int id_Phi = 0 ; id_Phi < num_points->phi; id_Phi++) {
                    for (int id_Theta = 0 ; id_Theta < num_points->theta; id_Theta++) {
                        index = (id_E*num_points->f + id_f) * num_points->phi + id_Phi;
                        index = index * num_points->theta + id_Theta;
                        
                        switch (M1_State.Max_Ent_Solution_Type) {
                            case CLOSING_FLUX:
                                cout << "Eddington factor must be computed in the 1D case" << endl;
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
        }
    }
    
    Exiting:;
    
    if (x_SH != NULL) {
        delete[] x_SH;
    }
    if (y_SH != NULL) {
        delete[] y_SH;
    }
    if (z_SH != NULL) {
        delete[] z_SH;
    }
    
    delete[] rec_PM_local;
        
    if (M1_State.id_proc == DISPLAY_ID) {
        delete[] rec_PM_global_temp;
    }
    
    delete[] Sk_final;
}

void OPTIM_NON_GRAY_M1_3D(const M1_Var_Num_Points *num_points, M1_State_Param &M1_State, const int &N_pts_Mob_Scale, const long double *Coefficients_Fit_Mob_Scale, fstream &in_Partial_Moments_out, const Algebraic_Mapping_N1 *struc_map_N1) {
    M1_Var_Num_Points Var_index;
    
    // Declare and allocate array for storing precomputed solutions of the maximum entropy problem
    // for the partial angular moments
    record_Partial_Moments *rec_PM_global;
    if (M1_State.id_proc == PRIMARY_ID) {
        switch (M1_State.Max_Ent_Solution_Type) {
            case CLOSING_FLUX:
                cout << "Eddington factor must be computed in the 1D case" << endl;
                exit(0);
                break;
            case PARTIAL_MOMENTS:
                rec_PM_global = new record_Partial_Moments[num_points->E*num_points->f*num_points->phi*num_points->theta];
                break;
            default:
                cout << "Maximum entropy solution type not specified" << endl;
                exit(0);
                break;
        }
    }
    
    // Solve the numerical maximum entropy problem to compute the partial angular moments for a given
    // set of known angular moments up to first-order in 3D
    OPTIM_NON_GRAY_M1_3D_Array(NULL, rec_PM_global, num_points, M1_State, N_pts_Mob_Scale, Coefficients_Fit_Mob_Scale, struc_map_N1);
    
    
    int index;
    if (M1_State.id_proc == DISPLAY_ID) {
        for (int id_E = 0; id_E < num_points->E; id_E++) {
            for (int id_f = 0; id_f < num_points->f; id_f++) {
                for (int i_Phi = 0 ; i_Phi < num_points->phi; i_Phi++) {
                    for (int i_Theta = 0 ; i_Theta < num_points->theta; i_Theta++) {
                        index = (id_E*num_points->f + id_f) * num_points->phi + i_Phi;
                        index = index * num_points->theta + i_Theta;
                        
                        // cout << "E_NON_GRAY = " << rec_PM_global[index].I0 << "  " << "norm_f = " << rec_PM_global[index].N1 << "  " << "E_Plus_NON_GRAY = " << rec_PM_global[index].I0_Plus<< "  " << "N1_1_Plus_NON_GRAY = " << rec_PM_global[index].N1_1_Plus << endl;
                
                        switch (M1_State.Max_Ent_Solution_Type) {
                            case CLOSING_FLUX:
                                cout << "Eddington factor must be computed in the 1D case" << endl;
                                exit(0);
                                break;
                            case PARTIAL_MOMENTS:
                                write<record_Partial_Moments>(in_Partial_Moments_out, rec_PM_global[index]);
                                break;
                            default:
                                cout << "Maximum entropy solution type not specified" << endl;
                                exit(0);
                                break;
                        }
                    }
                }
            }
        }
        if (M1_State.display) {
            cout << "***************************** Writing to file Completed *****************************" << endl;
        }
    }
    
    Exiting:;
    
    if (M1_State.id_proc == PRIMARY_ID) {
        switch (M1_State.Max_Ent_Solution_Type) {
            case CLOSING_FLUX:
                cout << "Eddington factor must be computed in the 1D case" << endl;
                exit(0);
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

void generate_polynomials_3D(long double* poly, const int &Id_angle, const M1_State_Param &M1_State) {
    poly[0] = 1.0;
    poly[1] = M1_State.Omega1[Id_angle];
    poly[2] = M1_State.Omega2[Id_angle];
    poly[3] = M1_State.Omega3[Id_angle];
}

void generate_polynomials_3D(long double* poly, const long double &Omega1, const long double &Omega2, const long double &Omega3) {
    poly[0] = 1.0;
    poly[1] = Omega1;
    poly[2] = Omega2;
    poly[3] = Omega3;
}

long double generate_partial_moments_monomials_basis(const int &Id_angle, const M1_State_Param &M1_State) {
    long double poly;
    //     cout << "Omega1 = " << Omega1 << "     " << "Omega2 = " << Omega2 << "     " << "Omega3 = " << Omega3 << endl;
    
    switch (M1_State.Index_Partial_Moments_Derivatives) {
        case 0:
            poly = 1.0;
            break;
        case 1:
            poly = M1_State.Omega1[Id_angle];
            break;
        case 2:
            poly = M1_State.Omega2[Id_angle];
            break;
        case 3:
            poly = pow(M1_State.Omega1[Id_angle],2);
            break;
        case 4:
            poly = M1_State.Omega1[Id_angle]*M1_State.Omega2[Id_angle];
            break;
        case 5:
            poly = pow(M1_State.Omega2[Id_angle],2);
            break;
        default:
            cout << "Invalid value for Index_Partial_Moments_Derivatives !!!!!!!!!" << endl;
            exit(0);
            break;
    };
    
    return poly;
}

long double generate_partial_moments_monomials_basis(const long double &Omega1, const long double &Omega2, const long double &Omega3, const int &index) {
    long double poly;
    //     cout << "Omega1 = " << Omega1 << "     " << "Omega2 = " << Omega2 << "     " << "Omega3 = " << Omega3 << endl;
    
    switch (index) {
        case 0:
            poly = 1.0;
            break;
        case 1:
            poly = Omega1;
            break;
        case 2:
            poly = Omega2;
            break;
        case 3:
            poly = pow(Omega1,2);
            break;
        case 4:
            poly = Omega1*Omega2;
            break;
        case 5:
            poly = pow(Omega2,2);
            break;
    };
    
    return poly;
}

void Set_Moments_3D(M1_State_Param *M1_State, const M1_Var_Num_Points *Var_index, const M1_Var_Num_Points *num_points, const long double *x_SH, const long double *y_SH, const long double *z_SH, Mobius_Scale_Parameters *Mobius_Scale_Params, const Algebraic_Mapping_N1 *struc_map_N1) {
    long double f_test;
    long double fx_test, fy_test, fz_test;
    long double Theta, Phi;
    
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
        } else {
            if (num_points->f == 1) {
                f_test = 0.0;
            } else {
                // if we want to construct a polynomial of degree n, we need to consider
                // the root of the (n+1)^th order polynomial
                // if (M1_State->N1_1_VAL_DOMAIN == N1_1_EQUAL_ZERO) {
                // f_test = zeros_shifted( Var_index->f + num_points->f - 1, 2*(num_points->f - 1) + 1, -1.0, 1.0, M1_State->Node_Dist_f);
                // } else {
                // f_test = zeros_shifted( Var_index->f, num_points->f, 0.0, 1.0, M1_State->Node_Dist_f);
                //
                f_test = zeros_shifted( Var_index->f + num_points->f - 1, 2*(num_points->f - 1) + 1, -1.0, 1.0, M1_State->Node_Dist_f);
            }
        }
    }
    
    if (M1_State->N1_1_VAL_DOMAIN == N1_1_EQUAL_N1) {
        fx_test = f_test;
        fy_test = 0.0;
        fz_test = 0.0;
        
        M1_State->x_val = 1.0;
        M1_State->y_val = 0.0;
        M1_State->z_val = 0.0;
    } else if (M1_State->N1_1_VAL_DOMAIN == N1_1_EQUAL_ZERO) {
        fx_test = 0.0;
        fy_test = f_test;
        fz_test = 0.0;
        
        M1_State->x_val = 0.0;
        M1_State->y_val = 1.0;
        M1_State->z_val = 0.0;
    } else if (M1_State->N1_1_VAL_DOMAIN == FULL_3D) {
        if (num_points->Symmetry_Type == SYMMETRY_N1_1) {
            fx_test = f_test * x_SH[Var_index->phi * num_points->theta + Var_index->theta];
            fy_test = f_test * y_SH[Var_index->phi * num_points->theta + Var_index->theta];
            fz_test = f_test * z_SH[Var_index->phi * num_points->theta + Var_index->theta];
            
            M1_State->x_val = x_SH[Var_index->phi * num_points->theta + Var_index->theta];
            M1_State->y_val = y_SH[Var_index->phi * num_points->theta + Var_index->theta];
            M1_State->z_val = z_SH[Var_index->phi * num_points->theta + Var_index->theta];
        } else if (M1_State->Node_Dist_Phi == UNIFORM_DISTRIBUTION) {
            Phi = Uniform_Distribution(Var_index->phi, num_points->phi - 1, 0.0, 2.0*PI);
            Theta = Uniform_Distribution(Var_index->theta, num_points->theta - 1, 0.0, PI);
            fx_test = f_test*cos(Theta);
            fy_test = f_test*sin(Theta)*cos(Phi);
            fz_test = f_test*sin(Theta)*sin(Phi);
            
            M1_State->x_val = cos(Theta);
            M1_State->y_val = sin(Theta)*cos(Phi);
            M1_State->z_val = sin(Theta)*sin(Phi);
        } else if (M1_State->Node_Dist_Phi == SPHERICAL_HARMONIC_DISTRIBUTION) {
            fx_test = f_test * x_SH[Var_index->phi * num_points->theta + Var_index->theta];
            fy_test = f_test * y_SH[Var_index->phi * num_points->theta + Var_index->theta];
            fz_test = f_test * z_SH[Var_index->phi * num_points->theta + Var_index->theta];
            
            M1_State->x_val = x_SH[Var_index->phi * num_points->theta + Var_index->theta];
            M1_State->y_val = y_SH[Var_index->phi * num_points->theta + Var_index->theta];
            M1_State->z_val = z_SH[Var_index->phi * num_points->theta + Var_index->theta];
        }
    } else {
        cout << "N1_1 Domain Type not Specified!!!!!!!!!!!!!!!!!!!!!" << endl;
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
                M1_State->ratio_I0 = Uniform_Distribution(Var_index->E, num_points->E - 1, -1.0, 1.0);
            } else {
                M1_State->ratio_I0 = zeros_shifted(Var_index->E, num_points->E, -1.0, 1.0, M1_State->Node_Dist_E);
            }
        }
    }
    
    // Now compute the energy density based on the type of mapping and nodal distribution adopted
    long double r_N1;
    if (M1_State->Node_Dist_E != DISCRETE_SET_ENERGY) {
        if (M1_State->flag_finite_difference_N1) {
            // In the case where we are performing finite differencing, the other
            // variables must remain constant
            M1_State->MOM_VEC[0] = Inverse_Mobius_Transformation(M1_State->ratio_I0, Mobius_Scale_Params->Evaluate_Length_Scale(M1_State->f_test_knot));
        } else {
            if (M1_State->Node_Dist_f == UNIFORM_ALGEBRAIC_MAPPING_N1 || M1_State->Node_Dist_f == CHEBYSHEV_ALGEBRAIC_MAPPING_N1) {
                // M1_State->MOM_VEC[0] = Inverse_Mobius_Transformation(M1_State->ratio_I0, Mobius_Scale_Params->Evaluate_Length_Scale_r_N1(r_N1));
            } else if (struc_map_N1->flag_Algebraic_Map_N1) {
                // r_N1 = struc_map_N1->Mapping_L_Chi2(struc_map_N1->f_L_Chi2, f_test);
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
    
//     if (M1_State->MOM_VEC[0] < 0.0) {
//         M1_State->MOM_VEC[0] = 1.0;
//     }
    
    Set_Regime(M1_State);
    
    M1_State->I0 = M1_State->MOM_VEC[0];
    M1_State->N1 = f_test;
    
    M1_State->N1_1 = fx_test;
    M1_State->N1_2 = fy_test;
    M1_State->N1_3 = fz_test;
    
    M1_State->flag_Realizability = Check_Realizability(M1_State);
    
    M1_State->MOM_VEC[1] = fx_test*M1_State->MOM_VEC[0];
    M1_State->MOM_VEC[2] = fy_test*M1_State->MOM_VEC[0];
    M1_State->MOM_VEC[3] = fz_test*M1_State->MOM_VEC[0];
    
    if (M1_State->display) {
        cout  << "ratio_I0 = " << M1_State->ratio_I0 << "    "  << "I0 = " << M1_State->I0 << "    "  << "N1 = " << M1_State->N1 << "   " << "f_test = " << f_test << "   " << "fx_test = " << fx_test << "   " << "fy_test = " << fy_test << "   " << "fz_test = " << fz_test << endl;
    }
}

void add_constraints(nlopt_opt &opt, my_constraint_data *data, M1_State_Param &M1_State) {
    long double *poly;
    poly = new long double[M1_State.NVARS];
    int Id_Angle;
    
    for (int Id_angle_Phi = 0; Id_angle_Phi < 2*M1_State.order_quad; Id_angle_Phi++) {
        for (int Id_angle_mu = 0; Id_angle_mu < M1_State.order_quad; Id_angle_mu++) {
            Id_Angle = M1_State.Id_Angle(Id_angle_mu, Id_angle_Phi);
            generate_polynomials_3D(poly, Id_Angle, M1_State);
            
            data[Id_Angle].a = poly[0];
            data[Id_Angle].b = poly[1];
            data[Id_Angle].c = poly[2];
            data[Id_Angle].f = poly[3];
            nlopt_add_inequality_constraint_orthog(opt, myconstraint, &data[Id_Angle], 0.0);
        }
    }
    delete[] poly;
 }
 
void Quadrature_3D_Func(long double *Vals, const int &NFUN, max_ent_obj_grad_hess func, void *fdata) {
    long double Temp_vals[NFUN];
    long double epsi;
    long double mu_start, mu_end, phi_start, phi_end;
    int Index_Angle;
    M1_State_Param *M1_State = (M1_State_Param *) fdata;
    
    for (int i = 0; i < NFUN; i++) {
        Vals[i] = 0.0;
    }
    
    for (int i_phi_quad = 0; i_phi_quad < M1_State->N_intervals_phi_quad; i_phi_quad++) {
        phi_start = M1_State->bounds_quad_phi_refin[i_phi_quad];
        phi_end = M1_State->bounds_quad_phi_refin[i_phi_quad+1];
        
        // if (M1_State->refinement_level > 1)
        //    cout << "i_phi_quad = " << i_phi_quad << "  " << "N_intervals_phi_quad = " << M1_State->N_intervals_phi_quad << "  " << "phi_start = " << phi_start << "  " << "phi_end = " << phi_end << "  " << endl;
        
        for (int i_mu_quad = 0; i_mu_quad < M1_State->N_intervals_mu_quad; i_mu_quad++) {
            mu_start = M1_State->bounds_quad_mu_refin[i_mu_quad];
            mu_end = M1_State->bounds_quad_mu_refin[i_mu_quad+1];
            
            M1_State->compute_Omegas(mu_start, mu_end, phi_start, phi_end);
            
            // if (M1_State->refinement_level > 2)
            //   cout << "i_mu_quad = " << i_mu_quad << "  " << "N_intervals_mu_quad = " << M1_State->N_intervals_mu_quad << "  " << "mu_start = " << mu_start << "  " << "mu_end = " << mu_end << "  " << endl;
            
            for ( int Id_Phi = 0; Id_Phi < 2*M1_State->order_quad; Id_Phi++ ) {
                for ( int Id_mu = 0; Id_mu < M1_State->order_quad; Id_mu++ ) {
                    Index_Angle = M1_State->Id_Angle(Id_mu, Id_Phi);
                    func(Temp_vals, NFUN, Index_Angle, fdata);
                    for (int i = 0; i < NFUN; i++) {
                        Vals[i] = Vals[i] + M1_State->w_quad_3D[Index_Angle] * Temp_vals[i];
                    }
                }
            }
        }
    }
}

void Partial_Quadrature_3D_Func(long double *Vals, const int &NFUN, max_ent_obj_grad_hess func, void *fdata) {
    long double Temp_vals[NFUN];
    long double mu_start, mu_end, phi_start, phi_end;
    int Index_Angle;
    M1_State_Param *M1_State = (M1_State_Param *) fdata;
    
    for (int i = 0; i < NFUN; i++) {
        Vals[i] = 0.0;
    }
    
    for (int i_phi_quad = 0; i_phi_quad < M1_State->N_intervals_phi_quad; i_phi_quad++) {
        phi_start = M1_State->bounds_quad_phi_refin[i_phi_quad];
        phi_end = M1_State->bounds_quad_phi_refin[i_phi_quad+1];
        
        // if (M1_State->refinement_level > 1)
        //    cout << "i_phi_quad = " << i_phi_quad << "  " << "N_intervals_phi_quad = " << M1_State->N_intervals_phi_quad << "  " << "phi_start = " << phi_start << "  " << "phi_end = " << phi_end << "  " << endl;
        
        for (int i_mu_quad = 0; i_mu_quad < M1_State->N_intervals_mu_quad; i_mu_quad++) {
            mu_start = M1_State->bounds_quad_mu_refin[i_mu_quad];
            mu_end = M1_State->bounds_quad_mu_refin[i_mu_quad+1];
            
            M1_State->compute_Omegas(mu_start, mu_end, phi_start, phi_end);
            
            // if (M1_State->refinement_level > 2)
            //   cout << "i_mu_quad = " << i_mu_quad << "  " << "N_intervals_mu_quad = " << M1_State->N_intervals_mu_quad << "  " << "mu_start = " << mu_start << "  " << "mu_end = " << mu_end << "  " << endl;
            
            for ( int Id_Phi = 0; Id_Phi < 2*M1_State->order_quad; Id_Phi++ ) {
                for ( int Id_mu = 0; Id_mu < M1_State->order_quad; Id_mu++ ) {
                    Index_Angle = M1_State->Id_Angle(Id_mu, Id_Phi);
                    if (M1_State->Omega1[M1_State->Id_Angle(Id_mu, Id_Phi)] >= 0.0) {
                        func(Temp_vals, NFUN, Index_Angle, fdata);
                        for (int i = 0; i < NFUN; i++) {
                            if (M1_State->Omega1[M1_State->Id_Angle(Id_mu, Id_Phi)] <= 1.0e-12) {
                                Vals[i] = Vals[i] + 0.5 * M1_State->w_quad_3D[Id_Phi] * Temp_vals[i];
                            } else {
                                Vals[i] = Vals[i] + M1_State->w_quad_3D[Index_Angle] * Temp_vals[i];
                            }
                        }
                    }
                }
            }
        }
    }
}

void setup_spherical_harmonics_data (const int &N_points_Leg, long double *x_SH, long double *y_SH, long double *z_SH, long double *w_quad) {
    long double *phi_quad, *w_phi_quad, *mu_quad, *w_Leg;
    int index;
    phi_quad = new long double[2*N_points_Leg];
    w_phi_quad = new long double[2*N_points_Leg];
    mu_quad = new long double[N_points_Leg];
    w_Leg = new long double[N_points_Leg];
        
    circle_rule ( 2*N_points_Leg, w_phi_quad, phi_quad );
    cout << "Update this" << endl;
    // legendre_set ( N_points_Leg, mu_quad, w_Leg );
        
    for (int Id_Phi = 0; Id_Phi < 2*N_points_Leg; Id_Phi++) {
        for (int Id_mu = 0; Id_mu < N_points_Leg; Id_mu++) {
            index = Id_Phi*N_points_Leg + Id_mu;
            x_SH[index] = mu_quad[Id_mu];
            y_SH[index] = sqrt(1.0 - pow(mu_quad[Id_mu],2))*cos(phi_quad[Id_Phi]);
            z_SH[index] = sqrt(1.0 - pow(mu_quad[Id_mu],2))*sin(phi_quad[Id_Phi]);
            if (w_quad != NULL) {
                w_quad[index] = 2.0*PI*w_phi_quad[Id_Phi]*w_Leg[Id_mu];
            }
        }
    }
    
    delete[] phi_quad;
    delete[] mu_quad;
    delete[] w_phi_quad;
    delete[] w_Leg;
}
