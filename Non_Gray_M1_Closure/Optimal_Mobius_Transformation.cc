#include "../M1_Closure_Interp/M1_Data_1D_Cheby.h"

int main(int argc, char *argv[]) {
    int N_pts_Mob_Scale, N_points_E, N_points_f;
    int N_Coeffs_E, N_Coeffs_f, N_pts_L_Chi2_LS;
    
    ofstream output_Opt_Coefficients;
    char path_Chi2_out[256], extension[256];
    
    int num_proc, id_proc, ierr, done_already;
    MPI_Initialized(&done_already);
    if (!done_already) {
        //     Initialize MPI.
        ierr = MPI_Init ( &argc, &argv );
        //     Get the number of processes.
        ierr = MPI_Comm_size ( MPI_COMM_WORLD, &num_proc );
        //     Get the ID of this process.
        ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id_proc );
    }
    
    if (id_proc == 0) {
        strcpy(path_Chi2_out, getenv(PATHVAR));
        strcat(path_Chi2_out, "/M1_Model/Non_Gray_M1_Closure/Coefficients_Chi2_BE.dat");
        output_Opt_Coefficients.open(path_Chi2_out);  
        
        if (!output_Opt_Coefficients) {
            cout << "Optimal_Coefficients_Chi2_NG_M1.dat could not be accessed!" << endl;    
        }
    }
    
    // Setup number of points
    N_pts_Mob_Scale = 400;
    N_points_E = 20;
    N_points_f = 20;
    
    // Setup static variables for M1_Data_1D_Cheby data class
    M1_Data_1D_Cheby :: Problem_Type = NON_GRAY;
    M1_Data_1D_Cheby :: flag_test_Implementation = false;
    M1_Data_1D_Cheby :: Implementation_type = IMPLEMENTATION_DERIVATIVES_N1;
    // M1_Data_1D_Cheby :: Implementation_type = IMPLEMENTATION_DERIVATIVES_N1_DERIVATIVES_R_I0;
    
    M1_Data_1D_Cheby Chi2_M1_1D_Uniform_f_L_Chi2(N_points_E, N_points_f, N_pts_Mob_Scale);
    M1_Data_1D_Cheby Chi2_M1_1D_Uniform_HL(1, N_points_f, 1);
    M1_Data_1D_Cheby Chi2_M1_1D_Uniform_LL(1, N_points_f, 1);
    
//     // Read Input for Hyperbolic limit
//     Chi2_M1_1D_Uniform_HL.OpenInputFile("Uniform_Nodes_M1_HL");
//     Chi2_M1_1D_Uniform_HL.ReadInputData();
//     Chi2_M1_1D_Uniform_HL.CloseInputFile();
//     Chi2_M1_1D_Uniform_HL.SetupInterpolant_Values_HL_LL();
//     
//     // Read Input for Logarithmic limit
//     Chi2_M1_1D_Uniform_LL.OpenInputFile("Uniform_Nodes_M1_LL");
//     Chi2_M1_1D_Uniform_LL.ReadInputData();
//     Chi2_M1_1D_Uniform_LL.CloseInputFile();
//     Chi2_M1_1D_Uniform_LL.SetupInterpolant_Values_HL_LL();
    
    Chi2_M1_1D_Uniform_f_L_Chi2.OpenInputFile("Uniform_Nodes_M1_f_L_Chi2");
    Chi2_M1_1D_Uniform_f_L_Chi2.ReadInputData();
    Chi2_M1_1D_Uniform_f_L_Chi2.CloseInputFile();
    Chi2_M1_1D_Uniform_f_L_Chi2.SetupInterpolant_Values_BE(Chi2_M1_1D_Uniform_HL, Chi2_M1_1D_Uniform_LL);
    Chi2_M1_1D_Uniform_f_L_Chi2.Write_Max_Ent_Data_Matlab();
    
    // exit(0);
    
    // Setup number of coefficients for the non-gray M1 closure polynomial interpolation
    N_Coeffs_E = 5;
    N_Coeffs_f = 5;
    
    // Setup number of coeffcients for interpolative-based approximation of optimal L_Chi2
    N_pts_L_Chi2_LS = N_Coeffs_f;
    
    // Declare data structure M1_Data_1D_Cheby for the purpose of the non-gray M1 closure 
    // polynomial interpolation
    M1_Data_1D_Cheby Chi2_M1_1D_Cheby(N_Coeffs_E, N_Coeffs_f, 1, N_pts_L_Chi2_LS);
    
    if (id_proc == 0) {
        cout  << "************************************************************************************"<< endl;
        cout  << "*************************** Performing Interpolation ************************"<< endl;
        cout  << "************************************************************************************"<< endl;
    }
    
    // Setup number of processors and id of current processor
    Chi2_M1_1D_Cheby.id_proc = id_proc;
    Chi2_M1_1D_Cheby.num_proc = num_proc;
    
    M1_1D_Data_Pointer M1_1D_Data;
    M1_1D_Data.M1_Data_Uniform_BE = &Chi2_M1_1D_Uniform_f_L_Chi2;
    M1_1D_Data.M1_Data_Uniform_HL = &Chi2_M1_1D_Uniform_HL;
    M1_1D_Data.M1_Data_Uniform_LL = &Chi2_M1_1D_Uniform_LL;
    
    // Compute optimal value for f_L_Chi2 as well as the coefficients D_L_Chi2 and D_Chi2
    Chi2_M1_1D_Cheby.Polynomial_Interpolation_Non_Gray_M1(M1_1D_Data, output_Opt_Coefficients);
    
    // Close output file for the coefficients
    if (id_proc == 0) {
        output_Opt_Coefficients.close();
    }
    
    Chi2_M1_1D_Uniform_f_L_Chi2.Write_Max_Ent_Data_Matlab();
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (!done_already) {
        // Termminate MPI
        MPI_Finalize(); // The code was not working when I did not include this
        //  Terminate
    }
}
