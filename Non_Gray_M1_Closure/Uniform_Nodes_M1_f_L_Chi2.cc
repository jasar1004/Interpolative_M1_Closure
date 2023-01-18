#include "../M1_Optimization/M1_Optimization.h"
#include "../M1_Closure_Interp/M1_Data_1D_Cheby.h"

using namespace nlopt;
using namespace std;

int main(int argc, char *argv[]) {
    char path1_out[256], prefix[256], extension[256];
    int Problem_Type;
    int id_max = 12;
    long double Length_Scale_Mobius = 0.0;
    double f_L_Chi2 = 0.0;
    int Node_Distribution_E, Node_Distribution_f;
    M1_Var_Num_Points Num_points;
    
    strcpy(path1_out, getenv(PATHVAR));
    strcat(path1_out, "/M1_Model/Non_Gray_M1_Closure/Uniform_Nodes_M1_f_L_Chi2");
    
    sprintf(extension, "_%.6d", 0);
   
    strcat(extension, ".dat");
    strcat(path1_out, extension);
    
    fstream in1_out;
    in1_out.open(path1_out, ios::out|ios::binary);
    
    if (!in1_out) {
        cout << "Uniform_Nodes_M1_f_L_Chi2.dat could not be accessed!" << endl;
    }
    
    int flag_MPI_Initialize = 0, flag_MPI_Finalize = 0;
    
    int ierr = 0, num_proc = 0, id_proc = 0;
    //     Initialize MPI.
    ierr = MPI_Init ( &argc, &argv );
    //     Get the number of processes.
    ierr = MPI_Comm_size ( MPI_COMM_WORLD, &num_proc );
    //     Get the ID of this process.
    ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id_proc );
    
    Problem_Type = NON_GRAY;
    Node_Distribution_E = CHEBYSHEV_SECOND_KIND_DISTRIBUTION;
    Node_Distribution_f = CHEBYSHEV_SECOND_KIND_DISTRIBUTION; //UNIFORM_DISTRIBUTION;
    
    Num_points.E = 20;
    Num_points.f = 20;
    
    Num_points.Length_Scale_Dist_Type = LENGTH_SCALE_DIST_UNIF;
    
    int N_pts_Mob_Scale = 400;
    int N_pts_f_L_Chi2 = 1;
    
    Algebraic_Mapping_N1 struct_map_N1;
    
    // Setup parameters for algebraic mapping for the first-order normalized angular moment
    struct_map_N1.flag_Algebraic_Map_N1 = false;
    struct_map_N1.f_L_Chi2 = 0.0;
    
    // Compute the nodal distribution for the length scale L_Chi2
    long double alpha_Lag;
    alpha_Lag = 0;
    long double x_Lag[N_pts_Mob_Scale], w_Lag[N_pts_Mob_Scale];
    gen_laguerre_ek_compute ( N_pts_Mob_Scale, alpha_Lag, x_Lag, w_Lag );
    
    if (in1_out.good()) {
        cout << "Number of Points ..........................." << "N_Points_E = " << Num_points.E << "     "  << "N_Points_f = " << Num_points.f << endl;
        
        for (int id_Length_Scale_Mobius = 0; id_Length_Scale_Mobius < N_pts_Mob_Scale; id_Length_Scale_Mobius++) {
            if (id_proc == DISPLAY_ID) {
                cout << "id_Length_Scale_Mobius = " << id_Length_Scale_Mobius << endl;
            }
            
            Length_Scale_Mobius = /*1000.0 +*/ x_Lag[id_Length_Scale_Mobius];
            
            // cout << "id_Length_Scale_Mobius = " << id_Length_Scale_Mobius << "  " << "Length_Scale_Mobius = " << Length_Scale_Mobius << endl;
            
            OPTIM_NON_GRAY_M1(&Num_points, CLOSING_FLUX, Problem_Type, ONE_DIMENSIONAL, Node_Distribution_E, Node_Distribution_f, id_proc, num_proc, 1, &Length_Scale_Mobius, in1_out, &struct_map_N1, true, true);
        }
        
        in1_out.close();
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
//     if (!done_already) {
        // Termminate MPI
        MPI_Finalize(); // The code was not working when I did not include this
        //  Terminate
//     }
}
