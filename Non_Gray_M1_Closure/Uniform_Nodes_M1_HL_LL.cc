#include "../M1_Optimization/M1_Optimization.h"
#include "../M1_Closure_Interp/M1_Data_1D_Cheby.h"

using namespace nlopt;
using namespace std;

int main(int argc, char *argv[]) {
    char path_HL_out[256], path_LL_out[256], prefix[256], extension[256];
    record_Ncoeffs rec_Ncoeffs;
    int Problem_Type;
    double Length_Scale_Mobius = 0.0;
    double f_L_Chi2 = 0.0;
    int Node_Distribution_E, Node_Distribution_f;
    M1_Var_Num_Points Num_points;
    
    strcpy(path_HL_out, getenv(PATHVAR));
    strcat(path_HL_out, "/M1_Model/Non_Gray_M1_Closure/Uniform_Nodes_M1_HL");
    
    strcpy(path_LL_out, getenv(PATHVAR));
    strcat(path_LL_out, "/M1_Model/Non_Gray_M1_Closure/Uniform_Nodes_M1_LL");
    
    sprintf(extension, "_%.6d", 0);
   
    strcat(extension, ".dat");
    
    strcat(path_HL_out, extension);
    
    strcat(path_LL_out, extension);
    
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
    Node_Distribution_f = CHEBYSHEV_SECOND_KIND_DISTRIBUTION; //UNIFORM_ALGEBRAIC_MAPPING_N1;
    
    Num_points.E = 1;
    Num_points.f = 20;
    
    int N_pts_f_L_Chi2 = 1;
    
    Algebraic_Mapping_N1 struct_map_N1;
    
    // Setup parameters for algebraic mapping for the first-order normalized angular moment
    struct_map_N1.flag_Algebraic_Map_N1 = false;
    
    fstream in_HL_out;
    in_HL_out.open(path_HL_out, ios::out|ios::binary);
    
    if (!in_HL_out) {
        cout << "Uniform_Nodes_M1_HL.dat could not be accessed!" << endl;
    }
    
    Num_points.Maximum_Entropy_Solution_Regime = HYPERBOLIC_LIMIT;
    if (in_HL_out.good()) {
        cout << "Number of Points ..........................." << "N_Points_E = " << Num_points.E << "     "  << "N_Points_f = " << Num_points.f << endl;
            
        rec_Ncoeffs.N_Points_E = Num_points.E;
        rec_Ncoeffs.N_Points_f = Num_points.f;
        
        for (int id_f_L_Chi2 = 0; id_f_L_Chi2 < N_pts_f_L_Chi2; id_f_L_Chi2++) {
            f_L_Chi2 = f_L_Chi2_vals[id_f_L_Chi2];
            struct_map_N1.f_L_Chi2 = f_L_Chi2;
            
            OPTIM_NON_GRAY_M1(&Num_points, CLOSING_FLUX, Problem_Type, ONE_DIMENSIONAL, Node_Distribution_E, Node_Distribution_f, id_proc, num_proc, 1, &Length_Scale_Mobius, in_HL_out, &struct_map_N1, Inverse_mapping_L_Chi2, mapping_L_Chi2, true, false);
        }
        
        in_HL_out.close();
    }
         
    fstream in_LL_out;
    in_LL_out.open(path_LL_out, ios::out|ios::binary);
    
    if (!in_LL_out) {
        cout << "Uniform_Nodes_M1_LL.dat could not be accessed!" << endl;
    }
    
    Num_points.Maximum_Entropy_Solution_Regime = LOGARITHMIC_LIMIT;
    if (in_LL_out.good()) {
        cout << "Number of Points ..........................." << "N_Points_E = " << Num_points.E << "     "  << "N_Points_f = " << Num_points.f << endl;
            
        rec_Ncoeffs.N_Points_E = Num_points.E;
        rec_Ncoeffs.N_Points_f = Num_points.f;
        
        for (int id_f_L_Chi2 = 0; id_f_L_Chi2 < N_pts_f_L_Chi2; id_f_L_Chi2++) {
            f_L_Chi2 = f_L_Chi2_vals[id_f_L_Chi2];
            struct_map_N1.f_L_Chi2 = f_L_Chi2;
            
            OPTIM_NON_GRAY_M1(&Num_points, CLOSING_FLUX, Problem_Type, ONE_DIMENSIONAL, Node_Distribution_E, Node_Distribution_f, id_proc, num_proc, 1, &Length_Scale_Mobius, in_LL_out, &struct_map_N1, Inverse_mapping_L_Chi2, mapping_L_Chi2, true, false);
        }
        
        in_LL_out.close();
    }
}
