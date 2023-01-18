#include "../../Non_Gray_Model/Chebyshev_Interpolation/M1_Model_1D.h"

using namespace nlopt;
using namespace std;

int main() {
    char path1_out[256], prefix[256], extension[256];
    record_Ncoeffs rec_Ncoeffs;
    int N_Coeffs_E, N_Coeffs_f;
    int id_max = 15;
    int Problem_Type;
    int Node_Distribution_E, Node_Distribution_f;
    M1_Var_Num_Points Num_points;
    
    int id = 0;
    
    strcpy(path1_out, getenv(PATHVAR));
    strcat(path1_out, "/M1_Model/Gray_Model/Chebyshev_Interpolation/Cheby_Nodes_M1");
    
    sprintf(extension, "_%.6d", id);
   
    strcat(extension, ".dat");
    strcat(path1_out, extension);
    
    fstream in1_out;
    in1_out.open(path1_out, ios::out|ios::binary);
    
    if (!in1_out) {
        cout << "Cheby_Nodes_M1.dat could not be accessed!" << endl;
    }
    
    Problem_Type = GRAY;
    Node_Distribution_f = CHEBYSHEV_FIRST_KIND_DISTRIBUTION;
    
    if (in1_out.good()) {
        for (int id_E = 1; id_E <= 1; id_E++) {
            for (int id_f = 1; id_f <= id_max; id_f++) {
                Num_points.E = id_E;
                Num_points.f = id_f;
                cout << "Set of Coefficients ..........................." << "N_Coeffs_f = " << Num_points.f << "     "  << "N_Coeffs_E = " << Num_points.E << endl;
                
                rec_Ncoeffs.N_Coeffs_E = id_E;
                rec_Ncoeffs.N_Coeffs_f = id_f;
                
                write<record_Ncoeffs>(in1_out, rec_Ncoeffs);
                
                OPTIM_NON_GRAY_M1(&Num_points, Problem_Type, Node_Distribution_E, Node_Distribution_f, id, in1_out);
            }
        }
        in1_out.close();
    }
}
