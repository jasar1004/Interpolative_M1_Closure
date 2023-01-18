#include "../../Non_Gray_Model/Chebyshev_Interpolation/M1_Model_1D.h"

using namespace nlopt;
using namespace std;

int main() {
    char path1_out[256], prefix[256], extension[256];
    int Problem_Type;
    int Node_Distribution_E, Node_Distribution_f;
    M1_Var_Num_Points Num_points;
    
    int id = 0;
    
    strcpy(path1_out, getenv(PATHVAR));
    strcat(path1_out, "/M1_Model/Gray_Model/Chebyshev_Interpolation/Uniform_Nodes_M1");
    
    sprintf(extension, "_%.6d", id);
   
    strcat(extension, ".dat");
    strcat(path1_out, extension);
    
    fstream in1_out;
    in1_out.open(path1_out, ios::out|ios::binary);
    
    if (!in1_out) {
        cout << "Cheby_Nodes_M1.dat could not be accessed!" << endl;
    }
    
    Problem_Type = GRAY;
    Node_Distribution_E = UNIFORM_DISTRIBUTION;
    Node_Distribution_f = UNIFORM_DISTRIBUTION;
    
    Num_points.E = 1;
    Num_points.f = 100;
    
    if (in1_out.good()) {
        cout << "Number of Points ..........................." << "N_Points_E = " << Num_points.E << "     "  << "N_Points_f = " << Num_points.f << endl;
        OPTIM_NON_GRAY_M1(&Num_points, Problem_Type, Node_Distribution_E, Node_Distribution_f, id, in1_out);
        in1_out.close();
    }
}
