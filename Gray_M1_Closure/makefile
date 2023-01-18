OTHERDIR_PACKAGES=../../../Packages
OTHERDIR_NG_M1=../../Non_Gray_Model

Uniform_Nodes_M1: Uniform_Nodes_M1.o $(OTHERDIR_NG_M1)/Chebyshev_Interpolation/M1_Model_1D.o $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.o $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.o	
	  mpicc -Wall -g -I $(CFFC_Path)/nlopt/include Uniform_Nodes_M1.cc -L $(CFFC_Path)/nlopt/lib -o Uniform_Nodes_M1 $(OTHERDIR_NG_M1)/Chebyshev_Interpolation/M1_Model_1D.o $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.o $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.o -lnlopt -lm

Cheby_Nodes_M1: Cheby_Nodes_M1.o $(OTHERDIR_NG_M1)/Chebyshev_Interpolation/M1_Model_1D.o $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.o $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.o	
	  mpicc -Wall -g -I $(CFFC_Path)/nlopt/include Cheby_Nodes_M1.cc -L $(CFFC_Path)/nlopt/lib -o Cheby_Nodes_M1 $(OTHERDIR_NG_M1)/Chebyshev_Interpolation/M1_Model_1D.o $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.o $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.o -lnlopt -lm
	  
M1_Fit_Convergence_Node_Dist: $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.o $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.o
	  mpicc -Wall -g -I $(CFFC_Path)/nlopt/include M1_Fit_Convergence_Node_Dist.cc -L $(CFFC_Path)/nlopt/lib -o M1_Fit_Convergence_Node_Dist $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.o $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.o -lnlopt -lm
	       
Uniform_Nodes_M1.o: Uniform_Nodes_M1.cc
	       mpicc -Wall -g -c Uniform_Nodes_M1.cc
	       
Cheby_Nodes_M1.o: Cheby_Nodes_M1.cc
	       mpicc -Wall -g -c Cheby_Nodes_M1.cc
	       
M1_Model_1D.o: $(OTHERDIR_NG_M1)/Chebyshev_Interpolation/M1_Model_1D.cc
	       mpicc -Wall -g -c $(OTHERDIR_NG_M1)/Chebyshev_Interpolation/M1_Model_1D.cc
	       
chebyshev.o: $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.cpp
	mpicc -Wall -g -c $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.cpp
		
quadrule.o: $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.cpp
	mpicc -Wall -g -c $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.cpp

clean:
	rm -rf *.o
	rm -rf $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/*.o
	rm -rf $(OTHERDIR_PACKAGES)/Chebyshev/*.o
