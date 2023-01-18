OTHERDIR_PACKAGES=../../Packages
OTHERDIR_M1_INTERP=../M1_Closure_Interp
M1_OPTIM_DIR=../M1_Optimization
LFLAGS = -Wall -Wextra -g -liomp5 -lpthread

OBJS_ALGLIB = $(OTHERDIR_PACKAGES)/Alglib/src/alglibinternal.o $(OTHERDIR_PACKAGES)/Alglib/src/alglibmisc.o $(OTHERDIR_PACKAGES)/Alglib/src/ap.o $(OTHERDIR_PACKAGES)/Alglib/src/dataanalysis.o $(OTHERDIR_PACKAGES)/Alglib/src/diffequations.o $(OTHERDIR_PACKAGES)/Alglib/src/fasttransforms.o $(OTHERDIR_PACKAGES)/Alglib/src/integration.o $(OTHERDIR_PACKAGES)/Alglib/src/interpolation.o $(OTHERDIR_PACKAGES)/Alglib/src/linalg.o $(OTHERDIR_PACKAGES)/Alglib/src/optimization.o $(OTHERDIR_PACKAGES)/Alglib/src/solvers.o $(OTHERDIR_PACKAGES)/Alglib/src/specialfunctions.o $(OTHERDIR_PACKAGES)/Alglib/src/statistics.o

OBJS_M1_OPT = $(M1_OPTIM_DIR)/M1_Optimization_1D.o $(M1_OPTIM_DIR)/M1_Optimization_3D.o $(M1_OPTIM_DIR)/M1_Model_1D_Utilities.o $(M1_OPTIM_DIR)/Eddington_Factor_1D_Computations.o $(M1_OPTIM_DIR)/Partial_Moments_3D_Computations.o

OBJS_M1_INTERP = $(OTHERDIR_M1_INTERP)/M1_Data_1D_Cheby.o $(OTHERDIR_M1_INTERP)/M1_Data_1D_Cheby_Nested_Least_Squares.o $(OTHERDIR_M1_INTERP)/Eigenstructure_3D.o

OBJS_OTHERDIR_PACKAGES = $(OTHERDIR_PACKAGES)/Permutations/Permutations_With_Order.o $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.o $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.o $(OTHERDIR_PACKAGES)/Lebedev_Quadrature/sphere_lebedev_rule.o $(OTHERDIR_PACKAGES)/Circle_Quadrature/circle_rule.o $(OTHERDIR_PACKAGES)/Finite_Difference/differ.o

Uniform_Nodes_M1_HL_LL: Uniform_Nodes_M1_HL_LL.o $(OBJS_M1_OPT) $(OBJS_M1_INTERP) $(OBJS_OTHERDIR_PACKAGES)  	
	  mpicc $(LFLAGS) -I $(CFFC_Path)/nlopt/include Uniform_Nodes_M1_HL_LL.cc -L $(CFFC_Path)/nlopt/lib -o Uniform_Nodes_M1_HL_LL $(OBJS_M1_OPT) $(OBJS_M1_INTERP) $(OBJS_OTHERDIR_PACKAGES) -lnlopt -lm
	  
Uniform_Nodes_M1_f_L_Chi2: Uniform_Nodes_M1_f_L_Chi2.o $(OBJS_M1_OPT) $(OBJS_M1_INTERP) $(OBJS_OTHERDIR_PACKAGES)	
	  mpicc $(LFLAGS) -I $(CFFC_Path)/nlopt/include Uniform_Nodes_M1_f_L_Chi2.cc -L $(CFFC_Path)/nlopt/lib -o Uniform_Nodes_M1_f_L_Chi2 $(OBJS_M1_OPT) $(OBJS_M1_INTERP) $(OBJS_OTHERDIR_PACKAGES) -lnlopt -lm
	  
Uniform_Nodes_M1_E_Chebyshev_N1: Uniform_Nodes_M1_E_Chebyshev_N1.o $(OBJS_M1_OPT) $(OBJS_M1_INTERP) $(OBJS_OTHERDIR_PACKAGES)
	  mpicc $(LFLAGS) -I $(CFFC_Path)/nlopt/include Uniform_Nodes_M1_E_Chebyshev_N1.cc -L $(CFFC_Path)/nlopt/lib -o Uniform_Nodes_M1_E_Chebyshev_N1 $(OBJS_M1_OPT) $(OBJS_M1_INTERP) $(OBJS_OTHERDIR_PACKAGES) -lnlopt -lm
	  
Cheby_Nodes_M1: Cheby_Nodes_M1.o $(OBJS_M1_OPT) $(OBJS_M1_INTERP) $(OBJS_OTHERDIR_PACKAGES)
	  mpicc $(LFLAGS) -I $(CFFC_Path)/nlopt/include Cheby_Nodes_M1.cc -L $(CFFC_Path)/nlopt/lib -o Cheby_Nodes_M1 $(OBJS_M1_OPT) $(OBJS_M1_INTERP) $(OBJS_OTHERDIR_PACKAGES) -lnlopt -lm
	
Optimal_Mobius_Transformation: Optimal_Mobius_Transformation.o $(OBJS_M1_OPT) $(OBJS_M1_INTERP) $(OBJS_OTHERDIR_PACKAGES)
	  mpicc $(LFLAGS) -I $(CFFC_Path)/nlopt/include Optimal_Mobius_Transformation.cc -L $(CFFC_Path)/nlopt/lib -o Optimal_Mobius_Transformation $(OBJS_M1_OPT) $(OBJS_M1_INTERP) $(OBJS_OTHERDIR_PACKAGES) -lnlopt -lm
		  
Illustrations_Length_Scale: Illustrations_Length_Scale.o $(OBJS_M1_OPT) $(OBJS_M1_INTERP) $(OBJS_OTHERDIR_PACKAGES)
	  mpicc $(LFLAGS) -I $(CFFC_Path)/nlopt/include Illustrations_Length_Scale.cc -L $(CFFC_Path)/nlopt/lib -o Illustrations_Length_Scale $(OBJS_M1_OPT) $(OBJS_M1_INTERP) $(OBJS_OTHERDIR_PACKAGES) -lnlopt -lm
		  
Illustrations_Energy_Density: Illustrations_Energy_Density.o $(OBJS_M1_OPT) $(OBJS_M1_INTERP) $(OBJS_OTHERDIR_PACKAGES)
	  mpicc $(LFLAGS) -I $(CFFC_Path)/nlopt/include Illustrations_Energy_Density.cc -L $(CFFC_Path)/nlopt/lib -o Illustrations_Energy_Density $(OBJS_M1_OPT) $(OBJS_M1_INTERP) $(OBJS_OTHERDIR_PACKAGES) -lnlopt -lm
	  	     
Permutations_With_Order: $(OTHERDIR_PACKAGES)/Permutations/Permutations_With_Order.o $(M1_OPTIM_DIR)/M1_Optimization_1D.o $(M1_OPTIM_DIR)/M1_Optimization_3D.o $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.o $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.o $(OTHERDIR_PACKAGES)/Lebedev_Quadrature/sphere_lebedev_rule.o $(OTHERDIR_PACKAGES)/Circle_Quadrature/circle_rule.o	
	  mpicc $(LFLAGS) -I $(CFFC_Path)/nlopt/include $(OTHERDIR_PACKAGES)/Permutations/Permutations_With_Order.cc -L $(CFFC_Path)/nlopt/lib -o $(OTHERDIR_PACKAGES)/Permutations/Permutations_With_Order $(M1_OPTIM_DIR)/M1_Optimization_1D.o $(M1_OPTIM_DIR)/M1_Optimization_3D.o $(OTHERDIR_PACKAGES)/Permutations/Permutations_With_Order.o $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.o $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.o $(OTHERDIR_PACKAGES)/Lebedev_Quadrature/sphere_lebedev_rule.o $(OTHERDIR_PACKAGES)/Circle_Quadrature/circle_rule.o -lnlopt -lm
	       
Illustrations_Length_Scale.o: Illustrations_Length_Scale.cc
	       mpicc $(LFLAGS) -c Illustrations_Length_Scale.cc
	       
Illustrations_Energy_Density.o: Illustrations_Energy_Density.cc
	       mpicc $(LFLAGS) -c Illustrations_Energy_Density.cc
	       	       
Uniform_Nodes_M1_HL_LL.o: Uniform_Nodes_M1_HL_LL.cc
	       mpicc $(LFLAGS) -c Uniform_Nodes_M1_HL_LL.cc
	
Uniform_Nodes_M1_f_L_Chi2.o: Uniform_Nodes_M1_f_L_Chi2.cc
	       mpicc $(LFLAGS) -c Uniform_Nodes_M1_f_L_Chi2.cc
	       
Uniform_Nodes_M1_E_Chebyshev_N1.o:
	       mpicc $(LFLAGS) -c Uniform_Nodes_M1_E_Chebyshev_N1.cc

Optimal_Mobius_Transformation.o: Optimal_Mobius_Transformation.cc
	       mpicc $(LFLAGS) -c Optimal_Mobius_Transformation.cc

M1_Data_1D_Cheby.o: $(OTHERDIR_M1_INTERP)/M1_Data_1D_Cheby.cc
	       mpicc $(LFLAGS) -c $(OTHERDIR_M1_INTERP)/M1_Data_1D_Cheby.cc
	       	       
M1_Data_1D_Cheby_Nested_Least_Squares.o: $(OTHERDIR_M1_INTERP)/M1_Data_1D_Cheby_Nested_Least_Squares.cc
	       mpicc $(LFLAGS) -c $(OTHERDIR_M1_INTERP)/M1_Data_1D_Cheby_Nested_Least_Squares.cc
	       	       
Eigenstructure_3D.o: $(OTHERDIR_M1_INTERP)/Eigenstructure_3D.cc
	       mpicc $(LFLAGS) -c $(OTHERDIR_M1_INTERP)/Eigenstructure_3D.cc
	       
Cheby_Nodes_M1.o: Cheby_Nodes_M1.cc
	       mpicc $(LFLAGS) -c Cheby_Nodes_M1.cc
	       
M1_Optimization_1D.o: $(M1_OPTIM_DIR)/M1_Optimization_1D.cc
	       mpicc $(LFLAGS) -c $(M1_OPTIM_DIR)/M1_Optimization_1D.cc
	       
M1_Optimization_3D.o: $(M1_OPTIM_DIR)/M1_Optimization_3D.cc
	       mpicc $(LFLAGS) -c $(M1_OPTIM_DIR)/M1_Optimization_3D.cc
	       
M1_Model_1D_Utilities.o: $(M1_OPTIM_DIR)/M1_Model_1D_Utilities.cc
	       mpicc $(LFLAGS) -c $(M1_OPTIM_DIR)/M1_Model_1D_Utilities.cc
	       
Eddington_Factor_1D_Computations.o: $(M1_OPTIM_DIR)/Eddington_Factor_1D_Computations.cc
	       mpicc $(LFLAGS) -c $(M1_OPTIM_DIR)/Eddington_Factor_1D_Computations.cc
	       
Partial_Moments_3D_Computations.o: $(M1_OPTIM_DIR)/Partial_Moments_3D_Computations.cc
	       mpicc $(LFLAGS) -c $(M1_OPTIM_DIR)/Partial_Moments_3D_Computations.cc
	       
Permutations_With_Order.o: $(OTHERDIR_PACKAGES)/Permutations/Permutations_With_Order.cc
	       mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Permutations/Permutations_With_Order.cc
	       
Permutations_With_Order_Tests.o: $(OTHERDIR_PACKAGES)/Permutations/Permutations_With_Order_Tests.cc
	       mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Permutations/Permutations_With_Order_Tests.cc
	       
differ.o: $(OTHERDIR_PACKAGES)/Finite_Difference/differ.cpp
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Finite_Difference/differ.cpp
	       
chebyshev.o: $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.cpp
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Chebyshev/chebyshev.cpp
		
sphere_lebedev_rule.o: $(OTHERDIR_PACKAGES)/Lebedev_Quadrature/sphere_lebedev_rule.cpp
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Lebedev_Quadrature/sphere_lebedev_rule.cpp
		
quadrule.o: $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.cpp
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/quadrule.cpp

alglibinternal.o : $(OTHERDIR_PACKAGES)/Alglib/src/alglibinternal.cpp $(OTHERDIR_PACKAGES)/Alglib/src/alglibinternal.h  $(OTHERDIR_PACKAGES)/Alglib/src/stdafx.h
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Alglib/src/alglibinternal.cpp

alglibmisc.o : $(OTHERDIR_PACKAGES)/Alglib/src/alglibmisc.cpp $(OTHERDIR_PACKAGES)/Alglib/src/alglibmisc.h $(OTHERDIR_PACKAGES)/Alglib/src/stdafx.h
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Alglib/src/alglibmisc.cpp

ap.o : $(OTHERDIR_PACKAGES)/Alglib/src/alglibmisc.cpp $(OTHERDIR_PACKAGES)/Alglib/src/ap.h $(OTHERDIR_PACKAGES)/Alglib/src/stdafx.h
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Alglib/src/ap.cpp
	
dataanalysis.o : $(OTHERDIR_PACKAGES)/Alglib/src/dataanalysis.cpp $(OTHERDIR_PACKAGES)/Alglib/src/dataanalysis.h $(OTHERDIR_PACKAGES)/Alglib/src/stdafx.h 
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Alglib/src/dataanalysis.cpp
	
diffequations.o : $(OTHERDIR_PACKAGES)/Alglib/src/diffequations.cpp $(OTHERDIR_PACKAGES)/Alglib/src/diffequations.h $(OTHERDIR_PACKAGES)/Alglib/src/stdafx.h
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Alglib/src/diffequations.cpp
	
fasttransforms.o : $(OTHERDIR_PACKAGES)/Alglib/src/fasttransforms.cpp $(OTHERDIR_PACKAGES)/Alglib/src/fasttransforms.h $(OTHERDIR_PACKAGES)/Alglib/src/stdafx.h
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Alglib/src/fasttransforms.cpp
	
integration.o : $(OTHERDIR_PACKAGES)/Alglib/src/integration.cpp $(OTHERDIR_PACKAGES)/Alglib/src/integration.h $(OTHERDIR_PACKAGES)/Alglib/src/stdafx.h
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Alglib/src/integration.cpp
	
interpolation.o : $(OTHERDIR_PACKAGES)/Alglib/src/interpolation.cpp $(OTHERDIR_PACKAGES)/Alglib/src/interpolation.h $(OTHERDIR_PACKAGES)/Alglib/src/stdafx.h
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Alglib/src/interpolation.cpp
	
linalg.o : $(OTHERDIR_PACKAGES)/Alglib/src/linalg.cpp $(OTHERDIR_PACKAGES)/Alglib/src/linalg.h $(OTHERDIR_PACKAGES)/Alglib/src/stdafx.h
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Alglib/src/linalg.cpp
	
optimization.o : $(OTHERDIR_PACKAGES)/Alglib/src/optimization.cpp $(OTHERDIR_PACKAGES)/Alglib/src/optimization.h $(OTHERDIR_PACKAGES)/Alglib/src/stdafx.h
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Alglib/src/optimization.cpp
	
solvers.o : $(OTHERDIR_PACKAGES)/Alglib/src/solvers.cpp $(OTHERDIR_PACKAGES)/Alglib/src/solvers.h $(OTHERDIR_PACKAGES)/Alglib/src/stdafx.h
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Alglib/src/solvers.cpp
	
specialfunctions.o : $(OTHERDIR_PACKAGES)/Alglib/src/specialfunctions.cpp $(OTHERDIR_PACKAGES)/Alglib/src/specialfunctions.h $(OTHERDIR_PACKAGES)/Alglib/src/stdafx.h
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Alglib/src/specialfunctions.cpp
	
statistics.o : $(OTHERDIR_PACKAGES)/Alglib/src/statistics.cpp $(OTHERDIR_PACKAGES)/Alglib/src/statistics.h $(OTHERDIR_PACKAGES)/Alglib/src/stdafx.h
	mpicc $(LFLAGS) -c $(OTHERDIR_PACKAGES)/Alglib/src/statistics.cpp
	
clean:
	rm -rf *.o
	rm -rf $(M1_OPTIM_DIR)/*.o
	rm -rf $(OTHERDIR_PACKAGES)/Lebedev_Quadrature/*.o
	rm -rf $(OTHERDIR_PACKAGES)/Circle_Quadrature/*.o
	rm -rf $(OTHERDIR_PACKAGES)/One_Dimensional_Quadratures/*.o
	rm -rf $(OTHERDIR_PACKAGES)/Chebyshev/*.o
	rm -rf $(OTHERDIR_PACKAGES)/Finite_Difference/*.o
	rm -rf $(OTHERDIR_PACKAGES)/Permutations/*.o
	rm -rf $(OTHERDIR_PACKAGES)/Alglib/src/*.o
