#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include "dir.hpp"
#include "fmm.hpp"

int main()
{
    /**********************************************************/
    /*                                                        */
    /*              Initializing the problem                  */
    /*                                                        */
    /**********************************************************/
    
    double L;                   // Length of simulation cell (assumed to be a cube)
    int interpolation_order;    // Number of interpolation nodes per dimension
    int Ns;                     // Number of sources in simulation cell
    int Nf;                     // Number of targes in simulation cell
    int nCols;
    int tree_level;
    double eps = 1e-5 ;
    int use_chebyshev = 1;
    
    string filenameMetadata = "./../input/metadata_test.txt";
    read_Metadata(filenameMetadata, L, interpolation_order, Ns, Nf, nCols, tree_level);
    std::vector<vector3> source(Ns); // Position array for the source points
    std::vector<vector3> target(Nf);  // Position array for the target points
    std::vector<double> weight(Ns*nCols); // Source array
   
    string filenameField  = "./../input/field_test.bin";
    string filenameSource = "./../input/source_test.bin";
    string filenameCharge = "./../input/charge_test.bin";
    read_Sources(filenameField,target,Nf,filenameField,source,Ns,filenameCharge,weight,nCols);
    std::vector<double> output(Nf*nCols);

  // begin solve
  Eigen::VectorXd b(Ns), x;
  b.setRandom();
  // Solve Ax = b using various iterative solver with direct calculation
  {
    // build tree
    kernel_Laplacian Atree(L, tree_level, interpolation_order, eps, use_chebyshev);
    Atree.buildFMMTree();
    MatrixDirect A;
    A.initialize(&Atree, target, source);
    Eigen::GMRES<MatrixDirect, Eigen::IdentityPreconditioner> gmres;
    gmres.compute(A);
    x = gmres.solve(b);
    std::cout << "Start GMRES with direct/exact mat-vec ..." << std::endl;
    std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
  }
  // Solve Ax = b using various iterative solver with FMM
  {
    MatrixFMM A;
    A.initialize(L, tree_level, interpolation_order, eps, use_chebyshev, target, source);
    Eigen::GMRES<MatrixFMM, Eigen::IdentityPreconditioner> gmres;
    gmres.compute(A);
    x = gmres.solve(b);
    std::cout << "Start GMRES with FMM/approximate mat-vec ..." << std::endl;
    std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
  }
}
