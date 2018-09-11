#include "bbfmm3d.hpp"
#include "stdint.h"
/*
 * Function: SetSources
 * -------------------------------------------------------------------
 * Distributes an equal number of positive and negative charges uniformly
 * in the simulation cell ([-0.5*L,0.5*L])^3 and take these same locations
 * as the target points.
 */

void SetSources(std::vector<vector3>& target, int Nf, std::vector<vector3>& source, int Ns, std::vector<double>& weight, int nCols,
                double L) {
	
	int l, i, j, k=0;
	
    // Distributes the sources randomly uniformly in a cubic cell
    for (l=0;l<nCols;l++) {
        for (i=0;i<Ns;i++) {
            weight[k] = frand(0,1);
            k += 1;
        }
    }
    
    for (i=0;i<Ns;i++) {
		source[i].x = frand(-0.5,0.5)*L;
		source[i].y = frand(-0.5,0.5)*L;
		source[i].z = frand(-0.5,0.5)*L;
	}
    
	for (i=0;i<Nf;i++) {
		target[i].x = source[i].x;
		target[i].y = source[i].y;
		target[i].z = source[i].z;
    }

}


class myKernel: public H2_3D_Tree {
public:
    myKernel(double L, int tree_level, int interpolation_order, double eps, int use_chebyshev):H2_3D_Tree(L, tree_level, interpolation_order, eps, use_chebyshev){};
    virtual void SetKernelProperty() {
        homogen = -1;
        symmetry = 1;
        kernelType = "myKernel";
    }
    virtual double EvaluateKernel(vector3& targetpos, vector3& sourcepos) {
        vector3 diff;        
        // Compute 1/r
        diff.x = sourcepos.x - targetpos.x;
        diff.y = sourcepos.y - targetpos.y;
        diff.z = sourcepos.z - targetpos.z;
        double r = sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
        return 1/r;
	}
};

int main(int argc, char *argv[]) {
    /**********************************************************/
    /*                                                        */
    /*              Initializing the problem                  */
    /*                                                        */
    /**********************************************************/
    double L;                   // Length of simulation cell (assumed to be a cube)
    int interpolation_order;    // Number of interpolation nodes per dimension
    int Ns;                     // Number of sources in simulation cell
    int Nf;                     // Number of targets in simulation cell
    int nCols;
    int tree_level;
    double eps;
    int use_chebyshev;

      // parse input arguments
    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            if (!strcmp(argv[i],"-n"))
                interpolation_order = atoi(argv[++i]);
            if (!strcmp(argv[i],"-N")){
                Ns = atof(argv[++i]);
                Nf = Ns;
            }
            if (!strcmp(argv[i],"-L"))
                L = atof(argv[++i]);
            if (!strcmp(argv[i],"-m"))
                nCols = atoi(argv[++i]);
            if (!strcmp(argv[i],"-level"))
                tree_level = atoi(argv[++i]);
            if (!strcmp(argv[i],"-eps"))
                eps = atof(argv[++i]);
            if (!strcmp(argv[i],"-cheb"))
                use_chebyshev = atoi(argv[++i]);
        }
    }

    std::vector<vector3> source(Ns); // Position array for the source points
    std::vector<vector3> target(Nf);  // Position array for the target points
    std::vector<double> weight(Ns*nCols); // Source array
    SetSources(target,Nf,source,Ns,weight,nCols,L);


    std::vector<double> output(Nf*nCols);   // Field array (BBFMM calculation)    
    
    for (int i = 0; i < Nf*nCols; i++)
        output[i] = 0;                          // TODO: check do we need to initialize
    
    /**********************************************************/
    /*                                                        */
    /*                 Fast matrix vector product             */
    /*                                                        */
    /**********************************************************/
    
    /*****      Pre Computation     ******/
    double t0 = omp_get_wtime();
    myKernel Atree(L,tree_level, interpolation_order, eps, use_chebyshev);
    Atree.buildFMMTree();
    double t1 = omp_get_wtime();
    double tPre = t1 - t0;
    
    /*****      FMM Computation     *******/
    t0 = omp_get_wtime();    
    H2_3D_Compute<myKernel> compute(Atree, target, source, weight, nCols, output);
    t1 = omp_get_wtime();    
    double tFMM = t1 - t0;

   /******      Check accuracy      *******/

    int num_rows = Nf;
    std::vector<double> output_dir(num_rows*nCols);// Field array (direct O(N^2) calculation)

    t0 = omp_get_wtime();   
    DirectCalc3D(&Atree, target, source, weight, nCols, output_dir, num_rows);
    t1 = omp_get_wtime();
    double tExact = t1 - t0;

    // Compute the 2-norm error
    double err = ComputeError(output,output_dir,Nf,nCols,num_rows);
    cout << "Relative Error for first " <<  num_rows  << " : " << err << endl;
  
    cout << "Pre-computation time: " << tPre << endl;
    cout << "FMM computing time:   " << tFMM  << endl;
    cout << "FMM total time:   " << tFMM+tPre  << endl;
   
    return 0;
}
