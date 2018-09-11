
#include "bbfmm3d.hpp"

class myKernel: public H2_3D_Tree {
public:
    myKernel(double L, int tree_level, int interpolation_order, double epsilon,int use_chebyshev):H2_3D_Tree(L,tree_level,interpolation_order,epsilon,use_chebyshev){};
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
        return 1. / r;

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
    doft dof;
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
    double err;
    std::vector<double> output(Nf*nCols);


    
    /**********************************************************/
    /*                                                        */
    /*                 Fast matrix vector product             */
    /*                                                        */
    /**********************************************************/
    
    /*****      Pre Computation     ******/
    double t0 = omp_get_wtime(); 
    myKernel Atree(L, tree_level, interpolation_order, eps, use_chebyshev);
    Atree.buildFMMTree();
    double t1 = omp_get_wtime(); 
    double tPre = t1 - t0;
    /*****      FMM Computation     ******/
    /*(this part can be repeated with different source, target and charge)*/
    
    t0 = omp_get_wtime(); 
    H2_3D_Compute<myKernel> compute(Atree, target, source, weight, nCols, output);
    t1 = omp_get_wtime(); 
    double tFMM = t1 - t0;


    /*****      output result to binary file    ******/
    string outputfilename = "../output/output.bin";
    write_Into_Binary_File(outputfilename, &output[0], nCols*Nf);
    /**********************************************************/
    /*                                                        */
    /*              Exact matrix vector product               */
    /*                                                        */
    /**********************************************************/

    int num_rows = Nf; // calculate first 100 rows.
    std::vector<double> output_dir(num_rows*nCols);

    t0 = omp_get_wtime(); 
    DirectCalc3D(&Atree, target, source, weight, nCols, output_dir, num_rows);
    t1 = omp_get_wtime(); 
    double tExact = t1 - t0;

    // Compute the 2-norm error
    err = ComputeError(output,output_dir,Nf,nCols, num_rows);
    cout << "Relative Error for first " <<  num_rows  << " : " << err << endl;

    
    cout << "Pre-computation time: " << tPre << endl;
    cout << "FMM computing time:   " << tFMM << endl;
    cout << "Exact computing time: " << tExact  << endl;    
  
    return 0;
}
