
#include "bbfmm3d.hpp"


int main(int argc, char *argv[]) {
    /**********************************************************/
    /*                                                        */
    /*              Initializing the problem                  */
    /*                                                        */
    /**********************************************************/
    
    double L;       // Length of simulation cell (assumed to be a cube)
    int n;          // Number of Chebyshev nodes per dimension
    doft dof;
    int Ns;         // Number of sources in simulation cell
    int Nf;         // Number of field points in simulation cell
    int m;
    int level;
    double alpha;
    double eps = 1e-9 ;
    int use_chebyshev = 1;
    
    string filenameMetadata = "./../input/metadata_test.txt";
    read_Metadata(filenameMetadata, L, n, dof, Ns, Nf, m, level,alpha);
    vector3 source[Ns];    // Position array for the source points
    vector3 field[Nf];     // Position array for the field points
    double q[Ns*dof.s*m];  // Source array

    string filenameField  = "./../input/field_test.bin";
    string filenameSource = "./../input/source_test.bin";
    string filenameCharge = "./../input/charge_test.bin";
    read_Sources(filenameField,field,Nf,filenameSource,source,Ns,filenameCharge,q,m,dof);
    double err;
    double *stress      =  new double[Nf*dof.f*m];// Field array (BBFMM calculation)
    double *stress_dir  =  new double[Nf*dof.f*m];// Field array (direct O(N^2) calculation)
    

    
    /**********************************************************/
    /*                                                        */
    /*                 Fast matrix vector product             */
    /*                                                        */
    /**********************************************************/
    
    /*****      Pre Computation     ******/
    clock_t  t0 = clock();
    kernel_LaplacianForce Atree(&dof,L,level, n, eps, use_chebyshev);
    Atree.buildFMMTree();
    clock_t t1 = clock();
    double tPre = t1 - t0;

    /*****      FMM Computation     *******/
    /* this part can be repeated with different source, field and charges*/
    
    t0 = clock();
    H2_3D_Compute<kernel_LaplacianForce> compute(&Atree, field, source, Ns, Nf, q,m, stress);
    t1 = clock();
    double tFMM = t1 - t0;

    
    /*****      output result to binary file    ******/
    string outputfilename = "../output/stress.bin";
    write_Into_Binary_File(outputfilename, stress, m*Nf*dof.f);
    
    /**********************************************************/
    /*                                                        */
    /*              Exact matrix vector product               */
    /*                                                        */
    /**********************************************************/
    t0 = clock();
    DirectCalc3D(&Atree, field, Nf, source, q, m, Ns, &dof,0 , L, stress_dir);
    t1 = clock();
    double tExact = t1 - t0;
    
    cout << "Pre-computation time: " << double(tPre) / double(CLOCKS_PER_SEC) << endl;
    cout << "FMM computing time:   " << double(tFMM) / double(CLOCKS_PER_SEC)  << endl;
    cout << "Exact computing time: " << double(tExact) / double(CLOCKS_PER_SEC)  << endl;
    
    // Compute the 2-norm error
    err = ComputeError(stress,stress_dir,Nf,&dof,m);
    cout << "Relative Error: "  << err << endl;
    
    /*******            Clean Up        *******/
    
    delete []stress;
    delete []stress_dir;
    return 0;
}
