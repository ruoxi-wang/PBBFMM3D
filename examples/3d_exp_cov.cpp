
#include "bbfmm3d.hpp"
#include <iostream>

void SetMetaData(double& L, int& interpolation_order, int& Ns, int& Nf, int& nCols, int& tree_level, double& eps, int& nx, int& ny, int& nz) {
//		read from parameters.txt
//    L   // Length of simulation cell (assumed to be a cube)
//    n   // Number of Chebyshev nodes per dimension
//    Ns  // Number of sources in simulation cell
//    Nf  // Number of target points in simulation cell

	ifstream fin;
	
	string filename = "../input/parameters.txt";

	fin.open(filename.c_str());
	
	if (!fin.good()){
		cerr << "Failed to open file " << filename << endl;
		throw runtime_error("Failed to open file!");
	}
    
	string line;
    getline(fin,line);
	getline(fin,line);
    
	line.erase(remove(line.begin(), line.end(), ' '),line.end());
    stringstream ss;
    ss << line;
    char comma;
    ss >> L >> comma >> interpolation_order >> comma >> nx >> comma >> ny 
		>> comma >> nz >> comma >> nCols >> comma >> tree_level >> comma >> eps;
    
	fin.close();

	Ns = nx*ny*nz;	// Number of sources in simulation cell
	Nf = nx*ny*nz;	// Number of target points in simulation cell
}

void read_xyz(const string& filenamex, int nx, const string& filenamey,int ny, const string& filenamez,int nz, vector<vector3>&xyz) {
	ifstream fin;
	double* x = new double[nx];
	double* y = new double[ny];
	double* z = new double[nz];
	
	/* Read source */
	fin.open(filenamex.c_str(),ios::binary);
	if (!fin.good()){
		cerr << "Failed to open file " << filenamex << endl;
		throw runtime_error("Failed to open file!");
	}
	for (int i=0;i<nx ;i++ )
	{
		fin >> x[i];
	}
	fin.close();

	fin.open(filenamey.c_str(),ios::binary);
	if (!fin.good()){
		cerr << "Failed to open file " << filenamex << endl;
		throw runtime_error("Failed to open file!");
	}
	for (int i=0;i<ny ;i++ )
	{
		fin >> y[i];
	}
	fin.close();

	fin.open(filenamez.c_str(),ios::binary);
	if (!fin.good()){
		cerr << "Failed to open file " << filenamex << endl;
		throw runtime_error("Failed to open file!");
	}
	for (int i=0;i<nz ;i++ )
	{
		fin >> z[i];
	}
	fin.close();

	for (int k=0;k < nz;k++)
	{
		for (int j=0;j < ny ;j++ )
		{
			for (int i=0;i < nx ;i++)
			{
				xyz[k*nx*ny + j*nx + i].x = x[i];
				xyz[k*nx*ny + j*nx + i].y = y[j];
				xyz[k*nx*ny + j*nx + i].z = z[k];
			}
		}
	}


}

/*
 * Function: SetSources
 * -------------------------------------------------------------------
 */

void SetSources(std::vector<vector3>& target, int Nf, std::vector<vector3>& source, int Ns, std::vector<double>& weight, int nCols,
                double L, int nx, int ny, int nz) {

	int l, i, k=0;
	
    // vector of ones for now
    for (l=0;l<nCols;l++) {
        for (i=0;i<Ns;i++, k++) {
                weight[k] = 1.0;
        }
    }

	read_xyz("../input/xcoord.txt",nx,"../input/ycoord.txt",ny,"../input/zcoord.txt",nz, source);

	for (i=0;i<Nf;i++) {
		target[i].x = source[i].x;
		target[i].y = source[i].y;
		target[i].z = source[i].z;
    }
}

class myKernel: public H2_3D_Tree {
public:
    myKernel(double L, int tree_level, int interpolation_order,  double epsilon, int use_chebyshev):H2_3D_Tree(L,tree_level,interpolation_order, epsilon, use_chebyshev){};
    virtual void SetKernelProperty() {
        homogen = 0;
        symmetry = 1;
        kernelType = "myKernel";
    }
    virtual double EvaluateKernel(vector3& targetpos, vector3& sourcepos) {
        vector3 diff;        
        // Compute exp(-r)
        diff.x = sourcepos.x - targetpos.x;
        diff.y = sourcepos.y - targetpos.y;
        diff.z = sourcepos.z - targetpos.z;
        double r = sqrt(diff.x*diff.x/100. + diff.y*diff.y/100. + diff.z*diff.z/2.25);
        return 0.75*exp(-r);
    }
};


int main(int argc, char *argv[]) {
    /**********************************************************/
    /*                                                        */
    /*              Initializing the problem                  */
    /*                                                        */
    /**********************************************************/
    
    double L;       // Length of simulation cell (assumed to be a cube)
    int interpolation_order;          // Number of Chebyshev nodes per dimension
    int Ns;         // Number of sources in simulation cell
    int Nf;         // Number of target points in simulation cell
    int nCols;
    int tree_level;
    double eps;
    int use_chebyshev = 0;
    int nx,ny,nz;

	// read data
    SetMetaData(L, interpolation_order, Ns, Nf, nCols, tree_level, eps, nx, ny, nz);

	std::vector<vector3> source(Ns); // Position array for the source points
	std::vector<vector3> target(Nf);  // Position array for the target points
	std::vector<double> weight(Ns*nCols); // Source array

	SetSources(target,Nf,source,Ns,weight,nCols, L,nx,ny,nz);



	//double err;
    std::vector<double> output(Nf*nCols);// Field array (BBFMM calculation)
    for (int i = 0; i < Nf*nCols; i++)
    	output[i] = 0; 				// TODO: check do we need to initialize
	cout << "L                  		  : " << L << endl;
	cout << "interpolation_order		  : " << interpolation_order << endl;
	cout << "Ns (=Nf)           		  : " << Ns << endl;
	cout << "Nf (=Nf)           		  : " << Nf << endl;
	cout << "nCols (# of cols in weight)  : " << nCols << endl;
	cout << "tree_level              	  : " << tree_level << endl;
    cout << "eps                		  : " << eps << endl;

    /**********************************************************/
    /*                                                        */
    /*                 Fast matrix vector product             */
    /*                                                        */
    /**********************************************************/

    /*****      Pre Computation     ******/
    // clock_t  t0 = clock();
	double t0 = omp_get_wtime(); 
	myKernel Atree(L,tree_level, interpolation_order, eps, use_chebyshev);
    Atree.buildFMMTree();
	// clock_t t1 = clock();
	double t1 = omp_get_wtime();
    
	double tPre = t1 - t0;
	cout << "Pre-computation finished" << endl;
    
    /*****      FMM Computation     *******/
    // t0 = clock();
    t0 = omp_get_wtime(); 
	H2_3D_Compute<myKernel> compute(Atree, target, source, weight, nCols, output);
	// t1 = clock();
	t1 = omp_get_wtime(); 
    double tFMM = t1 - t0;

	cout << "Fmm computation finished" << endl;

	/***  check accuracy ***/
	int num_rows = 100; // calculate first 100 rows.
    std::vector<double> output_dir(num_rows*nCols);

    t0 = omp_get_wtime(); 
    DirectCalc3D(&Atree, target, source, weight, nCols, output_dir, num_rows);
    t1 = omp_get_wtime(); 
    double tExact = t1 - t0;

    // Compute the 2-norm error
    double err = ComputeError(output,output_dir,Nf,nCols, num_rows);
    cout << "Relative Error for first " <<  num_rows  << " : " << err << endl;


	/*****      output result to binary file    ******/
    string outputfilename = "../output/output.bin";
    write_Into_Binary_File(outputfilename, &output[0], nCols*Nf);
    
	cout << "Pre-computation time: " << tPre << endl;
    cout << "FMM computing time:   " << tFMM  << endl;
	cout << "FMM total time:   "  << tPre+tFMM  << endl;
    
    /*******            Clean Up        *******/
    
 
    return 0;
}
