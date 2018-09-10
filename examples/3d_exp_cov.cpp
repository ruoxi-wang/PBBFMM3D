
#include "bbfmm3d.hpp"
#include<iostream>

void SetMetaData(double& L, int& n, int& Ns, int& Nf, int& m, int& level, double& eps, int& nx, int& ny, int& nz) {
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
    ss >> L >> comma >> n >> comma >> nx >> comma >> ny 
		>> comma >> nz >> comma >> m >> comma >> level >> comma >> eps;
    
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

void SetSources(std::vector<vector3>& target, int Nf, std::vector<vector3>& source, int Ns, std::vector<double>& q, int m,
                double L, int nx, int ny, int nz) {

	int l, i, k=0;
	
    // vector of ones for now
    for (l=0;l<m;l++) {
        for (i=0;i<Ns;i++, k++) {
                q[k] = 1.0;
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
    myKernel(double L, int level, int n,  double epsilon, int use_chebyshev):H2_3D_Tree(L,level,n, epsilon, use_chebyshev){};
    virtual void setKernelProperty() {
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
    int n;          // Number of Chebyshev nodes per dimension
    int Ns;         // Number of sources in simulation cell
    int Nf;         // Number of target points in simulation cell
    int m;
    int level;
    double eps;
    int use_chebyshev = 0;
    int nx,ny,nz;

	// read data
    SetMetaData(L, n, Ns, Nf, m, level, eps, nx, ny, nz);

	std::vector<vector3> source(Ns); // Position array for the source points
	std::vector<vector3> target(Nf);  // Position array for the target points
	std::vector<double> q(Ns*m); // Source array

	SetSources(target,Nf,source,Ns,q,m, L,nx,ny,nz);
	//cout << "q[0] : " << q[0] << " q[end] : " << q[Ns-1] << endl;



	//double err;
    std::vector<double> output(Nf*m);// Field array (BBFMM calculation)
    for (int i = 0; i < Nf*m; i++)
    	output[i] = 0; 				// TODO: check do we need to initialize
	cout << "L                  : " << L << endl;
	cout << "n (# chebyshev)    : " << n << endl;
	cout << "Ns (=Nf)           : " << Ns << endl;
	cout << "Nf (=Nf)           : " << Nf << endl;
	cout << "m (# of cols in q) : " << m << endl;
	cout << "level              : " << level << endl;
    cout << "eps                : " << eps << endl;

    /**********************************************************/
    /*                                                        */
    /*                 Fast matrix vector product             */
    /*                                                        */
    /**********************************************************/

    /*****      Pre Computation     ******/
    // clock_t  t0 = clock();
	double t0 = omp_get_wtime(); 
	myKernel Atree(L,level, n, eps, use_chebyshev);
    Atree.buildFMMTree();
	// clock_t t1 = clock();
	double t1 = omp_get_wtime();
    
	double tPre = t1 - t0;
	cout << "Pre-computation finished" << endl;
    
    /*****      FMM Computation     *******/
    // t0 = clock();
    t0 = omp_get_wtime(); 
	H2_3D_Compute<myKernel> compute(Atree, target, source, Ns, Nf, q, m, output);
	// t1 = clock();
	t1 = omp_get_wtime(); 
    double tFMM = t1 - t0;

	cout << "Fmm computation finished" << endl;

	/***  check accuracy ***/
	int num_rows = 10;
	cout << "Checking accuracy of first " << num_rows << " values" << endl;

	double output_exact[num_rows*m];
	for (int k = 0; k < m; k++)
		for (int i = 0; i < num_rows; i++) {
			output_exact[k*num_rows+i] = 0;
			for (int j = 0; j < Ns; j++) {
				double val = Atree.EvaluateKernel(target[i], source[j]);
				output_exact[k*num_rows+i] += val * q[k*Nf+j];
			}
		}

	double diff = 0;
	double sum_Ax = 0;
	for (int k = 0; k < m; k++)
		for (int i=0;i<num_rows ;i++ )
		{	
			diff += (output[k*Nf+i] - output_exact[k*num_rows+i])*(output[k*Nf+i] - output_exact[k*num_rows+i]); 
			sum_Ax += output_exact[k*num_rows+i] * output_exact[k*num_rows+i];
		}
	cout << "diff of first 10 values = " << sqrt(diff) / sqrt(sum_Ax) << endl;

	/*****      output result to binary file    ******/
    string outputfilename = "../output/output.bin";
    write_Into_Binary_File(outputfilename, &output[0], m*Nf);
    
	cout << "Pre-computation time: " << tPre << endl;
    cout << "FMM computing time:   " << tFMM  << endl;
	cout << "FMM total time:   "  << tPre+tFMM  << endl;
    
    /*******            Clean Up        *******/
    
 
    return 0;
}
