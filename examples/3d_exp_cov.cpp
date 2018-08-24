
#include "bbfmm3d.hpp"
#include<iostream>

void SetMetaData(double& L, int& n, doft& dof, int& Ns, int& Nf, int& m, int& level, double& eps, int& nx, int& ny, int& nz) {
//		read from parameters.txt
//    L   // Length of simulation cell (assumed to be a cube)
//    n   // Number of Chebyshev nodes per dimension
//    Ns  // Number of sources in simulation cell
//    Nf  // Number of field points in simulation cell

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
    ss >> L >> comma >> n >> comma >> dof.s >> comma >> dof.f >> comma >> nx >> comma >> ny 
		>> comma >> nz >> comma >> m >> comma >> level >> comma >> eps;
    
	fin.close();

	Ns = nx*ny*nz;	// Number of sources in simulation cell
	Nf = nx*ny*nz;	// Number of field points in simulation cell
}

void read_xyz(const string& filenamex, int nx, const string& filenamey,int ny, const string& filenamez,int nz, vector3 *xyz) {
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

void SetSources(vector3 *field, int Nf, vector3 *source, int Ns, double *q, int m,
                doft *dof, double L, int nx, int ny, int nz) {

	int l, i, j, k=0;
	
    // vector of ones for now
    for (l=0;l<m;l++) {
        for (i=0;i<Ns;i++) {
            for (j=0; j<dof->s; j++, k++){
                q[k] = 1.0;
            }
        }
        
    }

//	// read from binary file
//	string filename = "../input/input.bin";
//	
//	ifstream fin;
//	
//	fin.open(filename.c_str(),ios::binary);
//	if (!fin.good()){
//		cerr << "Failed to open file " << filename << endl;
//		throw runtime_error("Failed to open file!");
//	}
//
//	fin.read((char*) q, m*Nf*dof->s*sizeof(double));
//	fin.close();	
//
	read_xyz("../input/xcoord_small.txt",nx,"../input/xcoord_small.txt",ny,"../input/xcoord_small.txt",nz, source);

	for (i=0;i<Nf;i++) {
		field[i].x = source[i].x;
		field[i].y = source[i].y;
		field[i].z = source[i].z;
    }
}

class myKernel: public H2_3D_Tree {
public:
    myKernel(double L, int level, int n,  double epsilon, int use_chebyshev):H2_3D_Tree(L,level,n, epsilon, use_chebyshev){};
    virtual void setHomogen(string& kernelType,doft *dof) {
        homogen = 0;
        symmetry = 1;
        kernelType = "myKernel";
        dof->s = 1;
        dof->f = 1;
    }
    virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos, // TODO: change copy to reference.
                                double *K, doft *dof) {
        vector3 diff;        
        // Compute exp(-r)
        diff.x = sourcepos.x - fieldpos.x;
        diff.y = sourcepos.y - fieldpos.y;
        diff.z = sourcepos.z - fieldpos.z;
        double r = sqrt(diff.x*diff.x/100. + diff.y*diff.y/100. + diff.z*diff.z/2.25);
        *K = 0.75*exp(-r);
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
    doft dof;
    int Ns;         // Number of sources in simulation cell
    int Nf;         // Number of field points in simulation cell
    int m;
    int level;
    double eps;
    int use_chebyshev = 1;
    int nx,ny,nz;

	// read data
    SetMetaData(L, n, dof, Ns, Nf, m, level, eps, nx, ny, nz);

	vector3* source = new vector3[Ns]; // Position array for the source points
	vector3* field = new vector3[Nf];  // Position array for the field points
	double* q = new double[Ns*dof.s*m]; // Source array

	SetSources(field,Nf,source,Ns,q,m,&dof,L,nx,ny,nz);
	//cout << "q[0] : " << q[0] << " q[end] : " << q[Ns-1] << endl;



	//double err;
    double *stress      =  new double[Nf*dof.f*m];// Field array (BBFMM calculation)
	
	cout << "L                  : " << L << endl;
	cout << "n (# chebyshev)    : " << n << endl;
	cout << "Ns (=Nf)           : " << Ns << endl;
	cout << "Nf (=Nf)           : " << Nf << endl;
	cout << "dof.s              : " << dof.s << endl;
	cout << "dof.f              : " << dof.f << endl;
	cout << "m (# of cols in q) : " << m << endl;
	cout << "level              : " << level << endl;
    cout << "eps                : " << eps << endl;

    
    /**********************************************************/
    /*                                                        */
    /*                 Fast matrix vector product             */
    /*                                                        */
    /**********************************************************/
    
    /*****      Pre Computation     ******/
    clock_t  t0 = clock();
    
	myKernel Atree(L,level, n, eps, use_chebyshev);
    Atree.buildFMMTree();
	clock_t t1 = clock();
    
	double tPre = t1 - t0;
	cout << "Pre-computation finished" << endl;
    
    /*****      FMM Computation     *******/
    t0 = clock();
	H2_3D_Compute<myKernel> compute(&Atree, field, source, Ns, Nf, q,m, stress);
	t1 = clock();
    double tFMM = t1 - t0;

	cout << "Fmm computation finished" << endl;

	/***  check accuracy ***/
	cout << "Checking accuracy" << endl;

	double stress_exact[10];
	for (int i = 0; i < 10; i++) {
		stress_exact[i] = 0;
		for (int j = 0; j < Ns; j++) {
			double val = 0;
			Atree.EvaluateKernel(field[i], source[j], &val, &dof);
			stress_exact[i] += val * q[j];
		}
	}
	double diff = 0;
	double sum_Ax = 0;
	for (int i=0;i<10 ;i++ )
	{	
		diff += (stress[i] - stress_exact[i])*(stress[i] - stress_exact[i]); 
		sum_Ax += stress_exact[i] * stress_exact[i];
	}
	cout << "diff of first 10 values = " << sqrt(diff) / sqrt(sum_Ax) << endl;

	/*****      output result to binary file    ******/
    string outputfilename = "../output/stress.bin";
    write_Into_Binary_File(outputfilename, stress, m*Nf*dof.f);
    
    //skip Exact matrix vector product

	cout << "Pre-computation time: " << double(tPre) / double(CLOCKS_PER_SEC) << endl;
    cout << "FMM computing time:   " << double(tFMM) / double(CLOCKS_PER_SEC)  << endl;
	cout << "FMM total time:   "  << double(tPre+tFMM) / double(CLOCKS_PER_SEC)  << endl;
    
    /*******            Clean Up        *******/
    
    delete []stress;
	delete []source;
	delete []field;
    return 0;
}
