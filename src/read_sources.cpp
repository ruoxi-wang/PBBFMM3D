/*!\file	read_sources.cpp
 source file to target, source and charge information from binary files.
*/

#include"read_sources.hpp"
#include"environment.hpp"
using namespace std;

void read_Sources(const string& filenameField, vector3 *target, const int& Nf, const string& filenameSource, vector3 *source, const int& Ns, const string& filenameCharge, double *q, const int& m, const doft& dof) {
    ifstream fin;

    /* Read source */
	fin.open(filenameField.c_str(),ios::binary);
	
	if (!fin.good()){
		cerr << "Failed to open file " << filenameField << endl;
		throw runtime_error("Failed to open file!");
	}
    for (int k = 0; k < 3; k++) {
        for (int i = 0; i < Nf; i++) {
            if (k == 0 ) {
                fin.read((char*) &target[i].x, sizeof(double));
            }else if ( k == 1 ) {
                fin.read((char*) &target[i].y, sizeof(double));
            }else {
                fin.read((char*) &target[i].z, sizeof(double));
            }
            
        }
    }
	fin.close();
    
    /* Read target */
    fin.open(filenameSource.c_str(),ios::binary);
	if (!fin.good()){
		cerr << "Failed to open file " << filenameSource << endl;
		throw runtime_error("Failed to open file!");
	}
    for (int k = 0; k < 3; k++) {
        for (int i = 0; i < Ns; i++) {
            if (k == 0 ) {
                fin.read((char*) &source[i].x, sizeof(double));
            }else if ( k == 1 ) {
                fin.read((char*) &source[i].y, sizeof(double));
            }else {
                fin.read((char*) &source[i].z, sizeof(double));
            }
            
        }
    }
	fin.close();
    
    /* Read charge */
    fin.open(filenameCharge.c_str(),ios::binary);
	if (!fin.good()){
		cerr << "Failed to open file " << filenameCharge << endl;
		throw runtime_error("Failed to open file!");
	}
    fin.read((char*) q, m*Ns*dof.s*sizeof(double));
    fin.close();
}