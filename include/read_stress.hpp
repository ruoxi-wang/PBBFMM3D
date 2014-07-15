/*!\file	read_sources.cpp
 source file to field, source and charge information from binary files.
*/

#include"environment.hpp"
using namespace std;

void read_Stress(const string& filenameStress, double *stress_dir, const int& N) {
    ifstream fin;

    /* Read source */
	fin.open(filenameStress.c_str(),ios::binary);
	
	if (!fin.good()){
		cerr << "Failed to open file " << filenameStress << endl;
		throw runtime_error("Failed to open file!");
	}
    
   
    fin.read((char*) stress_dir, N*sizeof(double));
    fin.close();
}