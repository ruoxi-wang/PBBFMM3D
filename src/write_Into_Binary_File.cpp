
/*!\file	write_Into_Binary_File.cpp
*/

#include"write_Into_Binary_File.hpp"

using namespace std;

void write_Into_Binary_File(const string& filename, double* data, int numOfElems) {
    
    ofstream outdata;
	outdata.open(filename.c_str(),ios::binary);
    
    outdata.write((char *)data, numOfElems*sizeof(double));
	
	outdata.close();
}