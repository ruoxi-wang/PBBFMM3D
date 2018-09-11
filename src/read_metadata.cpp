/*!\file read_metadata.cpp
 */

#include"read_metadata.hpp"



void read_Metadata(const string& filenameMetadata,double& L, int& interpolation_order, doft& dof, int& Ns, int& Nf, int& nCols, int& level) {
    ifstream fin;
	fin.open(filenameMetadata.c_str());
	if (!fin.good()){
		cerr << "Failed to open file " << filenameMetadata << endl;
		throw runtime_error("Failed to open file!");
	}
    string line;
    getline(fin,line);
    line.erase(remove(line.begin(), line.end(), ' '),
               line.end());
    stringstream ss;
    ss << line;
    char comma;
    ss >> L >> comma >> interpolation_order >> comma >> dof.s >> comma >> dof.f >> comma >>
    Ns >> comma >> Nf >> comma >> nCols >> comma >> level;
    fin.close();
}

