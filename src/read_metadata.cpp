/*!\file read_metadata.cpp
 */

#include"read_metadata.hpp"



void read_Metadata(const string& filenameMetadata,double& L, int& n, doft& dof, int& Ns, int& Nf, int& m, int& level, double& alpha) {
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
    ss >> L >> comma >> n >> comma >> dof.s >> comma >> dof.f >> comma >>
    Ns >> comma >> Nf >> comma >> m >> comma >> level >> comma >> alpha;
    fin.close();
}

