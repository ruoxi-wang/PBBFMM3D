/*!\file read_sources.hpp
 read sources from binary files
*/
#ifndef __read_sources_hpp__
#define __read_sources_hpp__

#include"environment.hpp"
#include"bbfmm.h"
using namespace std;

void read_Sources(const string& filenameField, std::vector<vector3>& target, const int& Nf, const string& filenameSource, 
	std::vector<vector3>& source, const int& Ns, const string& filenameCharge, std::vector<double>& weight, const int& nCols);

#endif //(__read_sources_hpp__)
