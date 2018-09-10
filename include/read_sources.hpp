/*!\file read_sources.hpp
 read sources from binary files
*/
#ifndef __read_sources_hpp__
#define __read_sources_hpp__

#include"environment.hpp"
#include"bbfmm.h"
using namespace std;

void read_Sources(const string& filenameField, vector3 *target, const int& Nf, const string& filenameSource, vector3 *source, const int& Ns, const string& filenameCharge, double *q, const int& m, const doft& dof);

#endif //(__read_sources_hpp__)
