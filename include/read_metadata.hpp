/*! \file	read_metadata.hpp
*/
#ifndef __read_metadata_hpp__
#define __read_metadata_hpp__

#include"environment.hpp"
#include"bbfmm.h"


using namespace std;


void read_Metadata(const string& filenameMetadata,double& L, int& n, doft& dof, int& Ns, int& Nf, int& m, int& level, double& alpha);


#endif //(__read_metadata_hpp__)
