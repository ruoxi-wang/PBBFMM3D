/*! \file	read_metadata.hpp
*/
#ifndef __read_metadata_hpp__
#define __read_metadata_hpp__

#include"environment.hpp"
#include"bbfmm.h"


using namespace std;


void read_Metadata(const string& filenameMetadata,double& L, int& interpolation_order, int& Ns, int& Nf, int& nCols, int& tree_level);


#endif //(__read_metadata_hpp__)
