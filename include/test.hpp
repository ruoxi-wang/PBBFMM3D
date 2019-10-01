
/*!	\file test.hpp
*/
#ifndef __test_hpp__
#define __test_hpp__

#include "bbfmm.h"
#include"environment.hpp"
#include"compute.hpp"

//#include<iostream>
/*
 * Function: DirectCalc3D
 * ---------------------------------------------------------------------
 * Computes the potential at the first target point and returns
 * the result in output.
 */
template<typename T>
void DirectCalc3D(T* FMMtree, const std::vector<vector3>& target, const std::vector<vector3>& source, std::vector<double>& weight, int nCols,
                  std::vector<double>& output_exact, int num_rows) {

    int Ns = source.size();
    int Nf = target.size();
    for (int k = 0; k < nCols; k++)
        for (int i = 0; i < num_rows; i++) {
            output_exact[k*num_rows+i] = 0;
            for (int j = 0; j < Ns; j++) {
                double val = FMMtree->EvaluateKernel(target[i], source[j]);
                if (isinf(val))
                    val = 0;
                output_exact[k*num_rows+i] += val * weight[k*Nf+j];
            }
        }
}

double ComputeError(std::vector<double>& output, std::vector<double>& outputdir, int Nf, int nCols, int num_rows) {
    
    double diff = 0;
    double sum_Ax = 0;
    for (int k = 0; k < nCols; k++)
        for (int i=0;i<num_rows ;i++ )
        {   
            diff += (output[k*Nf+i] - outputdir[k*num_rows+i])*(output[k*Nf+i] - outputdir[k*num_rows+i]); 
            sum_Ax += outputdir[k*num_rows+i] * outputdir[k*num_rows+i];
        }
    return sqrt(diff) / sqrt(sum_Ax);   
}

template <typename T>
double testInterplationErr(T* Atree, int Ns, int Nf) {

    Atree->buildFMMTree();
    int  i, j, k=0;
    vector3 source[Ns];    // Position array for the source points
    vector3 target[Nf];     // Position array for the target points
    double weight[Ns*Atree->dof->s];  // Source array
    
    for (i=0;i<Ns;i++) {
        for (j=0; j<Atree->dof->s; j++, k++){
            weight[k] = frand(0,1);
        }
    }
    double size = Atree->L/4;
    //generate random points in well seperated box.

    for (i=0;i<Ns;i++) {
		
		source[i].x = frand(-0.5,0.5)*size - 1.5*size;
		source[i].y = frand(-0.5,0.5)*size - 1.5*size;
		source[i].z = frand(-0.5,0.5)*size - 1.5*size;
	}
    
    // Randomly set target points
	for (i=0;i<Nf;i++) {
		target[i].x = frand(-0.5,0.5)*size + 0.5*size;
		target[i].y = frand(-0.5,0.5)*size + 0.5*size;
		target[i].z = frand(-0.5,0.5)*size + 0.5*size;
    }
    
    std::vector<double> output(Nf*Atree->dof->f);// Field array (BBFMM calculation)
    std::vector<double> output_dir(Nf*Atree->dof->f);// Field array (direct O(N^2) calculation)

    
    H2_3D_Compute<T> compute(Atree, target, source, Ns, Nf, weight,1, output);

    DirectCalc3D(Atree, target, source, weight, 1, output_dir, Nf);

    double err = ComputeError(output,output_dir,Nf, 1, Nf);
   
    cout << "Interplation error is: " << err << endl;
    return err;
}
#endif //(__test_hpp__)
