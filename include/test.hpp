
/*!	\file test.hpp
  Defines different types of kernel
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
 * Computes the potential at the first field point and returns
 * the result in phi.
 */
template<typename T>
void DirectCalc3D(T* FMMtree, vector3 *field, int Nf, vector3 *source, double *q, int m,
                  int Ns, doft *dof, int lpbc, double L,
                  double *phi) {
    
    int i, j, l1, l2, l3;
    int dof_f = dof->f, dof_s = dof->s;
    vector3 sourcepos, cshift;
    char trans[] = "n";
    double alpha = 1;
    int incr = 1;
    double* Kij = new double[dof_f*dof_s];
    
    
    int dofNf = dof_f*Nf, dofNs = dof_s*Ns, begGlobal, k,l, count;
    double *KPBC = new double[dofNf*dofNs];
    
    int start = floor(0.5 + (pow(3,lpbc)-1)/2.);
    
    
    int end;
    if (lpbc == 0)
        end = -1;
    else if (lpbc == 1)
        end = 0;
    else
        end = 1;
    
    
    printf("Direct calculation starts from: %d to %d.\n", end+1, start);
    
    // Compute the interactions due to periodic image cells
    for (i=0;i<Nf;i++) {
        for (j=0;j<Ns;j++) {
            
            begGlobal = j*dofNf*dof_s + i*dof_f;
            
            for (l1=-start;l1<start+1;l1++) {
                cshift.x = (double)l1*L;
                for (l2=-start;l2<start+1;l2++) {
                    cshift.y = (double)l2*L;
                    for (l3=-start;l3<start+1;l3++) {
                        cshift.z = (double)l3*L;
                        if (abs(l1) > end || abs(l2) > end || abs(l3) > end) {
                            sourcepos.x = source[j].x + cshift.x;
                            sourcepos.y = source[j].y + cshift.y;
                            sourcepos.z = source[j].z + cshift.z;
                            FMMtree->EvaluateKernel(field[i],sourcepos, Kij,dof);
                            count = 0;
                            for (l=0; l<dof_s; l++)
                                for (k=0; k<dof_f; k++, count++)
                                    KPBC[begGlobal + l*dofNf + k] += Kij[count];
                            
                        }
                    }
                }
            }
        }
    }
    
    double beta = 0;
    for (int i = 0; i < m; i++) {
        dgemv_(trans,&dofNf,&dofNs,&alpha,KPBC,&dofNf,&q[i*Ns*dof->s],&incr,&beta,&phi[i*Nf*dof->f],&incr);
    }
    delete []Kij;
    Kij = NULL;
    //free(KPBC);
    delete []KPBC;
    KPBC=NULL;
}

double ComputeError(double *phi, double *phidir, int Nf, doft *dof, int m) {
    
    // L2 relative error
    int i, j, dof_f = dof->f;
    double sum1, sum2, tmp;
    
    sum1 = 0;
    sum2 = 0;
    
    for (i=0;i<Nf*m;i++) {
        for (j=0;j<dof_f;j++) {
            tmp = phi[i*dof_f+j]-phidir[i*dof_f+j];
            sum1 += tmp*tmp;
            sum2 += phidir[i*dof_f+j]*phidir[i*dof_f+j];
        }
    }
    
    return sqrt(sum1)/sqrt(sum2);
}

template <typename T>
double testInterplationErr(T* Atree, int Ns, int Nf) {
    Atree->buildFMMTree();
    int  i, j, k=0;
    
    vector3 source[Ns];    // Position array for the source points
    vector3 field[Nf];     // Position array for the field points
    double q[Ns*Atree->dof->s];  // Source array
    
    for (i=0;i<Ns;i++) {
        for (j=0; j<Atree->dof->s; j++, k++){
            q[k] = frand(0,1);
        }
    }
    double size = Atree->L/4;
    //generate random points in well seperated box.

    for (i=0;i<Ns;i++) {
		
		source[i].x = frand(-0.5,0.5)*size - 1.5*size;
		source[i].y = frand(-0.5,0.5)*size - 1.5*size;
		source[i].z = frand(-0.5,0.5)*size - 1.5*size;
	}
    
    // Randomly set field points
	for (i=0;i<Nf;i++) {
		field[i].x = frand(-0.5,0.5)*size + 0.5*size;
		field[i].y = frand(-0.5,0.5)*size + 0.5*size;
		field[i].z = frand(-0.5,0.5)*size + 0.5*size;
    }
    
    double *stress      =  new double[Nf*Atree->dof->f];// Field array (BBFMM calculation)
    double *stress_dir  =  new double[Nf*Atree->dof->f];// Field array (direct O(N^2) calculation)
  
    
    H2_3D_Compute<T> compute(Atree, field, source, Ns, Nf, q,1, stress);
    DirectCalc3D(Atree, field, Nf, source, q, 1, Ns, Atree->dof,0 ,0, stress_dir);
    double err = ComputeError(stress,stress_dir,Nf,Atree->dof,1);
    delete []stress;
    delete []stress_dir;
    cout << "Interplation error is: " << err << endl;
    return err;
}
#endif //(__test_hpp__)
