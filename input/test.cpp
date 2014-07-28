//
//  test.cpp
//  
//
//  Created by Sivaram Ambikasaran on 7/22/14.
//
//

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;
int main() {
        srand(time(NULL));
        double L        =       1.0;
        unsigned nCheb  =       4;
        unsigned dof_s  =       3;
        unsigned dof_f  =       3;
        unsigned Ns     =       5000;
        unsigned Nf     =       5000;
        unsigned m      =       1;
        unsigned level  =       3;


        unsigned dim    =       3;
        double *q       = new double[Ns*dof_s*m];
        double *source  = new double[Ns*dim];
        double *field   = new double[Nf*dim];
        double randmax  =       L/RAND_MAX;

        for (unsigned k=0; k<dim*Ns; ++k) {
                source[k] = rand()*randmax - 0.5;
        }
        for (unsigned k=0; k<dim*Nf; ++k) {
                field[k] = rand()*randmax - 0.5;
        }
        for (unsigned k=0; k<Ns*dof_s*m; ++k) {
                q[k] = rand()*randmax;
        }
        std::ofstream myfile;

        myfile.open("source_test.bin",ios::binary);
        myfile.write((char *)source, Ns*dim*sizeof(double));
        myfile.close();

        myfile.open("field_test.bin",ios::binary);
        myfile.write((char *)field, Nf*dim*sizeof(double));
        myfile.close();

        myfile.open("charge_test.bin",ios::binary);
        myfile.write((char *)q, Ns*dof_s*m*sizeof(double));
        myfile.close();
        
        myfile.open("metadata_test.txt", std::ios::out | std::ios::binary);
        myfile << L << "," << nCheb << "," << dof_s << "," << dof_f << "," << Ns << "," << Nf << "," << m << "," << level << std::endl;
        myfile.close();

        delete []source;
        delete []field;
        delete []q;
}