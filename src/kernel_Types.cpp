
/*!\file	kernel_Types.cpp
*/

#include"kernel_Types.hpp"

void kernel_Laplacian::setHomogen(string& kernelType, doft* dof) {
    homogen = -1;
    symmetry = 1;
    kernelType =  "Laplacian";
    dof->f = 1;
    dof->s = 1;
}

void kernel_Laplacian::EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                    double *K, doft *dof){
    vector3 diff;
	double rinv;

	
    // Compute 1/r
	diff.x = sourcepos.x - fieldpos.x;
	diff.y = sourcepos.y - fieldpos.y;
	diff.z = sourcepos.z - fieldpos.z;
	rinv   = 1./sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
	
    // Result = 1/r
	*K = rinv;
}



void kernel_LaplacianForce::setHomogen(string& kernelType, doft* dof) {
    homogen = -2;
    symmetry = -1;
    kernelType = "LaplacianForce";
    dof->f = 6;
    dof->s = 9;
}

void kernel_LaplacianForce::EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                       double *K, doft *dof){
    vector3 diff;
	double rinv;
	
    // Compute 1/r
	diff.x = sourcepos.x - fieldpos.x;
	diff.y = sourcepos.y - fieldpos.y;
	diff.z = sourcepos.z - fieldpos.z;
	rinv   = 1./sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);

    int idof, dof2 = dof->f*dof->s;
    double axis;
    for (idof=0; idof < dof2; idof++) {
        if (idof%3 == 0)
            axis = diff.x;
        else if (idof%3 == 1)
            axis = diff.y;
        else
            axis = diff.z;
        
        K[idof] = axis*axis*axis * rinv*rinv*rinv*rinv*rinv;
    }

}


void kernel_OneOverR4::setHomogen(string& kernelType, doft* dof) {
    homogen = -4;
    symmetry = 1;
    kernelType =  "OneOverR4";
    dof->f = 6;
    dof->s = 9;
}

void kernel_OneOverR4::EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                       double *K, doft *dof){
    vector3 diff;
	double rinv;
	
    // Compute 1/r
	diff.x = sourcepos.x - fieldpos.x;
	diff.y = sourcepos.y - fieldpos.y;
	diff.z = sourcepos.z - fieldpos.z;
	rinv   = 1./sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
	
    // Result = 1/r^4
	*K = rinv*rinv*rinv*rinv;
}

void kernel_Gaussian::setHomogen(string& kernelType, doft* dof) {
    homogen = 0;
    symmetry = 1;
    kernelType =  "Gaussian";
    dof->f = 1;
    dof->s = 1;
}

void kernel_Gaussian::EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                      double *K, doft *dof){
    vector3 diff;
	
    // Compute r
    diff.x = sourcepos.x - fieldpos.x;
    diff.y = sourcepos.y - fieldpos.y;
    diff.z = sourcepos.z - fieldpos.z;
    
    double r = sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
    *K = exp(-r*r);
}


void kernel_Logarithm::setHomogen(string& kernelType, doft* dof) {
    homogen = 0;
    symmetry = 1;
    kernelType =  "Logarithm";
    dof->f = 1;
    dof->s = 1;
}

void kernel_Logarithm::EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                     double *K, doft *dof){
    vector3 diff;
	
    diff.x = sourcepos.x - fieldpos.x;
    diff.y = sourcepos.y - fieldpos.y;
    diff.z = sourcepos.z - fieldpos.z;
	double rSquare = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
    if (rSquare == 0){
        *K = 0;
    }
    else{
        *K = 0.5*log(rSquare);
    }
    //*K = (sourcepos.x + 1) * sourcepos.y * sourcepos.x;
}

void kernel_OneOverR2::setHomogen(string& kernelType, doft* dof) {
    homogen = -2;
    symmetry = 1;
    kernelType =  "OneOverR2";
    dof->f = 1;
    dof->s = 1;
}

void kernel_OneOverR2::EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                      double *K, doft *dof){
    vector3 diff;
	double rinv;
	
    // Compute 1/r
	diff.x = sourcepos.x - fieldpos.x;
	diff.y = sourcepos.y - fieldpos.y;
	diff.z = sourcepos.z - fieldpos.z;
	rinv   = 1./sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
	
    // Result = 1/r^2
	*K = rinv*rinv;
}

void kernel_Quadric::setHomogen(string& kernelType, doft* dof) {
    homogen = 0;
    symmetry = 1;
    kernelType =  "Quadric";
    dof->f = 1;
    dof->s = 1;
}

void kernel_Quadric::EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                      double *K, doft *dof){
    vector3 diff;
	diff.x = sourcepos.x - fieldpos.x;
    diff.y = sourcepos.y - fieldpos.y;
    diff.z = sourcepos.z - fieldpos.z;
	double rSquare = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
    *K = 1.0+rSquare;
}

void kernel_InverseQuadric::setHomogen(string& kernelType, doft* dof) {
    homogen = 0;
    symmetry = 1;
    kernelType =  "InverseQuadric";
    dof->f = 1;
    dof->s = 1;
}

void kernel_InverseQuadric::EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                    double *K, doft *dof){
    vector3 diff;
	diff.x = sourcepos.x - fieldpos.x;
    diff.y = sourcepos.y - fieldpos.y;
    diff.z = sourcepos.z - fieldpos.z;
	double rSquare = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
    *K = 1.0/(1.0+rSquare);
}

void kernel_ThinPlateSpline::setHomogen(string& kernelType, doft* dof) {
    homogen = 0;
    symmetry = 1;
    kernelType =  "ThinPlateSpline";
    dof->f = 1;
    dof->s = 1;
}

void kernel_ThinPlateSpline::EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                           double *K, doft *dof){
    vector3 diff;
	diff.x = sourcepos.x - fieldpos.x;
    diff.y = sourcepos.y - fieldpos.y;
    diff.z = sourcepos.z - fieldpos.z;
	double rSquare = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
    if (rSquare == 0){
        *K = 0;
    }
    else{
        *K = 0.5*rSquare*log(rSquare);
    }

}


void kernel_Stokes::setHomogen(string& kernelType, doft* dof) {
    homogen = -1;
    symmetry = 1;
    kernelType = "Stokes";
    dof->f = 3;
    dof->s = 3;

}

void kernel_Stokes::EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                       double *K, doft *dof){
    double r[3];
    double rinv;
    
    r[0] = sourcepos.x - fieldpos.x;
    r[1] = sourcepos.y - fieldpos.y;
    r[2] = sourcepos.z - fieldpos.z;
    rinv   = 1./sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);

    // K is stored colum-wise    
    for (int col=0; col < dof->s; col++) {
        for (int row=0; row < dof->f; row++) {
            if (row == col)
                K[col*dof->f+row] = rinv + r[row]*r[col]*pow(rinv,3.0);
            else
                K[col*dof->f+row] = r[row]*r[col]*pow(rinv,3.0);
        }
    }
}






