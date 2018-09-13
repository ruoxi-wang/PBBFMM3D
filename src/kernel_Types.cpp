
/*!\file	kernel_Types.cpp
*/

#include"kernel_Types.hpp"

void kernel_Laplacian::SetKernelProperty( ) {
    homogen = -1;
    symmetry = 1;
    kernelType =  "Laplacian";
    dof->f = 1;
    dof->s = 1;
}

double kernel_Laplacian::EvaluateKernel(vector3& targetpos, vector3& sourcepos){
    vector3 diff;
	double rinv;

	
    // Compute 1/r
	diff.x = sourcepos.x - targetpos.x;
	diff.y = sourcepos.y - targetpos.y;
	diff.z = sourcepos.z - targetpos.z;
	rinv   = 1./sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
	
    // Result = 1/r
	return rinv;
}





void kernel_OneOverR4::SetKernelProperty( ) {
    homogen = -4;
    symmetry = 1;
    kernelType =  "OneOverR4";
    dof->f = 6;
    dof->s = 9;
}

double kernel_OneOverR4::EvaluateKernel(vector3& targetpos, vector3& sourcepos){
    vector3 diff;
	double rinv;
	
    // Compute 1/r
	diff.x = sourcepos.x - targetpos.x;
	diff.y = sourcepos.y - targetpos.y;
	diff.z = sourcepos.z - targetpos.z;
	rinv   = 1./sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
	
    // Result = 1/r^4
	return rinv*rinv*rinv*rinv;
}

void kernel_Gaussian::SetKernelProperty( ) {
    homogen = 0;
    symmetry = 1;
    kernelType =  "Gaussian";
    dof->f = 1;
    dof->s = 1;
}

double kernel_Gaussian::EvaluateKernel(vector3& targetpos, vector3& sourcepos){
    vector3 diff;
	
    // Compute r
    diff.x = sourcepos.x - targetpos.x;
    diff.y = sourcepos.y - targetpos.y;
    diff.z = sourcepos.z - targetpos.z;
    
    double r = sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
    return exp(-r*r);
}


void kernel_Logarithm::SetKernelProperty( ) {
    homogen = 0;
    symmetry = 1;
    kernelType =  "Logarithm";
    dof->f = 1;
    dof->s = 1;
}

double kernel_Logarithm::EvaluateKernel(vector3& targetpos, vector3& sourcepos){
    vector3 diff;
	
    diff.x = sourcepos.x - targetpos.x;
    diff.y = sourcepos.y - targetpos.y;
    diff.z = sourcepos.z - targetpos.z;
	double rSquare = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
    if (rSquare == 0){
        return 0;
    }
    else{
        return 0.5*log(rSquare);
    }
    //*K = (sourcepos.x + 1) * sourcepos.y * sourcepos.x;
}

void kernel_OneOverR2::SetKernelProperty( ) {
    homogen = -2;
    symmetry = 1;
    kernelType =  "OneOverR2";
    dof->f = 1;
    dof->s = 1;
}

double kernel_OneOverR2::EvaluateKernel(vector3& targetpos, vector3& sourcepos){
    vector3 diff;
	double rinv;
	
    // Compute 1/r
	diff.x = sourcepos.x - targetpos.x;
	diff.y = sourcepos.y - targetpos.y;
	diff.z = sourcepos.z - targetpos.z;
	rinv   = 1./sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
	
    // Result = 1/r^2
	return rinv*rinv;
}

void kernel_Quadric::SetKernelProperty( ) {
    homogen = 0;
    symmetry = 1;
    kernelType =  "Quadric";
    dof->f = 1;
    dof->s = 1;
}

double kernel_Quadric::EvaluateKernel(vector3& targetpos, vector3& sourcepos){
    vector3 diff;
	diff.x = sourcepos.x - targetpos.x;
    diff.y = sourcepos.y - targetpos.y;
    diff.z = sourcepos.z - targetpos.z;
	double rSquare = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
    return 1.0+rSquare;
}

void kernel_InverseQuadric::SetKernelProperty( ) {
    homogen = 0;
    symmetry = 1;
    kernelType =  "InverseQuadric";
    dof->f = 1;
    dof->s = 1;
}

double kernel_InverseQuadric::EvaluateKernel(vector3& targetpos, vector3& sourcepos){
    vector3 diff;
	diff.x = sourcepos.x - targetpos.x;
    diff.y = sourcepos.y - targetpos.y;
    diff.z = sourcepos.z - targetpos.z;
	double rSquare = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
    return 1.0/(1.0+rSquare);
}

void kernel_ThinPlateSpline::SetKernelProperty( ) {
    homogen = 0;
    symmetry = 1;
    kernelType =  "ThinPlateSpline";
    dof->f = 1;
    dof->s = 1;
}

double kernel_ThinPlateSpline::EvaluateKernel(vector3& targetpos, vector3& sourcepos){
    vector3 diff;
	diff.x = sourcepos.x - targetpos.x;
    diff.y = sourcepos.y - targetpos.y;
    diff.z = sourcepos.z - targetpos.z;
	double rSquare = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
    if (rSquare == 0){
        return 0;
    }
    else{
        return 0.5*rSquare*log(rSquare);
    }

}







