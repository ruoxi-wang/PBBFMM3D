
/*!	\file kernel_Types.hpp
  Defines different types of kernel
*/
#ifndef __kernel_Types_hpp__
#define __kernel_Types_hpp__

#include"H2_3D_Tree.hpp"
#include"environment.hpp"
/*! Laplacian kernel */
class kernel_Laplacian: public H2_3D_Tree {
public:
    kernel_Laplacian(double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(L,level,n, epsilon, use_chebyshev){};
     void setHomogen(string& kernelType,doft *dof);
     void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};


/*! LaplacianForce kernel */
class kernel_LaplacianForce: public H2_3D_Tree {
public:
    kernel_LaplacianForce(double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(L,level,n, epsilon, use_chebyshev){};
     void setHomogen(string& kernelType,doft *dof);
     void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};

/*! OneOver4 kernel */
class kernel_OneOverR4: public H2_3D_Tree {
public:
    kernel_OneOverR4(double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(L,level,n, epsilon, use_chebyshev){};
     void setHomogen(string& kernelType,doft *dof);
     void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};

/*! Gaussian kernel */
class kernel_Gaussian: public H2_3D_Tree {
public:
    kernel_Gaussian(double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(L,level,n,epsilon, use_chebyshev){};
     void setHomogen(string& kernelType,doft *dof);
     void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};

/*! Polynomial kernel */
class kernel_Logarithm: public H2_3D_Tree {
public:
    kernel_Logarithm(double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(L,level,n,epsilon, use_chebyshev){};
     void setHomogen(string& kernelType,doft *dof);
     void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};

/*! OneOverR2 kernel */
class kernel_OneOverR2: public H2_3D_Tree {
public:
    kernel_OneOverR2(double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(L,level,n,epsilon, use_chebyshev){};
     void setHomogen(string& kernelType,doft *dof);
     void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};

/*! Quadric kernel */
class kernel_Quadric: public H2_3D_Tree {
public:
    kernel_Quadric(double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(L,level,n,epsilon, use_chebyshev){};
     void setHomogen(string& kernelType,doft *dof);
     void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};

/*! InverseQuadric kernel */
class kernel_InverseQuadric: public H2_3D_Tree {
public:
    kernel_InverseQuadric(double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(L,level,n,epsilon, use_chebyshev){};
     void setHomogen(string& kernelType,doft *dof);
     void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};

/*! ThinPlateSpline kernel */
class kernel_ThinPlateSpline: public H2_3D_Tree {
public:
    kernel_ThinPlateSpline(double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(L,level,n,epsilon, use_chebyshev){};
     void setHomogen(string& kernelType,doft *dof);
     void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};

/*! Stokes kernel */
class kernel_Stokes: public H2_3D_Tree {
public:
    kernel_Stokes(double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(L,level,n, epsilon, use_chebyshev){};
     void setHomogen(string& kernelType,doft *dof);
     void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};



#endif //(__kerbel_Types_hpp__)
