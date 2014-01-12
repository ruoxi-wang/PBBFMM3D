
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
    kernel_Laplacian(doft* dof, double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(dof,L,level,n, epsilon, use_chebyshev){};
    virtual void setHomogen(string& kernelType);
    virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};


/*! LaplacianForce kernel */
class kernel_LaplacianForce: public H2_3D_Tree {
public:
    kernel_LaplacianForce(doft* dof, double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(dof,L,level,n, epsilon, use_chebyshev){};
    virtual void setHomogen(string& kernelType);
    virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};

/*! OneOver4 kernel */
class kernel_OneOverR4: public H2_3D_Tree {
public:
    kernel_OneOverR4(doft* dof, double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(dof,L,level,n, epsilon, use_chebyshev){};
    virtual void setHomogen(string& kernelType);
    virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};

/*! Gaussian kernel */
class kernel_Gaussian: public H2_3D_Tree {
public:
    kernel_Gaussian(doft* dof, double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(dof,L,level,n,epsilon, use_chebyshev){};
    virtual void setHomogen(string& kernelType);
    virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};

/*! Polynomial kernel */
class kernel_Logarithm: public H2_3D_Tree {
public:
    kernel_Logarithm(doft* dof, double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(dof,L,level,n,epsilon, use_chebyshev){};
    virtual void setHomogen(string& kernelType);
    virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};

/*! OneOverR2 kernel */
class kernel_OneOverR2: public H2_3D_Tree {
public:
    kernel_OneOverR2(doft* dof, double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(dof,L,level,n,epsilon, use_chebyshev){};
    virtual void setHomogen(string& kernelType);
    virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};

/*! Quadric kernel */
class kernel_Quadric: public H2_3D_Tree {
public:
    kernel_Quadric(doft* dof, double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(dof,L,level,n,epsilon, use_chebyshev){};
    virtual void setHomogen(string& kernelType);
    virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};

/*! InverseQuadric kernel */
class kernel_InverseQuadric: public H2_3D_Tree {
public:
    kernel_InverseQuadric(doft* dof, double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(dof,L,level,n,epsilon, use_chebyshev){};
    virtual void setHomogen(string& kernelType);
    virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};

/*! ThinPlateSpline kernel */
class kernel_ThinPlateSpline: public H2_3D_Tree {
public:
    kernel_ThinPlateSpline(doft* dof, double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(dof,L,level,n,epsilon, use_chebyshev){};
    virtual void setHomogen(string& kernelType);
    virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                double *K, doft *dof);
};



#endif //(__kerbel_Types_hpp__)
