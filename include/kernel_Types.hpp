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
    kernel_Laplacian(double L, int tree_level, int interpolation_order, double epsilon, int use_chebyshev):H2_3D_Tree(L,tree_level,interpolation_order, epsilon, use_chebyshev){};
     void SetKernelProperty();
     double EvaluateKernel(const vector3& targetpos, const vector3& sourcepos);
};


/*! LaplacianForce kernel */
class kernel_LaplacianForce: public H2_3D_Tree {
public:
    kernel_LaplacianForce(double L, int tree_level, int interpolation_order, double epsilon, int use_chebyshev):H2_3D_Tree(L,tree_level,interpolation_order, epsilon, use_chebyshev){};
     void SetKernelProperty();
     double EvaluateKernel(vector3& targetpos, vector3& sourcepos);
};

/*! OneOver4 kernel */
class kernel_OneOverR4: public H2_3D_Tree {
public:
    kernel_OneOverR4(double L, int tree_level, int interpolation_order, double epsilon, int use_chebyshev):H2_3D_Tree(L,tree_level,interpolation_order, epsilon, use_chebyshev){};
     void SetKernelProperty();
     double EvaluateKernel(vector3& targetpos, vector3& sourcepos);
};

/*! Gaussian kernel */
class kernel_Gaussian: public H2_3D_Tree {
public:
    kernel_Gaussian(double L, int tree_level, int interpolation_order, double epsilon, int use_chebyshev):H2_3D_Tree(L,tree_level,interpolation_order,epsilon, use_chebyshev){};
     void SetKernelProperty();
     double EvaluateKernel(vector3& targetpos, vector3& sourcepos);
};

/*! Polynomial kernel */
class kernel_Logarithm: public H2_3D_Tree {
public:
    kernel_Logarithm(double L, int tree_level, int interpolation_order, double epsilon, int use_chebyshev):H2_3D_Tree(L,tree_level,interpolation_order,epsilon, use_chebyshev){};
     void SetKernelProperty();
     double EvaluateKernel(vector3& targetpos, vector3& sourcepos);
};

/*! OneOverR2 kernel */
class kernel_OneOverR2: public H2_3D_Tree {
public:
    kernel_OneOverR2(double L, int tree_level, int interpolation_order, double epsilon, int use_chebyshev):H2_3D_Tree(L,tree_level,interpolation_order,epsilon, use_chebyshev){};
     void SetKernelProperty();
     double EvaluateKernel(vector3& targetpos, vector3& sourcepos);
};

/*! Quadric kernel */
class kernel_Quadric: public H2_3D_Tree {
public:
    kernel_Quadric(double L, int tree_level, int interpolation_order, double epsilon, int use_chebyshev):H2_3D_Tree(L,tree_level,interpolation_order,epsilon, use_chebyshev){};
     void SetKernelProperty();
     double EvaluateKernel(vector3& targetpos, vector3& sourcepos);
};

/*! InverseQuadric kernel */
class kernel_InverseQuadric: public H2_3D_Tree {
public:
    kernel_InverseQuadric(double L, int tree_level, int interpolation_order, double epsilon, int use_chebyshev):H2_3D_Tree(L,tree_level,interpolation_order,epsilon, use_chebyshev){};
     void SetKernelProperty();
     double EvaluateKernel(vector3& targetpos, vector3& sourcepos);
};

/*! ThinPlateSpline kernel */
class kernel_ThinPlateSpline: public H2_3D_Tree {
public:
    kernel_ThinPlateSpline(double L, int tree_level, int interpolation_order, double epsilon, int use_chebyshev):H2_3D_Tree(L,tree_level,interpolation_order,epsilon, use_chebyshev){};
     void SetKernelProperty();
     double EvaluateKernel(vector3& targetpos, vector3& sourcepos);
};

/*! Stokes kernel */
class kernel_Stokes: public H2_3D_Tree {
public:
    kernel_Stokes(double L, int tree_level, int interpolation_order, double epsilon, int use_chebyshev):H2_3D_Tree(L,tree_level,interpolation_order, epsilon, use_chebyshev){};
     void SetKernelProperty();
     double EvaluateKernel(vector3& targetpos, vector3& sourcepos);
};



#endif //(__kerbel_Types_hpp__)
