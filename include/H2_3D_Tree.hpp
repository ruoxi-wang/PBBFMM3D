
/*! \file H2_3D_Tree.hpp
 The FMM tree
 */
#ifndef __H2_3D_Tree_hpp__
#define __H2_3D_Tree_hpp__

#include"bbfmm.h"
#include"environment.hpp"
#include "rfftw.h"

using namespace std;
/*! The fmm tree */
class H2_3D_Tree{
public:
    /*! Constructor of class H2_3D_Tree */
    H2_3D_Tree( double L,  int level, int n, double
               epsilon, int use_chebyshev);
    doft* dof;
    double L;
    int level;
    int n;
    double epsilon;
	nodeT *tree;    // Octree for FMM hierarchy
	double homogen; // Order of kernel homogeneity
    int symmetry;
    double alpha;
    int use_chebyshev;
    
	doft cutoff;    // Compression index of SVD
	char Kmat[50];  // File name for kernel interaction matrix K
	char Umat[50];  // File name for matrix of singular vectors U
	char Vmat[50];  // File name for matrix of singular vectors V
    string kernelType;
	
    double Ucomp, Vcomp;
	
	int n2;
	int n3;
	int dofn3_s;
	int dofn3_f;
	
    rfftw_plan p_r2c;
    rfftw_plan p_c2r;
    
    bool computed;
	double *Kweights, *Cweights, *Tkz;
    int skipLevel;
    
	// All the cells in consideration: 7^3 = 343
	int Ktable[343];
	
    //void get_Compression_Rate(double* Ucomp, double* Vcomp);//Compression rate of the singular vector matrices
	
    // Set up for FMM
	void FMMSetup(nodeT **A, double *Tkz, double *Kweights,
                  double *Cweights, double L, doft *cutoff,
                  int n,  double epsilon, doft * dof,  int level, char *Kmat, char *Umat, char *Vmat,
                  double *Ucomp, double *Vcomp, int& skipLevel, double alpha, int use_chebyshev,rfftw_plan p_r2c);
    
    void buildFMMTree();
    void ComputeKernelSVD(double *Kweights, int n,double epsilon, doft *dof, char*Kmat, char *Umat, char *Vmat, int symm,  double *Ucomp,double *Vcomp, double alphaAdjust, double boxLen);
    void ComputeKernelUniformGrid(double *Kweights, int n, doft *dof,  char *Kmat, double alphaAdjust, rfftw_plan p_r2c);
    
    void ComputeWeights(double *Tkz, int *Ktable, double *Kweights,
                        double *Cweights, int n,double alpha, int use_chebyshev);
    
    void CalculateNodeLocations(int n, double* nodes, int use_chebyshev);
    
    void ComputeSn(vector3 *point, double *Tkz, int n, int N, vector3 *Sn, int use_chebyshev);
    
    void ComputeTk(double x, int n, double *vec);
    
    void EvaluateKernelCell(vector3 *field, vector3 *source, int Nf,
                            int Ns, doft *dof, double *kernel);
	
    // Initialize arrays, K is the compressed M2L operator C^{(i)}
	double *K, *U, *VT;
	
	/* U: U^k_r p. 8719; downward pass; field
	 * V: S^K_r p. 8718; upward pass; source
	 */
	
    // Read kernel interaction matrix K and singular vectors U and VT
	void FMMReadMatrices(double *K, double *U, double *VT, doft *cutoff, int n, doft *dof,char *Kmat, char *Umat, char *Vmat, int treeLevel, double homogen, int skipLevel, int use_chebyshev);
    
    void BuildFMMHierarchy(nodeT **A, int level, int n, doft *cutoff, doft *dof);
    
    void NewNode(nodeT **A, vector3 center, double L, int n);
    
    void FreeNode(nodeT *A);
        
    void GetPosition(int n, int idx, double *fieldpos, double *sourcepos, double *nodepos);
    
    virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                                 double *K, doft *dof){};
    void EvaluateKernelMulti(vector3 fieldpos, vector3 sourcepos,
                        double *K, doft *dof, int m);
    virtual void setHomogen(string& kernelType, doft* dof){};

    /*! Destructor of class H2_3D_Tree */
    ~H2_3D_Tree();
	

};


#endif //(__H2_3D_Tree_hpp__)
