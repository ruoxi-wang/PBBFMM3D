/*
 * File: bbfmm.h
 * Description: Header file for bbfmm.c which contains functions used
 * in the implementation of the black-box fast multiple method.
 * ----------------------------------------------------------------------
 * 
 * Black-Box Fast Multipole Method (BBFMM)
 * William Fong
 * Stanford University
 *
 */

#ifndef _BBFMM_H
#define _BBFMM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include<sys/time.h>

// 3d point
typedef struct _vector3 {
	double x,y,z;
} vector3;

// 3d segment
typedef struct _point2 {
	vector3 p1;
	vector3 p2;
} point2;

typedef double timeType;

/* Uniform random number generator */
#define frand(xmin,xmax) ((double)xmin+(double)(xmax-xmin)*rand()/ \
(double)RAND_MAX) 

/*
 * Function: Timer
 * ----------------------------------------
 * Returns the time in seconds.
 *
 */
extern timeType Timer(void);
#define evaltime(timeval_time) (double)timeval_time.tv_sec \
+ (double)timeval_time.tv_usec*1e-6

/* Struct: nodeT
 * -------------------------------------------------------------------
 * This struct is a node in the octree.
 */
typedef struct _nodeT {
	struct _nodeT *leaves[8];
	struct _nodeT *parent;
	struct _nodeT *neighbors[27];
	struct _nodeT *interaction[189];
	vector3 center, cshiftneigh[27], cshiftinter[189];
	double length;
	double *fieldval, *sourceval, *proxysval,*sourcefre;
	int *fieldlist, *sourcelist;
	int Nf, Ns, ineigh, iinter;
} nodeT;


/* Struct: dof
 * -------------------------------------------------------------------
 * The FMM can handle cases where the input are output variables
 * are vectors (not just scalars).
 * This struct stores information about the size of the input
 * vector (e.g., 1, 3, 9, etc) and the output vector (e.g., 1, 3, 6).
 * This is different from the number of source and field points.
 */
typedef struct _dof_struct {
	int s; // Size of source vector
	int f; // Size of field vector
} doft; // Type of variable dof


/* Struct: fmmparam
 * -------------------------------------------------------------------
 * This struct stores all of the FMM parameters.
 */
typedef struct _fmmparam {
	int Ns;       // Number of sources
	int Nf;       // Number of field points
	doft dof;        // Number of degrees of freedom
	double L;        // Length of one side of simulation cube
	int n;        // Number of Chebyshev nodes in one direction
	int levels;        // Maximum number of levels in octree
	int PBClevels;        // Number of PBC levels (images = 27^PBClevels)
	int PBCshells;        // Number of PBC shells (for direct calculation)
	int precomp;        // Turn on (1) or off (0) pre-computation
	doft cutoff;       // Number of singular values to keep
	double homogen;        // Order of homogeneity of kernel
	char filesval[80];     // File name for storing singular values
	char filesvec[80];     // File name for storing singular vectors
} fmmparam;

// Declaration for LAPACK's LU decomposition routine
extern void dgetrf_(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);

// Declaration for LAPACK's matrix inverse using LU decomposition routine
extern void dgetri_(int *N, double *A, int *lda, int *IPIV, double *WORK, 
        int *lwork, int *INFO);

extern void MatrixInverse(double *A, int N);

        // Declaration for LAPACK's eigen routine
extern void dgeev_(char *jobvl, char *jobvr, int *n, double *A, int *lda,
        double *wr, double *wi, double *vl, int *ldvl, double *vr,
        int *ldvr, double *work, int *lwork, int *info);

	// Declaration for LAPACK's SVD routine
extern void dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *A,
					int *lda, double *S, double *U, int *ldu, double *VT,
					int *ldvt, double *work, int *lwork, int *info);

	// Declaration for BLAS matrix-matrix multiply
extern void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
				   double *alpha, double *A, int *lda, double *B,
				   int *ldb, double *beta, double *C, int *ldc);

	// Declaration for BLAS matrix-vector multiply
extern void dgemv_(char *trans, int *m, int *n, double *alpha, double *A,
				   int *lda, double *x, int *incx, double *beta, 
				   double *y, int *incy);

	// Declaration for BLAS dot product
extern double ddot_(int *n, double *dx, int *incx, double *dy, int *incy);

	// Declaration for BLAS daxpy
extern double daxpy_(int *n, double *da, double *dx, int *incx, double *dy,
					 int *incy);

	// Declaration for BLAS vector copy
extern double dcopy_(int *n, double *dx, int *incx, double *dy, int *incy);

    
#ifdef __cplusplus
}
#endif

#endif
