
/*! \file	H2_3D_Tree.cpp
*/

#include"H2_3D_Tree.hpp"
#include"bbfmm.h"

#define HOMO_THRESHOLD  1e-1          // threshold for homogeneous kenrel
#define FFTW_FLAG       FFTW_ESTIMATE // option for fftw plan type


H2_3D_Tree::H2_3D_Tree(double L, int tree_level, int interpolation_order,  double epsilon, int use_chebyshev){
    this->dof = new doft;
    this->L     =   L;
    this->tree_level =   tree_level;
    this->interpolation_order    =   interpolation_order;
    this->epsilon   =   epsilon;
    this->use_chebyshev = use_chebyshev;
    alpha = 0;
    n2 = interpolation_order*interpolation_order;        
	  n3 = n2*interpolation_order;      
    this->dof->s = 1;
    this->dof->f = 1;
    this->indexToLeafPointer = (nodeT**) malloc(pow(8, this->tree_level) * sizeof(nodeT*));

    // Omega matrix
	Kweights = (double *)malloc(n3 * sizeof(double));    
	// Chebyshev interpolation coefficients: Sn (page 8715)
	Cweights = (double *)malloc(2 * n2 * sizeof(double));

    Tkz      = (double *)malloc(n2 * sizeof(double));
    K       =   NULL;
    U       =   NULL;
    VT      =   NULL;
    tree    =   NULL;
    skipLevel = 0;
    computed   =   false;
}

void H2_3D_Tree::buildFMMTree() {


    vector3 center;
    SetKernelProperty();
    homogen = -homogen;
    // Compute the Chebyshev weights and sets up the lookup table

    ComputeWeights(Tkz,Ktable,Kweights,Cweights,interpolation_order,alpha,use_chebyshev); 

    // Precompute the SVD of the kernel interaction matrix (if
    // necessary)

    PrecomputeM2L(Kweights, L, alpha, &cutoff, interpolation_order, homogen, epsilon, *dof,
        tree_level, &K, &U, &VT, use_chebyshev);

        // Builds the FMM hierarchy
    center.x = 0;
    center.y = 0;
    center.z = 0;

    NewNode(&tree,center,L,interpolation_order);

    BuildFMMHierarchy(&tree,tree_level,interpolation_order,&cutoff,dof, 0, indexToLeafPointer, cellPointers);


}

/*
 * Function: ComputeM2L
 * -----------------------------------------------------------------
 * Prepare for the FMM calculation by pre-computing the SVD (if necessary),
 * and reading in the necessary matrices.
 */
void H2_3D_Tree::PrecomputeM2L(double *Kweights, double boxLen, double alpha,
        doft *cutoff, int interpolation_order, int homogen,
        double epsilon, doft dof, int treeLevel,
        double **K, double **U, double **VT, int use_chebyshev) {

  char Kmat[50];
  char Umat[50];
  char Vmat[50];

  if (use_chebyshev) {
    sprintf(Kmat,"./../output/%sChebK%d.bin", kernelType.c_str(), interpolation_order);
    sprintf(Umat,"./../output/%sChebU%d.bin", kernelType.c_str(), interpolation_order);
    sprintf(Vmat,"./../output/%sChebV%d.bin", kernelType.c_str(), interpolation_order);
  } else {
    sprintf(Kmat,"./../output/%sUnifK%d.bin", kernelType.c_str(), interpolation_order);
    sprintf(Umat,"bbfmm.c"); // uniform grid does not have U or V file,
    sprintf(Vmat,"bbfmm.c"); // so just make sure these two files exist
  }
    
  if ( !PrecomputeAvailable(Kmat, Umat, Vmat, homogen, boxLen, treeLevel, use_chebyshev) ) {
    StartPrecompute( boxLen, treeLevel, interpolation_order, dof, homogen, symmetry, Kmat, Umat, Vmat, alpha, Kweights, epsilon, use_chebyshev );
  }
  // Read kernel interaction matrix K and singular vectors U and VT
  FMMReadMatrices(K,U,VT, cutoff,interpolation_order,dof,Kmat,Umat,Vmat, treeLevel, homogen, use_chebyshev);

}

// check if the precompute files exist or if the existing files are usable
bool H2_3D_Tree::PrecomputeAvailable( char *Kmat, char *Umat, char *Vmat,
              double homogen, double boxLen,
              int treeLevel, int grid_type ) {

  bool avail = true;
  
  FILE *fK, *fU, *fV;
  fK = fopen(Kmat, "rb");
  fU = fopen(Umat, "rb");
  fV = fopen(Vmat, "rb");
     
  if (fK == NULL || fU == NULL || fV == NULL) { // files do not exist    
    avail = false;
  } else if ( !IsHomoKernel(homogen) ) { // non-homogeneous kernel
    double len_file;
    int    lvl_file;    
    READ_CHECK( fread(&len_file, sizeof(double), 1, fK), 1 );
    READ_CHECK( fread(&lvl_file, sizeof(int),    1, fK), 1 );    
    if ( lvl_file != treeLevel || (fabs(len_file - boxLen) >= EQUAL_DOUBLE_EPS) )  
      avail = false;
  }

  if (fK!=NULL) fclose(fK);
  if (fU!=NULL) fclose(fU);
  if (fV!=NULL) fclose(fV);
  if (avail) printf("Precompute files exist.\n");
  else printf("Precompute files do NOT exist.\n");
  return avail;
}

bool H2_3D_Tree::IsHomoKernel( double homogen ) {
  return fabs(homogen) > HOMO_THRESHOLD;
}

void H2_3D_Tree::StartPrecompute(double boxLen, int treeLevel, int interpolation_order, doft dof, int homogen, int symmetry, char *Kmat, char *Umat, char *Vmat, double alpha, double *Kweights, double epsilon, int use_chebyshev) {

  printf("Generate precompute file ...\n");
  FILE *fK = fopen(Kmat, "wb");
  fwrite(&boxLen, sizeof(double), 1, fK); // write box size
  fwrite(&treeLevel, sizeof(int), 1, fK); // write tree level
  fclose(fK);                
  bool first_time_call = true;

  if ( IsHomoKernel(homogen) ) {  // homogeneours kernel

    double unit_len = 1.0;
    compute_m2l_operator(interpolation_order, dof, symmetry, Kmat, Umat, Vmat, unit_len, alpha, Kweights, epsilon, use_chebyshev, first_time_call);
    
  } else {                                // non-homogeneous kernel
    int j;
    double boxLenLevel = boxLen/4;        // FMM starts from the second level
    for (j=2; j<=treeLevel; j++) {
      compute_m2l_operator(interpolation_order, dof, symmetry, Kmat, Umat, Vmat, boxLenLevel, alpha, Kweights, epsilon, use_chebyshev, first_time_call);
      boxLenLevel/=2;
      first_time_call = false;
    }

  }                                       // end non-homogeneous kernel
}


void H2_3D_Tree::compute_m2l_operator (int interpolation_order, doft dof, int symmetry, char *Kmat, char *Umat, char *Vmat, double l, double alpha, double *Kweights, double epsilon, int grid_type, bool first_time_call) {

  switch (use_chebyshev) {
  case 0:
    ComputeKernelUnif(interpolation_order, dof, Kmat, alpha, l);
    break;
  case 1:
    ComputeKernelCheb(Kweights, interpolation_order, epsilon, &dof,
                     Kmat, Umat, Vmat, symmetry, alpha,
                     l, first_time_call);
    break;
  }
}


/*
 * Function: FMMReadMatrices
 * ------------------------------------------------------------------
 * Read in the kernel interaction matrix M and the matrix of singular
 * vectors U.
 */
void H2_3D_Tree::FMMReadMatrices(double **K, double **U, double **VT, doft *cutoff,
             int interpolation_order, doft dof, char *Kmat, char *Umat,
             char *Vmat, int treeLevel, double homogen,
             int use_chebyshev) {

  int preCompLevel;
  if ( !IsHomoKernel(homogen) )
    preCompLevel = treeLevel - 1;
  else
    preCompLevel = 1;
  

  unsigned long long Ksize;
  unsigned int Usize;
  unsigned int Vsize;
  
  if (!use_chebyshev) {
       
    cutoff->f = interpolation_order*interpolation_order*interpolation_order;
    cutoff->s = interpolation_order*interpolation_order*interpolation_order;

    Ksize = 316*(2*interpolation_order-1)*(2*interpolation_order-1)*(2*interpolation_order-1)*dof.s*dof.f * preCompLevel;
    
  } else if (use_chebyshev) { // read 'U' and 'V'
    
    FILE *fU = fopen(Umat,"rb"); 
    FILE *fV = fopen(Vmat,"rb");

    READ_CHECK( fread(&(cutoff->f), sizeof(int), 1, fU), 1 );
    READ_CHECK( fread(&(cutoff->s), sizeof(int), 1, fV), 1 );

    Usize = cutoff->f * interpolation_order*interpolation_order*interpolation_order * preCompLevel;
    Vsize = cutoff->s * interpolation_order*interpolation_order*interpolation_order * preCompLevel;
    Ksize = 316 * cutoff->f*cutoff->s * preCompLevel;
    
    *U  = (double *)malloc(Usize *sizeof(double)); 
    *VT = (double *)malloc(Vsize *sizeof(double));
    READ_CHECK( fread(*U,  sizeof(double), Usize, fU), Usize );
    READ_CHECK( fread(*VT, sizeof(double), Vsize, fV), Vsize );

    fclose(fU);
    fclose(fV);
    
  } else
    assert(false);

  
  FILE *fK  = fopen(Kmat, "rb");
  assert(Ksize > 0); // check over flow
  *K = (double*)malloc(Ksize *sizeof(double));    
  fseek(fK, 1*sizeof(int) + 1*sizeof(double), SEEK_SET); // skip 'tree level' and 'box length'
  READ_CHECK( fread(*K, sizeof(double), Ksize, fK), Ksize );
  fclose(fK);
}

/*
 * Function: ComputeKernelCheb
 * ---------------------------------------------------------------------
 * Computes the kernel for 316n^6 interactions between Chebyshev nodes
 * and then computes the SVD of the kernel matrix.
 * symmetry = 0 no symmetric property, 1 symmetric kernel or -1 anti-
 * symmetric kernel
 */


void H2_3D_Tree::ComputeKernelCheb(double *Kweights, int interpolation_order,double epsilon, doft *dof, char*Kmat, char *Umat, char *Vmat, int symm, double alphaAdjust, double boxLen,  bool first_time_call) {
    
    int i, j, l, m, k1, k2, k3, l1, l2, l3, z;
    int count, count1;
    vector3 vtmp;
    double pi = M_PI;
	
    int interpolation_order3 = interpolation_order*interpolation_order*interpolation_order;            // n3 = n^3
    int dofn3_s = n3;
    int dofn3_f = n3;
    size_t dof2n6 = dofn3_s * dofn3_f; // Total size
    int Sigma_size;
    doft cutoff;
	
    double *K0, *U0, *Sigma, *VT0;
    double *nodes, *work;
    vector3 *targetpos;
	
    K0 = (double *)malloc(316 * dof2n6 * sizeof(double));
    targetpos  = (vector3 *)malloc(n3 * sizeof(vector3));
	
    // 316 M2L operators
    double* Kcell[316];
    for (z = 0; z < 316; ++z)
        Kcell[z] = (double *) malloc(dof2n6 * sizeof(double));
	
	
    // Compute Chebyshev nodes of T_n(x)
    double scale = 1+alphaAdjust;
    nodes = (double *)malloc(interpolation_order * sizeof(double));
    for (m=0;m<interpolation_order;m++)
        nodes[m] = cos(pi*((double)m+0.5)/(double)interpolation_order) * scale;
	
    // Compute the locations of the target points in a unit cube
    count = 0;
    for (l1=0;l1<interpolation_order;l1++) {
        vtmp.x = 0.5 * nodes[l1] * boxLen;
        for (l2=0;l2<interpolation_order;l2++) {
            vtmp.y = 0.5 * nodes[l2] * boxLen;
            for (l3=0;l3<interpolation_order;l3++) {
                targetpos[count].x = vtmp.x;
                targetpos[count].y = vtmp.y;
                targetpos[count].z = 0.5 * nodes[l3] * boxLen;
                count++;
            }
        }
    }
    
    // Compute the kernel values for interactions with all 316 cells
    int countM2L=0;
    int symmNum = 158*(2-abs(symm)); // symmNum = 158, 316 for symm=pm1, 0

    int col, row, idx1, idx2;
    #pragma omp parallel for private(k1, k2, k3, countM2L, l1, l2, l3) collapse(3)
    for (k1=-3;k1<4;k1++) {
        for (k2=-3;k2<4;k2++) {
            for (k3=-3;k3<4;k3++) {
                if (abs(k1) > 1 || abs(k2) > 1 || abs(k3) > 1) {
                    countM2L = (k1+3)*49+(k2+3)* 7 + (k3+3) - (min(max(k1 + 1, 0), 3) *9 +  
                        (-1 <= k1 && k1 <= 1 ? 1 : 0)* min(max(k2+1, 0), 3)*3 + 
                        (-1 <= k1 && k1 <= 1 ? 1 : 0)*(-1 <= k2 && k2 <= 1 ? 1 : 0) * min(max(k3 + 1, 0), 3));
                    if (countM2L < symmNum) {
                        vector3 * sourcepos = (vector3 *)malloc(n3 * sizeof(vector3));
                        double*   kernel = (double *)malloc(dof2n6 * sizeof(double));

                        // Compute the locations of the source points in the cell
                        for (l1=0;l1<interpolation_order;l1++) {
                            for (l2=0;l2<interpolation_order;l2++) {
                                for (l3=0; l3<interpolation_order; l3++) {
                                    int ind = l3 + l2*interpolation_order+l1*interpolation_order*interpolation_order;
                                    sourcepos[ind].x = (k1 + 0.5 * nodes[l1]) * boxLen;
                                    sourcepos[ind].y = (k2 + 0.5 * nodes[l2]) * boxLen;
                                    sourcepos[ind].z = (k3 + 0.5 * nodes[l3]) * boxLen;
                                }
                            }
                        }
                        
                        // Compute the kernel at each of the target Chebyshev nodes
                        EvaluateKernelCell(targetpos, sourcepos, n3, n3, dof, kernel);

                        
                        // Copy the kernel values to the appropriate location
                        int globle_offset = countM2L * n3*n3;
                        int local_offset;
                        for (l1=0; l1<n3; l1++) {
                            for (l2=0; l2<n3; l2++) {
                                local_offset = l1*n3+l2;
                                K0[globle_offset + local_offset] = kernel[local_offset]/( Kweights[l1]*Kweights[l2]);
                            }
                        }
                        free(sourcepos);
                        free(kernel);
                    }
                }
            }
        }
    }

    if (symm) {
        countM2L=0;
    	for (k1=-3;k1<4;k1++)
            for (k2=-3;k2<4;k2++)
                for (k3=-3;k3<4;k3++) 
                    if (abs(k1) > 1 || abs(k2) > 1 || abs(k3) > 1) {
                        if (countM2L >= symmNum) {
                            for (col=0; col<dofn3_s; col++)
                                for (row=0; row<dofn3_f; row++) {
                                    idx1 = (315-countM2L)*dof2n6 + col*dofn3_f + row;
                                    idx2 = countM2L*dof2n6 + row *dofn3_f + col;
                                    // symm=-1,1 for anti-symmetry and symmetry
                                    K0[idx2] = symm * K0[idx1];
                                }
                        }
                        countM2L++;
                    }
    }

    // Extract the submatrix for each of the 316 cells
    count = 0;
    for (i=0;i<316;i++) {
        for (j=0;j<(int)dof2n6;j++) {
            Kcell[i][j] = K0[count];
            count++;
        }
    }
	
    /****
     * Compute the SVD of K_fat
     ****/
	
    // Compute the SVD of K0
    char save[]="S", nosave[]="N";
    int nosavedim=1;
    int info, lwork;
    int cols_s = 316*dofn3_s;
	
    /* See dgesvd documentation:
     *          LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code):
     *             - PATH 1  (M much larger than N, JOBU='N')  - our case for K_thin
     *             - PATH 1t (N much larger than M, JOBVT='N') - our case for K_fat
     */
    lwork = 5*n3; // Change
    work  = (double *) malloc(lwork * sizeof(double));
	
    int U0size;
    Sigma_size = dofn3_f;
    U0size     = dofn3_f * dofn3_f;
    Sigma      = (double *)malloc(Sigma_size * sizeof(double));
    U0         = (double *)malloc(U0size * sizeof(double));
    VT0        = NULL;
    dgesvd_(save, nosave, &dofn3_f, &cols_s, K0, &dofn3_f, Sigma, U0,
            &dofn3_f, VT0, &nosavedim, work, &lwork, &info);

    FILE *ptr_file;

    double sumsqr, sum, epsqr;
    if (first_time_call) { // The first time called
        
        // Determine the number of singular values to keep. We use epsilon for this.
        sumsqr = 0;
        for (i=Sigma_size-1; i>=0; --i)
            sumsqr += Sigma[i] * Sigma[i];
        
        sum = 0, epsqr = sumsqr * epsilon * epsilon;
        cutoff.f = Sigma_size;
        for (i=Sigma_size-1; i>=0; --i) {
            sum += Sigma[i] * Sigma[i];
            if (sum < epsqr)
                --cutoff.f;
            else
                break;
        }
        
        // Extract the needed columns from U0 and write out to a file
        ptr_file = fopen(Umat,"wb");

        fwrite(&cutoff.f, sizeof(int), 1, ptr_file);
        fclose(ptr_file);
    }
    else {
        ptr_file = fopen(Umat, "rb");
        if (ptr_file == NULL)
            printf("Can't read cutoff.f!\n");
        i = fread(&cutoff.f, sizeof(int), 1, ptr_file);
        if (i != 1)
            printf("fread error in ComputeKernelCheb().\n");
        fclose(ptr_file);
    }
        
    double trancatedSize = dofn3_f * cutoff.f;
    ptr_file = fopen(Umat, "ab");
    fwrite(U0, sizeof(double), trancatedSize, ptr_file);
    fclose(ptr_file);
	

    free(Sigma); Sigma = NULL;
	
	
    /****
     * Compute the SVD of K_thin
     ****/
	
    // Form K_thin using all of 316 M2L operators stored in Kcell
    count = 0;
    for (j=0;j<dofn3_s;++j) {
        for (i=0;i<316;i++) {
            for (l=0;l<dofn3_f;l++) {
                K0[count] = Kcell[i][l+j*dofn3_f];
                count++;
            }
        }
    }
	
    // save = "S"; nosave = "N"
    double *U1 = NULL;
    int rows_f = 316*dofn3_f;
    Sigma_size = dofn3_s;
    Sigma      = (double *)malloc(Sigma_size * sizeof(double));
    VT0        = (double *)malloc(dofn3_s * dofn3_s * sizeof(double));
    
    // Compute the SVD of the K_thin
    dgesvd_(nosave,save,&rows_f,&dofn3_s,K0,&rows_f,Sigma,U1,&nosavedim,VT0,&dofn3_s,
            work,&lwork,&info);
    
    if (first_time_call) {
        
        // Determine the number of singular values to keep. We use
        // epsilon for this.
        
        sumsqr = 0;
        for (i=Sigma_size-1; i>=0; --i)
            sumsqr += Sigma[i] * Sigma[i];
        
        sum = 0, epsqr = sumsqr * epsilon * epsilon;
        cutoff.s = Sigma_size;
        for (i=Sigma_size-1; i>=0; --i) {
            sum += Sigma[i] * Sigma[i];
            if (sum < epsqr)
                --cutoff.s;
            else
                break;
        }
        
        // Extract trancated VT[cutoff.s][dofn3_s] from
        // VT0[dofn3_s][dofn3_s] and write out to a file
        
        ptr_file = fopen(Vmat,"wb");
        fwrite(&cutoff.s, sizeof(int), 1, ptr_file);
        fclose(ptr_file);
    }
    else {
        ptr_file = fopen(Vmat, "rb");
        if (ptr_file == NULL)
            printf("Can't read cutoff.s!\n");
        i = fread(&cutoff.s, sizeof(int), 1, ptr_file);
        if (i != 1)
            printf("fread error in ComputeKernelCheb().\n");
        fclose(ptr_file);
    }
        
    ptr_file = fopen(Vmat, "ab");	    
    for (j=0;j<dofn3_s;j++) { // column
        count1 = j*dofn3_s;
        fwrite(VT0+count1, sizeof(double), cutoff.s, ptr_file);
    }
    fclose(ptr_file);
    
    
    /** Computing the compressed kernel using the orthonormal basis U and VT.
     **/
    // char *transa, *transb;
    // transa = new char[2];
    // transb = new char[2];
    int cutoff2 = cutoff.f * cutoff.s;	
    double alpha=1, beta=0;
    double *Ecell;
    Ecell = (double *)malloc(cutoff2 *316* sizeof(double));	
	
    ptr_file = fopen(Kmat,"ab");
    #pragma omp parallel for private(i)
    for (i=0;i<316;i++) {       
        /* Compute K V:
         * K  is dofn3_f  x dofn3_s
         * VT is cutoff.s x dofn3_s
         * V  is dofn3_s  x cutoff.s 
         * KV is dofn3_f  x cutoff.s
         * (Notice that VT is a submatrix of VT0)
         */
        double* KV    = (double *)malloc(dofn3_f * cutoff.s * sizeof(double));
        char transa = 'n', transb = 't';

        dgemm_(&transa, &transb, &dofn3_f, &cutoff.s, &dofn3_s, &alpha, 
               Kcell[i], &dofn3_f, VT0, &dofn3_s, &beta, KV, &dofn3_f);
		
        /* Compute U^T K V:
         * KV is dofn3_f x cutoff.s
         * U  is dofn3_f x cutoff.f
         * UT is cutoff.f x dofn3_f
         * U^T K V is cutoff.f x cutoff.s
         * (Notice that U is a submatrix of U0)
         */
        transa = 't';
        transb = 'n';
        dgemm_(&transa, &transb, &cutoff.f, &cutoff.s, &dofn3_f, &alpha, 
               U0, &dofn3_f, KV, &dofn3_f, &beta, Ecell+i*cutoff2, &cutoff.f);
        free(KV);	
    }

    fwrite(Ecell, sizeof(double), cutoff2*316, ptr_file);

    fclose(ptr_file);
    
    free(K0);
    free(VT0);
    free(U0);
    free(Sigma);
    free(nodes);
    free(targetpos);
    free(Ecell);
    free(work);

    for (z = 0; z < 316; ++z) {
        free(Kcell[z]);
    }
    
}



/*
 * Function: ComputeKernelUnif
 * ------------------------------------------------------------------
 * Computes the kernel for 316(2n-1)^3 interactions between Uniform
 * Grid nodes.
 */
void H2_3D_Tree::ComputeKernelUnif(int interpolation_order, doft dof, char *Kmat,
               double alphaAdjust, double len) {
  int i, k1, k2, k3, l1, l2, l3;    
  double nodes[interpolation_order];
  GridPos1d(alphaAdjust, len, interpolation_order, 0, nodes);

  int vecSize = 2*interpolation_order-1, reducedMatSize = pow(vecSize, 3);
  int M2LSize = reducedMatSize;
  double *freqMat = fftw_alloc_real( 316*M2LSize );
  double *MatM2L = fftw_alloc_real(316*M2LSize);
  // Compute the kernel values for interactions with all 316 cells

  int countM2L;
  fftw_plan p[316];
  for (i = 0; i < 316; i++) {
    p[i] = fftw_plan_r2r_1d(reducedMatSize, MatM2L+ i*M2LSize, freqMat + i*M2LSize , FFTW_R2HC, FFTW_FLAG);
  }

  #pragma omp parallel for private(k1, k2, k3, countM2L, l1, l2, l3) collapse(3)
  for (k1=-3;k1<4;k1++) {
    for (k2=-3;k2<4;k2++) {
      for (k3=-3;k3<4;k3++) {
        if ( abs(k1) > 1 || abs(k2) > 1 || abs(k3) > 1 ) {
            countM2L =  (k1+3)*49+(k2+3)* 7 + (k3+3) - (min(max(k1 + 1, 0), 3) *9 +  
                        (-1 <= k1 && k1 <= 1 ? 1 : 0)* min(max(k2+1, 0), 3)*3 + 
                        (-1 <= k1 && k1 <= 1 ? 1 : 0)*(-1 <= k2 && k2 <= 1 ? 1 : 0) * min(max(k3 + 1, 0), 3));
            vector3 scenter, targetpos, sourcepos;
            scenter.x = (double)k1*len;
            scenter.y = (double)k2*len;
            scenter.z = (double)k3*len;

            int count=0;
            for (l1=0; l1<vecSize; l1++) {
                GetPosition(interpolation_order, l1, &targetpos.x, &sourcepos.x, nodes);
                sourcepos.x += scenter.x;
                for (l2=0; l2<vecSize; l2++) {
                    GetPosition(interpolation_order, l2, &targetpos.y, &sourcepos.y, nodes);
                    sourcepos.y += scenter.y;
                    for (l3=0; l3<vecSize; l3++, count++) {
                        GetPosition(interpolation_order, l3, &targetpos.z, &sourcepos.z, nodes);
                        sourcepos.z += scenter.z;

                        MatM2L[countM2L*M2LSize+count] = EvaluateKernel(targetpos, sourcepos);

                    }
                }
            }       

           // FFT
            fftw_execute(p[countM2L]);
            fftw_destroy_plan( p[countM2L]);

        }
      }
    }
  }

  FILE *ptr_file;
  ptr_file = fopen(Kmat, "ab");
  fwrite(freqMat, sizeof(double), 316 *M2LSize, ptr_file);
  fclose(ptr_file);  
  fftw_free(freqMat);
  fftw_free(MatM2L);
}


/*
 * Function: ComputeGrids
 * ---------------------------------------------------------
 * Calculates node locations for grids between [1,-1] and
 * stores them in pre-allocated array nodes.
 */                          
void H2_3D_Tree::GridPos1d(double alpha, double len, int interpolation_order, int use_chebyshev,
           double* nodes) {
      
  int m;
  double L = len * (1+alpha);
  if (use_chebyshev) {
    double pi = M_PI;
    for (m=0; m<interpolation_order; m++)
      nodes[m] = cos(pi*((double)m+0.5)/(double)interpolation_order) * L;
  }
  else
    for (m=0; m<interpolation_order; m++)
      nodes[m] = 1 - 2*(double)m/((double)(interpolation_order-1)) * L;
}


/*
 * Given n node positions and returns corresponding target and source
 * positions with respect to the index
 */
void H2_3D_Tree::GetPosition(int interpolation_order, int idx, double *targetpos, double *sourcepos, double *nodepos) {
    
    if (idx < interpolation_order) {
        *targetpos  = nodepos[interpolation_order-1]/2;
        *sourcepos = nodepos[idx]/2;
    } else {
        *targetpos  = nodepos[2*(interpolation_order-1)-idx]/2;
        *sourcepos = nodepos[interpolation_order-1]/2;
    }
    
}


/*
 * Function: EvaluateKernelCell
 * -------------------------------------------------------------------
 * Evaluates the kernel for interactions between a pair of cells.
 * M2L operator initialization
 */
void H2_3D_Tree::EvaluateKernelCell(vector3 *target, vector3 *source, int Nf,
                        int Ns, doft *dof, double *kernel) {
	int i, j;	
	for (j=0;j<Ns;j++) {
        vector3 cur_source = source[j];
		for (i=0;i<Nf;i++) {
            kernel[i + Nf * j] = EvaluateKernel(target[i],cur_source);
		}
	}
}

/*
 * Function: ComputeWeights
 * ------------------------------------------------------------------
 * Computes the weights for the Chebyshev nodes of all children cells
 * (identical for all cells and all levels so just compute once and
 * store in memory) and set up the lookup table.
 */

void H2_3D_Tree::ComputeWeights(double *Tkz, int *Ktable, double *Kweights,
					double *Cweights, int interpolation_order, double alpha, int use_chebyshev) {
	int i, k, m, l1, l2, l3, count, ncell, ninteract;
	double tmp1, tmp2;
	vector3 vtmp;
	
	double *nodes, *vec;
	vector3 *targett, *Sn;
	nodes = (double *)malloc(interpolation_order * sizeof(double));
	vec   = (double *)malloc(interpolation_order * sizeof(double));
	int interpolation_order3 = interpolation_order*interpolation_order*interpolation_order;                     // n3 = n^3
	int Nc = 2*n3;                      // Number of child Chebyshev nodes
	targett = (vector3 *)malloc(Nc * sizeof(vector3));
    // Chebyshev-transformed coordinates
	Sn     = (vector3 *)malloc(interpolation_order * Nc * sizeof(vector3));
	//double pi = M_PI;
	
    // Initialize lookup table
	for (i=0;i<343;i++)
		Ktable[i] = -1;
    // Create lookup table
	ncell = 0;
	ninteract = 0;
	for (l1=-3;l1<4;l1++) {
		for (l2=-3;l2<4;l2++) {
			for (l3=-3;l3<4;l3++) {
				if (abs(l1) > 1 || abs(l2) > 1 || abs(l3) > 1) {
					Ktable[ncell] = ninteract;
					ninteract++;
				}
				ncell++;
			}
		}
	}
    GridPos1d(0, 1.0, interpolation_order, use_chebyshev, nodes);

    // Evaluate the Chebyshev polynomials of degree 0 to n-1 at the nodes
	for (m=0; m<interpolation_order; m++) {
		ComputeTk(nodes[m],interpolation_order,vec);
		i = m*interpolation_order;
		for (k=0;k<interpolation_order;k++) {
			Tkz[i+k] = vec[k];
        }
	}

    if (use_chebyshev) {
		// Compute the weights for the kernel matrix K
        count = 0;
        for (l1=0;l1<interpolation_order;l1++) {
            tmp1 = 1./sqrt(1-nodes[l1]*nodes[l1]);
            for (l2=0;l2<interpolation_order;l2++) {
                tmp2 = tmp1/sqrt(1-nodes[l2]*nodes[l2]);
                for(l3=0;l3<interpolation_order;l3++) {
                    Kweights[count] = tmp2/sqrt(1-nodes[l3]*nodes[l3]);
                    count++;
                }
            }
        }
	}

	// Map all Chebyshev nodes from the children cells to the parent
	k = 0;
    double scale = 1 / (1+alpha);
    //scale = 1;
	for (i=0;i<2;i++) {
        // Determine the mapping function for the specific child cell
		vtmp.x = -0.5;
		vtmp.y = -0.5;
		
		if (i == 0)
			vtmp.z = -0.5;
		else
			vtmp.z = 0.5;
		
		for (l1=0;l1<interpolation_order;l1++) {
			for (l2=0;l2<interpolation_order;l2++) {
				for (l3=0;l3<interpolation_order;l3++) {
                    targett[k].x = nodes[l1] * .5 + vtmp.x * scale;
                    targett[k].y = nodes[l2] * .5 + vtmp.y * scale;
                    targett[k].z = nodes[l3] * .5 + vtmp.z * scale;
                    k++;
				}
			}
		}
	}

    
    // Compute Sc, the mapping function for the target points
	ComputeSn(targett,Tkz,interpolation_order,Nc,Sn, use_chebyshev);

    // Extract out the Chebyshev weights
	count = 0;
	for (l1=0;l1<interpolation_order;l1++) {
		k = l1*Nc;
		for (l2=0;l2<interpolation_order;l2++) {
			Cweights[count] = Sn[k+l2].z;
			count++;
		}
	}
	for (l1=0;l1<interpolation_order;l1++) {
		k = l1*Nc;
		for (l2=0;l2<interpolation_order;l2++) {
			Cweights[count] = Sn[k+n3+l2].z;
			count++;
		}
	}

	free(nodes);
	free(vec);
	free(targett);
	free(Sn);

}

double ClenshawSum(int interpolation_order, double x, double *Tk);
double LagrangeWeight(int interpolation_order, double x, int m);
void H2_3D_Tree::ComputeSn(vector3 *point, double *Tkz, int interpolation_order, int N, vector3 *Sn, int use_chebyshev){
     
  int i, k, m;

  if (use_chebyshev) {
        
    double pfac = 2./interpolation_order;     // TODO: figure out why 
    // double pfac = 1;     

    double *Tkz_m;

    // Loop over Chebyshev node m
    for (m=0;m<interpolation_order;m++) {
      k    = m*N;
      Tkz_m = Tkz + m*interpolation_order;
      for (i=0;i<N;i++) {
    // Compute S_n for each direction using Clenshaw
    Sn[k+i].x = pfac * ClenshawSum(interpolation_order, point[i].x, Tkz_m);
    Sn[k+i].y = pfac * ClenshawSum(interpolation_order, point[i].y, Tkz_m);
    Sn[k+i].z = pfac * ClenshawSum(interpolation_order, point[i].z, Tkz_m);
      }
    }

        
  } else {
     
    for (m=0;m<interpolation_order;m++) {
      k = m*N;
      for (i=0;i<N;i++) {
    Sn[k+i].x = LagrangeWeight(interpolation_order, point[i].x, m);
    Sn[k+i].y = LagrangeWeight(interpolation_order, point[i].y, m);
    Sn[k+i].z = LagrangeWeight(interpolation_order, point[i].z, m);
      }
    }
        
  } // end else
} // end function

// summation using Clenshaw's recurrence relation
double ClenshawSum(int interpolation_order, double x, double *Tk) {
  int j;
  double d0, d1, d2;
  d0 = d1 = 0;
  for (j = interpolation_order - 1; j > 0; j--) {
    d2 = d0;
    d0 = 2.0 * x * d0 - d1 + Tk[j];
    d1 = d2;
  }
  return x * d0 - d1 + 0.5 * Tk[0];
}


double LagrangeWeight(int interpolation_order, double x, int m) {
  int j;
  double num = 1., denom = 1.;
  for (j=0;j<interpolation_order;j++) {
    if(m!=j) {
      num   *= x - (1-2*(double)j/(double)(interpolation_order-1));
      denom *= -2*((double)m-(double)j)/(double)(interpolation_order-1);
    }
  }
  return num/denom;
}



/*
 * Function: ComputeTk
 * ------------------------------------------------------------------
 * Computes T_k(x) for k between 0 and n-1 inclusive.
 */
void H2_3D_Tree::ComputeTk(double x, int interpolation_order, double *vec) {
	int k;
	
	vec[0] = 1;
	vec[1] = x;
	
	for (k=2;k<interpolation_order;k++)
		vec[k] = 2.0*x*vec[k-1] - vec[k-2];
}

/*
 * Function: BuildFMMHierarchy
 * ------------------------------------------------------------------
 * Builds the FMM hierarchy with l levels.
 *
 */
void H2_3D_Tree::BuildFMMHierarchy(nodeT **A, int tree_level, int interpolation_order, doft *cutoff, doft *dof, int leafIndex, nodeT** indexToLeafPointer, std::vector<nodeT*>& cellPointers) {
	int i;
	double L, halfL, quarterL;
	vector3 center, left, right;

  cellPointers.push_back(*A); // collect all pointers to cells
  (*A)->cur_level = this->tree_level - tree_level; // can only work with a single thread

	if (tree_level > 0) {
        // Compute the length and center coordinate of the children cells
		L        = (*A)->length;
		halfL    = 0.5*L;
		quarterL = 0.25*L;
		center   = (*A)->center;
		left.x   = center.x - quarterL;
		left.y   = center.y - quarterL;
		left.z   = center.z - quarterL;
		right.x  = center.x + quarterL;
		right.y  = center.y + quarterL;
		right.z  = center.z + quarterL;
		
        // Add children nodes to node A
		for (i=0;i<8;i++) {
            // Determine the center of the child cell
			if (i<4) {
				center.x = left.x;
				
				if (i<2)
					center.y = left.y;
				else
					center.y = right.y;
				
			} else {
				center.x = right.x;
				
				if (i<6)
					center.y = left.y;
				else
					center.y = right.y;
			}
			
			if (i%2 == 0)
				center.z = left.z;
			else
				center.z = right.z;
			
            // Create the child cell
			NewNode(&((*A)->leaves[i]),center,halfL,interpolation_order);
			(*A)->leaves[i]->parent = *A;
		}
		
    // Recursively build octree if there is a subsequent level
    // #pragma omp parallel for private(i)
		for (i=0;i<8;i++)
			BuildFMMHierarchy(&((*A)->leaves[i]),tree_level-1,interpolation_order,cutoff,dof,leafIndex + i*pow(8, tree_level-1), indexToLeafPointer, cellPointers);
	} else {
    // handle index of leaves
    (*A)->leafIndex = leafIndex;
    indexToLeafPointer[leafIndex] = *A;
  }
}

/*
 * Function: NewNode
 * ------------------------------------------------------------------
 * Dynamically allocates space for a new node A in the octree and
 * initializes the quantities in the node.
 *
 */
void H2_3D_Tree::NewNode(nodeT **A, vector3 center, double L, int interpolation_order) {

    // Initialization
	*A = (nodeT *)malloc(sizeof(nodeT));
	(*A)->leafIndex = -1;
    // Initializes the child, neighbor, interaction, and parent nodes to NULL
	int i;
	for (i=0;i<8;i++)
		(*A)->leaves[i] = NULL;
    for (i=0;i<27;i++) {
        (*A)->neighborComputed[i] = false;
        (*A)->neighbors[i] = NULL;
    }
	(*A)->parent = NULL;
	
    (*A)->targetval  = NULL;
    (*A)->sourceval = NULL;
    (*A)->proxysval = NULL;
    (*A)->sourcefre = NULL;
	(*A)->targetlist  = NULL;
	(*A)->sourcelist = NULL;
	(*A)->center     = center;
	(*A)->length     = L;
	(*A)->Nf         = 0;
	(*A)->Ns         = 0;
	(*A)->ineigh     = 0;
	(*A)->iinter     = 0;
    (*A)->chargeComputed = false;
    (*A)->locationComputed = false;
    (*A)->location = NULL;
    (*A)->charge = NULL;
    (*A)->max_neighbor_Ns = 0;
    (*A)->nodePhi = NULL;
}


void H2_3D_Tree::get_Charge(nodeT*& node, double* q, int N, int m){
    if(node->chargeComputed==true){
        return;
    }
    else{
        node->chargeComputed    =   true;
        node->charge = (double*) malloc(node->Ns*m*sizeof(double) );  // TODO: change 1 to m for later 
        for(int i=0;i<m;i++)
            for(int k=0;k<node->Ns;++k){
                node->charge[i*node->Ns+k] =  q[i*N+node->sourcelist[k]];
            }
    }
}

void H2_3D_Tree::get_Location(nodeT*& node, vector3 *source){
    if(node->locationComputed==true){
        return;
    }
    else{
        node->locationComputed    =   true;
        node->location = (vector3*) malloc(node->Ns * sizeof(vector3));
        int l;
        for(int k=0;k<node->Ns;++k){
            l = node->sourcelist[k];
            node->location[k].x = source[l].x;
            node->location[k].y = source[l].y;
            node->location[k].z = source[l].z;
        }
    }
}

void H2_3D_Tree::FreeNode(nodeT *A) {
	int i;
	
    // Free all child nodes first
	for (i=0;i<8;i++) {
		if (A->leaves[i] != NULL) {
			FreeNode(A->leaves[i]);
		}
    }
	
	if (A->targetval != NULL)
        free(A->targetval), A->targetval=NULL;
	if (A->sourceval != NULL)
        free(A->sourceval), A->sourceval=NULL;
	if (A->proxysval != NULL)
        free(A->proxysval), A->proxysval=NULL;
	if (A->sourcefre != NULL)
	    free(A->sourcefre), A->sourcefre=NULL;
	if (A->targetlist != NULL)
        free(A->targetlist), A->targetlist=NULL;
	if (A->sourcelist != NULL)
        free(A->sourcelist), A->sourcelist=NULL;
    if (A->location != NULL)
        free(A->location), A->location=NULL;
    if (A->charge != NULL)
        free(A->charge), A->charge=NULL;
    if (A->targetlist != NULL)
        free(A->targetlist), A->targetlist=NULL;
    if (A->nodePhi != NULL)
        free(A->nodePhi), A->nodePhi=NULL;
	free(A);
}


H2_3D_Tree::~H2_3D_Tree() {
    if (tree != NULL) {
        FreeNode(tree);
        tree = NULL;
    }
    if (Kweights!= NULL) {
        free(Kweights);
        Kweights = NULL;
    }
    if (Cweights!= NULL) {
        free(Cweights);
        Cweights = NULL;
    }
    if (Tkz!= NULL) {
        free(Tkz);
        Tkz = NULL;
    }
    if (K!= NULL) {
        free(K);
        K = NULL;
    }
    if (U!= NULL) {
        free(U);
        U = NULL;
    }
    if (VT!= NULL) {
        free(VT);
        VT = NULL;
    }
    dof = NULL;
}




