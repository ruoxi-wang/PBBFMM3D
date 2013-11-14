
/*! \file	H2_3D_Tree.cpp
*/

#include"H2_3D_Tree.hpp"
#include"bbfmm.h"

H2_3D_Tree::H2_3D_Tree(doft *dof, double L, int level, int n,  double epsilon){
    this->dof   =   dof;
    this->L     =   L;
    this->level =   level;
    this->n    =   n;
    this->epsilon   =   epsilon;
    alpha = 0;
    n2 = n*n;         // n2 = n^2
	n3 = n2*n;       // n3 = n^3
	dofn3_s = dof->s * n3;
	dofn3_f = dof->f * n3;
    
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
    FMMSetup(&tree,Tkz,Kweights,Cweights,L,&cutoff,
			 n,epsilon,dof,level,Kmat,
             Umat,Vmat,&Ucomp,&Vcomp, skipLevel, alpha);
     K  = (double *)malloc(cutoff.f * (316*cutoff.s) * sizeof(double));
     U  = (double *)malloc(dofn3_f * cutoff.f * sizeof(double));
     VT = (double *)malloc(cutoff.s * dofn3_s * sizeof(double));
    FMMReadMatrices(K,U,VT,&cutoff,n,dof,Kmat,Umat,Vmat, level, homogen, skipLevel);
}

/*
 * Function: FMMSetup
 * -----------------------------------------------------------------
 * Prepare for the FMM calculation by setting the parameters, computing
 * the weight matrices, pre-computing the SVD (if necessary), reading
 * in the necessary matrices, and building the FMM hierarchy.
 */
void H2_3D_Tree::FMMSetup(nodeT **A, double *Tkz,  double *Kweights,
                          double *Cweights, double boxLen, doft *cutoff,
                          int n, double epsilon, doft * dof,  int treeLevel, char *Kmat, char *Umat, char *Vmat,
                          double *Ucomp, double *Vcomp,int& skipLevel, double alpha) {
    
	vector3 center;
    setHomogen(kernelType);
    homogen = -homogen;
	sprintf(Kmat,"./../output/%sK%d.bin",kernelType.c_str(),n);
	sprintf(Umat,"./../output/%sU%d.bin",kernelType.c_str(),n);
	sprintf(Vmat,"./../output/%sV%d.bin",kernelType.c_str(),n);
    // Compute the Chebyshev weights and sets up the lookup table
	ComputeWeights(Tkz,Ktable,Kweights,Cweights,n);
    // Precompute the SVD of the kernel interaction matrix (if
    // necessary)
    int i;
    FILE *fK, *fU, *fV;
    fK = fopen(Kmat, "rb");
    fU = fopen(Umat, "rb");
    fV = fopen(Vmat, "rb");
    if (fK == NULL || fU == NULL || fV == NULL) { // Create files
        if (fK!=NULL)
            fclose(fK);
        if (fU!=NULL)
            fclose(fU);
        if (fV!=NULL)
            fclose(fV);
        
        printf("Pre-Compute files do not exit. Creating now ...\n");
        
        if (homogen > 1e-9) // homogeneous kernel
            ComputeKernelSVD(Kweights, n, epsilon, dof,
                             Kmat, Umat, Vmat, symmetry, Ucomp, Vcomp,alpha,
                             1);
        
        else { // non-homogeneous kernel
            
            // Open new files for writing boxLen and treeLevel
            
            fK = fopen(Kmat, "wb");
            fwrite(&boxLen, sizeof(double), 1, fK);
            fwrite(&treeLevel, sizeof(int), 1, fK);
            fclose(fK);
            
            
            double boxLenLevel = boxLen/4; // first FMM level
            for (i=2; i<=treeLevel; i++, boxLenLevel/=2)
                // FMM starts from the second level
                ComputeKernelSVD(Kweights, n,  epsilon, dof,
                                 Kmat, Umat, Vmat, symmetry, Ucomp, Vcomp,alpha,
                                 boxLenLevel);
        }
    }
    
    //non-homogen
    else if (homogen < 1e-9) { // check if the file is usable
        
        fK = fopen(Kmat, "rb");
        double fileBoxLen;
        int fileTreeLevel;
        i  = fread(&fileBoxLen, sizeof(double), 1, fK);
        i += fread(&fileTreeLevel, sizeof(int), 1, fK);
        if (i != 2)
            printf("fread error in FMMSetup().\n");
        
        int count = 0;
        while (fileBoxLen > boxLen +1e-9) {
            fileBoxLen /= 2;
            count ++;
        }
        if (fabs(boxLen-fileBoxLen) < 1e-9 && treeLevel + count <=
            fileTreeLevel) {
            skipLevel = count; // count * Ksize
            printf("Reading pre-compute files ...\n");
        }
        else { // Recreate the files
            
            printf("Recreating pre-compute files now ...\n");
            
            fK = fopen(Kmat, "wb");
            fwrite(&boxLen, sizeof(double), 1, fK);
            fwrite(&treeLevel, sizeof(int), 1, fK);
            fclose(fK);
            
            int i;
            double boxLenLevel = boxLen/4; // first FMM level
            for (i=2; i<=treeLevel; i++, boxLenLevel/=2)
                // FMM starts from the second level
                ComputeKernelSVD(Kweights, n, epsilon, dof,
                                 Kmat, Umat, Vmat, symmetry, Ucomp, Vcomp,alpha,
                                 boxLenLevel);
        }
    }
    else
        printf("Reading pre-compute files ...\n");
    
    fU = fopen(Umat,"rb");
    int j = fread(&(cutoff->f), sizeof(int), 1, fU);
    fclose(fU);

    fV = fopen(Vmat,"rb");
    j += fread(&(cutoff->s), sizeof(int), 1, fV);
    fclose(fV);

    if (j != 2)
        printf("fread() error in FMMSetup().\n");
    // Builds the FMM hierarchy
	center.x = 0;
	center.y = 0;
	center.z = 0;
	(*A) = NULL;

	NewNode(A,center,L,n);
	BuildFMMHierarchy(A,treeLevel,n,cutoff,dof);
}


/*
 * Function: FMMReadMatrices
 * ------------------------------------------------------------------
 * Read in the kernel interaction matrix M and the matrix of singular
 * vectors U.
 */
void H2_3D_Tree::FMMReadMatrices(double *K, double *U, double *VT, doft *cutoff, int n, doft *dof,char *Kmat, char *Umat, char *Vmat, int treeLevel, double homogen, int skipLevel) {
	int i = 0;
    FILE *ptr_file;
    int preCompLevel = (treeLevel-2)* (!homogen) + 1;
    
    int Ksize = cutoff->f * (316*cutoff->s);
    int Usize = cutoff->f * dof->f *n*n*n;
    int Vsize = cutoff->s * dof->s *n*n*n;
	
    // Read in kernel interaction matrix K
    ptr_file = fopen(Kmat,"rb");
    fseek(ptr_file, (1*sizeof(int) + 1*sizeof(double)) *(!homogen) +
          Ksize *skipLevel *sizeof(double), SEEK_SET);
    i += fread(K, sizeof(double), Ksize *preCompLevel, ptr_file);
    fclose(ptr_file);
    
    // Read in matrix of singular vectors U
    ptr_file = fopen(Umat,"rb");
    fseek(ptr_file, 1*sizeof(int) + Usize *skipLevel *sizeof(double),
          SEEK_SET);
    i += fread(U, sizeof(double), Usize *preCompLevel, ptr_file);
    fclose(ptr_file);
	
    // Read in matrix of singular vectors VT
    ptr_file = fopen(Vmat,"rb");
    fseek(ptr_file, 1*sizeof(int) + Vsize *skipLevel *sizeof(double),
          SEEK_SET);
    i += fread(VT, sizeof(double), Vsize *preCompLevel, ptr_file);
    fclose(ptr_file);
    
    int totleNum = (Ksize + Usize + Vsize) *preCompLevel;
    if (i != totleNum)
        printf("fread error in FMMReadMatrices!\n Expected numer:%d,"
               "numbers read: %d\n", totleNum, i);
}

/*
 * Function: ComputeKernelSVD
 * ---------------------------------------------------------------------
 * Computes the kernel for 316n^6 interactions between Chebyshev nodes
 * and then computes the SVD of the kernel matrix.
 * symmetry = 0 no symmetric property, 1 symmetric kernel or -1 anti-
 * symmetric kernel
 */


void H2_3D_Tree::ComputeKernelSVD(double *Kweights, int n,double epsilon, doft *dof, char*Kmat, char *Umat, char *Vmat, int symm,  double *Ucomp,double *Vcomp, double alphaAdjust, double boxLen) {
    
    static int callTime = -1;
    callTime += 1; // callTime = 0 for the first time called
    
    int i, j, l, m, m1, k1, k2, k3, l1, l2, l3, z;
    int count, count1, count2, count3;
    vector3 scenter, vtmp;
    double sweight, fweight;
    double pi = M_PI;
	
    int n3 = n*n*n;            // n3 = n^3
    int dofn3_s = dof->s * n3;
    int dofn3_f = dof->f * n3;
    int dof2n6 = dofn3_s * dofn3_f; // Total size
    int Sigma_size;
    doft cutoff;
	
    double *K0, *U0, *Sigma, *VT0;
    double *nodes, *kernel, *work;
    vector3 *fieldpos, *sourcepos;
	
    K0 = (double *)malloc(316 * dof2n6 * sizeof(double));
    kernel = (double *)malloc(dof2n6 * sizeof(double));
    fieldpos  = (vector3 *)malloc(n3 * sizeof(vector3));
    sourcepos = (vector3 *)malloc(n3 * sizeof(vector3));
	
    // 316 M2L operators
    double* Kcell[316];
    for (z = 0; z < 316; ++z)
        Kcell[z] = (double *) malloc(dof2n6 * sizeof(double));
	
	
    // Compute Chebyshev nodes of T_n(x)
    double scale = 1+alphaAdjust;
    nodes = (double *)malloc(n * sizeof(double));
    for (m=0;m<n;m++)
        nodes[m] = cos(pi*((double)m+0.5)/(double)n) * scale;
	
    // Compute the locations of the field points in a unit cube
    count = 0;
    for (l1=0;l1<n;l1++) {
        vtmp.x = 0.5 * nodes[l1] * boxLen;
        for (l2=0;l2<n;l2++) {
            vtmp.y = 0.5 * nodes[l2] * boxLen;
            for (l3=0;l3<n;l3++) {
                fieldpos[count].x = vtmp.x;
                fieldpos[count].y = vtmp.y;
                fieldpos[count].z = 0.5 * nodes[l3] * boxLen;
                count++;
            }
        }
    }
    
    // Compute the kernel values for interactions with all 316 cells
    int countM2L=0;
    int symmNum = 158*(2-abs(symm)); // symmNum = 158, 316 for symm=pm1, 0
    int col, row, idx1, idx2;
    for (count=0, k1=-3;k1<4;k1++) {
        scenter.x = (double)k1;
        for (k2=-3;k2<4;k2++) {
            scenter.y = (double)k2;
            for (k3=-3;k3<4;k3++) {
                scenter.z = (double)k3;
                if (abs(k1) > 1 || abs(k2) > 1 || abs(k3) > 1) {
                    if (countM2L < symmNum) {
                        
                        // Compute the locations of the source points in the cell
                        for (count1=0, l1=0;l1<n;l1++) {
                            vtmp.x = (scenter.x + 0.5 * nodes[l1]) * boxLen;
                            for (l2=0;l2<n;l2++) {
                                vtmp.y = (scenter.y + 0.5 * nodes[l2]) * boxLen;
                                for (l3=0; l3<n; l3++, count1++) {
                                    sourcepos[count1].x = vtmp.x;
                                    sourcepos[count1].y = vtmp.y;
                                    sourcepos[count1].z = (scenter.z + 0.5
                                                           * nodes[l3]) * boxLen;
                                }
                            }
                        }
                        
                        // Compute the kernel at each of the field Chebyshev nodes
                        EvaluateKernelCell(fieldpos, sourcepos, n3, n3, dof,
                                            kernel);
                        
                        // Copy the kernel values to the appropriate location
                        for (count1=0, count2=0, l1=0; l1<n3; l1++, count2++) {
                            sweight = Kweights[count2];
                            for (l=0;l<dof->s;l++) {
                                for (count3=0, m1=0; m1<n3; m1++, count3++) {
                                    fweight = Kweights[count3];
                                    for (m=0; m<dof->f; m++, count++, count1++) {
                                        K0[count] = kernel[count1]
                                        /(sweight*fweight);
                                    }
                                }
                            }
                        }
                        countM2L++;
                    }
                    else { // Use kernel symmetry
                        for (col=0; col<dofn3_s; col++)
                            for (row=0; row<dofn3_f; row++) {
                                
                                idx1 = (315-countM2L)*dof2n6 + col*dofn3_f + row;
                                idx2 = countM2L*dof2n6 +
                                (row/dof->f*dof->s + col%dof->s) *dofn3_f
                                + (col/dof->s*dof->f + row%dof->f);
                                
                                // symm=-1,1 for anti-symmetry and symmetry
                                K0[idx2] = symm * K0[idx1];
                                //K0[idx2] = K0[idx1];
                            }
                        countM2L++;
                    }
                }
            }
        }
    }
    
	
    // Extract the submatrix for each of the 316 cells
    count = 0;
    for (i=0;i<316;i++) {
        for (j=0;j<dof2n6;j++) {
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
    int max_dof = dof->s > dof->f ? dof->s : dof->f;
    lwork = 5*max_dof*n3; // Change
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
    if (callTime == 0) { // The first time called
        
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
            printf("fread error in ComputeKernelSVD().\n");
        fclose(ptr_file);
    }
    
    *Ucomp = ((double)cutoff.f)/Sigma_size;
    
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
    
    if (callTime == 0) {
        
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
            printf("fread error in ComputeKernelSVD().\n");
        fclose(ptr_file);
    }
    
    *Vcomp = ((double)cutoff.s)/Sigma_size;
    
    ptr_file = fopen(Vmat, "ab");	    
    for (j=0;j<dofn3_s;j++) { // column
        count1 = j*dofn3_s;
        fwrite(VT0+count1, sizeof(double), cutoff.s, ptr_file);
    }
    fclose(ptr_file);
    
    
    /** Computing the compressed kernel using the orthonormal basis U and VT.
     **/
    char *transa, *transb;
    transa = new char[2];
    transb = new char[2];
    int cutoff2 = cutoff.f * cutoff.s;	
    double alpha=1, beta=0;
    double *Ecell, *KV;
    Ecell = (double *)malloc(cutoff2 * sizeof(double));	
    KV    = (double *)malloc(dofn3_f * cutoff.s * sizeof(double));
	
    ptr_file = fopen(Kmat,"ab");
    for (i=0;i<316;i++) {       
        
        /* Compute K V:
         * K  is dofn3_f  x dofn3_s
         * VT is cutoff.s x dofn3_s
         * V  is dofn3_s  x cutoff.s 
         * KV is dofn3_f  x cutoff.s
         * (Notice that VT is a submatrix of VT0)
         */
        transa[0] = 'n';
        transb[0] = 't';
        dgemm_(transa, transb, &dofn3_f, &cutoff.s, &dofn3_s, &alpha, 
               Kcell[i], &dofn3_f, VT0, &dofn3_s, &beta, KV, &dofn3_f);
		
        /* Compute U^T K V:
         * KV is dofn3_f x cutoff.s
         * U  is dofn3_f x cutoff.f
         * UT is cutoff.f x dofn3_f
         * U^T K V is cutoff.f x cutoff.s
         * (Notice that U is a submatrix of U0)
         */
        transa[0] = 't';
        transb[0] = 'n';
        dgemm_(transa, transb, &cutoff.f, &cutoff.s, &dofn3_f, &alpha, 
               U0, &dofn3_f, KV, &dofn3_f, &beta, Ecell, &cutoff.f);
		
        fwrite(Ecell, sizeof(double), cutoff2, ptr_file);
    }
    fclose(ptr_file);
    
    free(K0);
    free(VT0);
    free(U0);
    free(Sigma);
    free(nodes);
    free(kernel);
    free(fieldpos);
    free(sourcepos);
    free(work);
    free(KV);
    free(Ecell);
    delete []transa;
    delete []transb;

    for (z = 0; z < 316; ++z)
        free(Kcell[z]);
    
}

/*
 * Function: EvaluateKernelCell
 * -------------------------------------------------------------------
 * Evaluates the kernel for interactions between a pair of cells.
 * M2L operator initialization
 */
void H2_3D_Tree::EvaluateKernelCell(vector3 *field, vector3 *source, int Nf,
                        int Ns, doft *dof, double *kernel) {
	int i, j, k, l, count, count_kernel;
	int dof_f = dof->f;
	int dof_s = dof->s;
	int dofNf = dof->f * Nf;
	int dof2  = dof->f * dof->s;
	double * Kij;
	
	int LDA_kernel = dofNf;
	
	Kij = (double*) malloc(dof2 * sizeof(double));
	
	for (j=0;j<Ns;j++) {
		for (i=0;i<Nf;i++) {
            EvaluateKernel(field[i],source[j],Kij, dof);
			
			count_kernel = dof_f * i + LDA_kernel * dof_s * j;
			count = 0;
			for (k=0;k<dof_s;k++)
				for (l=0;l<dof_f;l++, count++)
                /* Column-major storage */
					kernel[count_kernel + k * LDA_kernel + l] = Kij[count];
		}
	}
    
	free(Kij);
}

/*
 * Function: ComputeWeights
 * ------------------------------------------------------------------
 * Computes the weights for the Chebyshev nodes of all children cells
 * (identical for all cells and all levels so just compute once and
 * store in memory) and set up the lookup table.
 */

void H2_3D_Tree::ComputeWeights(double *Tkz, int *Ktable, double *Kweights,
					double *Cweights, int n) {
	int i, k, m, l1, l2, l3, count, ncell, ninteract;
	double tmp1, tmp2;
	vector3 vtmp;
	
	double *nodes, *vec;
	vector3 *fieldt, *Sn;
	nodes = (double *)malloc(n * sizeof(double));
	vec   = (double *)malloc(n * sizeof(double));
	int n3 = n*n*n;                     // n3 = n^3
	int Nc = 2*n3;                      // Number of child Chebyshev nodes
	fieldt = (vector3 *)malloc(Nc * sizeof(vector3));
    // Chebyshev-transformed coordinates
	Sn     = (vector3 *)malloc(n * Nc * sizeof(vector3));
	double pi = M_PI;
	
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

    // Compute the n Chebyshev nodes of T_n(x)
	for (m=0; m<n; m++)
		nodes[m] = cos(pi*((double)m+0.5)/(double)n);

    // Evaluate the Chebyshev polynomials of degree 0 to n-1 at the nodes
	for (m=0; m<n; m++) {
		ComputeTk(nodes[m],n,vec);
		i = m*n;
		for (k=0;k<n;k++) {
			Tkz[i+k] = vec[k];
        }
	}

    // Compute the weights for the kernel matrix K
	count = 0;
	for (l1=0;l1<n;l1++) {
		tmp1 = 1./sqrt(1-nodes[l1]*nodes[l1]);
		for (l2=0;l2<n;l2++) {
			tmp2 = tmp1/sqrt(1-nodes[l2]*nodes[l2]);
			for (l3=0;l3<n;l3++) {
				Kweights[count] = tmp2/sqrt(1-nodes[l3]*nodes[l3]);
				count++;
			}
		}
	}
	
    // Map all Chebyshev nodes from the children cells to the parent
	k = 0;
	for (i=0;i<2;i++) {
        // Determine the mapping function for the specific child cell
		vtmp.x = -1;
		vtmp.y = -1;
		
		if (i == 0)
			vtmp.z = -1;
		else
			vtmp.z = 1;
		
		for (l1=0;l1<n;l1++) {
			for (l2=0;l2<n;l2++) {
				for (l3=0;l3<n;l3++) {
					fieldt[k].x = 0.5*(nodes[l1] + vtmp.x);
					fieldt[k].y = 0.5*(nodes[l2] + vtmp.y);
					fieldt[k].z = 0.5*(nodes[l3] + vtmp.z);
					k++;
				}
			}
		}
	}
    
    // Compute Sc, the mapping function for the field points
	ComputeSn(fieldt,Tkz,n,Nc,Sn);
	
    // Extract out the Chebyshev weights
	count = 0;
	for (l1=0;l1<n;l1++) {
		k = l1*Nc;
		for (l2=0;l2<n;l2++) {
			Cweights[count] = Sn[k+l2].z;
			count++;
		}
	}
	for (l1=0;l1<n;l1++) {
		k = l1*Nc;
		for (l2=0;l2<n;l2++) {
			Cweights[count] = Sn[k+n3+l2].z;
			count++;
		}
	}
	free(nodes);
	free(vec);
	free(fieldt);
	free(Sn);
}




/*
 * Function: ComputeSn
 * ------------------------------------------------------------------
 * Computes S_n(x_m,x_i) for all Chebyshev node-point pairs using
 * Clenshaw's recurrence relation.
 * Sn(xm, xi) = -1/2 + sum_0^(n-1) Tk(xm) Tk(xi)
 * Notice there is another prefactor = 2/n that will be applied later
 */
void H2_3D_Tree::ComputeSn(vector3 *point, double *Tkz, int n, int N, vector3 *Sn) {
	int i, j, k, m;
	double vec[n], d[n+2], x;
	
	for (m=0;m<n;m++) {
        // Extract T_k for the Chebyshev node x_m
		k = m*n;
		for (j=0;j<n;j++)
			vec[j] = Tkz[k+j];
		
        // Compute S_n for each direction independently using Clenshaw
		k = m*N;
		for (i=0;i<N;i++) {
			x = point[i].x;
			d[n] = d[n+1] = 0.;
			for (j=n-1;j>0;j--)
				d[j] = 2.0*x*d[j+1] - d[j+2] + vec[j];
			Sn[k+i].x = x*d[1] - d[2] + 0.5*vec[0];
			
			x = point[i].y;
			d[n] = d[n+1] = 0.;
			for (j=n-1;j>0;j--)
				d[j] = 2.0*x*d[j+1] - d[j+2] + vec[j];
			Sn[k+i].y = x*d[1] - d[2] + 0.5*vec[0];
			
			x = point[i].z;
			d[n] = d[n+1] = 0.;
			for (j=n-1;j>0;j--)
				d[j] = 2.0*x*d[j+1] - d[j+2] + vec[j];
			Sn[k+i].z = x*d[1] - d[2] + 0.5*vec[0];
		}
	}
}
/*
 * Function: ComputeTk
 * ------------------------------------------------------------------
 * Computes T_k(x) for k between 0 and n-1 inclusive.
 */
void H2_3D_Tree::ComputeTk(double x, int n, double *vec) {
	int k;
	
	vec[0] = 1;
	vec[1] = x;
	
	for (k=2;k<n;k++)
		vec[k] = 2.0*x*vec[k-1] - vec[k-2];
}

/*
 * Function: BuildFMMHierarchy
 * ------------------------------------------------------------------
 * Builds the FMM hierarchy with l levels.
 *
 */
void H2_3D_Tree::BuildFMMHierarchy(nodeT **A, int level, int n, doft *cutoff, doft *dof) {
	int i;
	double L, halfL, quarterL;
	vector3 center, left, right;
	
	if (level > 0) {
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
			NewNode(&((*A)->leaves[i]),center,halfL,n);
			(*A)->leaves[i]->parent = *A;
		}
		
        // Recursively build octree if there is a subsequent level
		for (i=0;i<8;i++)
			BuildFMMHierarchy(&((*A)->leaves[i]),level-1,n,cutoff,dof);
	}
}

/*
 * Function: NewNode
 * ------------------------------------------------------------------
 * Dynamically allocates space for a new node A in the octree and
 * initializes the quantities in the node.
 *
 */
void H2_3D_Tree::NewNode(nodeT **A, vector3 center, double L, int n) {

    // Initialization
	*A = (nodeT *)malloc(sizeof(nodeT));
	
    // Initializes the child, neighbor, interaction, and parent nodes to NULL
	int i;
	for (i=0;i<8;i++)
		(*A)->leaves[i] = NULL;
	(*A)->parent = NULL;
	
    (*A)->fieldval  = NULL;
    (*A)->sourceval = NULL;
    (*A)->proxysval = NULL;
    
	(*A)->fieldlist  = NULL;
	(*A)->sourcelist = NULL;
	(*A)->center     = center;
	(*A)->length     = L;
	(*A)->Nf         = 0;
	(*A)->Ns         = 0;
	(*A)->ineigh     = 0;
	(*A)->iinter     = 0;
}

/*void H2_3D_Tree::get_Compression_Rate(double* Ucomp, double* Vcomp) {
    *Ucomp = this->Ucomp;
    *Vcomp = this->Vcomp;
}*/

/*void H2_3D_Tree::EvaluateKernelMulti(vector3 fieldpos, vector3 sourcepos,
                                     double *K, doft *dof, int m) {
    doft dof_kernel;
    dof_kernel.f = dof->f;
    dof_kernel.s = dof->s*m;
    EvaluateKernel(fieldpos, sourcepos, K, dof_kernel);
}*/

void H2_3D_Tree::FreeNode(nodeT *A) {
	int i;
	
    // Free all child nodes first
	for (i=0;i<8;i++) {
		if (A->leaves[i] != NULL) {
			FreeNode(A->leaves[i]);
		}
	}
	
    // Free the arrays for the field and source values
	if (A->fieldval != NULL)
        free(A->fieldval), A->fieldval=NULL;
	if (A->sourceval != NULL)
        free(A->sourceval), A->sourceval=NULL;
	if (A->proxysval != NULL)
        free(A->proxysval), A->proxysval=NULL;
	
    // Free the field and source lists
	if (A->fieldlist != NULL)
        free(A->fieldlist), A->fieldlist=NULL;
	if (A->sourcelist != NULL)
        free(A->sourcelist), A->sourcelist=NULL;
	
    // Last free the node
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


