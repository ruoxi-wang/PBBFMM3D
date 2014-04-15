
/*! \file compute.hpp
 */
#ifndef __compute_hpp__
#define __compute_hpp__

#include"bbfmm.h"
#include "environment.hpp" 
#include"kernel_Types.hpp"
#include "rfftw.h"

template <typename T>
class H2_3D_Compute {
public:
    H2_3D_Compute(T* FMMTree, vector3 * field, vector3 *source, int Ns, int Nf, double* charge,int m,double* stress);
    T* FMMTree;
    vector3 * field;
    vector3 * source;
    int Ns;
    int Nf;
    double* charge;
    int m;
    void FMMDistribute(nodeT **A, vector3 *field, vector3 *source, int Nf,
                       int Ns, int level);
    void FMMCompute(nodeT **A, vector3 *field, vector3 *source, double *charge,
                    double *K, double *U, double *VT, double *Tkz, int *Ktable,
                    double *Kweights, double *Cweights, double homogen,
                    doft *cutoff, int n, doft *dof,double*stress, int use_chebyshev, rfftw_plan p_r2c,
                    rfftw_plan p_c2r);
    void UpwardPass(nodeT **A, vector3 *source, double *q, double *Cweights, double *Tkz, double *VT, double *Kweights, doft *cutoff, int n, doft *dof,  int use_chebyshev, rfftw_plan p_r2c);
    void FMMInteraction(nodeT **A, double *E, int *Ktable, double *U,
                        double *VT, double *Kweights, int n, doft *dof,
                        doft *cutoff, double homogen, int curTreeLevel, int use_chebyshev, rfftw_plan p_c2r);
    void DownwardPass(nodeT **A, vector3 *field, vector3 *source,
                      double *Cweights, double *Tkz, double *q, doft *cutoff,
                      int n, doft *dof, double *stress, int use_chebyshev);
    void EvaluateField(vector3 *field, vector3 *source, double *q, int Nf,
                       int Ns, doft *dof, double *fieldval);
    void Local2Local(int n, double *r, double *F, double *Fchild, doft *dof, double *Cweights, int use_chebyshev);
    void Moment2Local(int n, double *R, double *cell_mpCoeff, double *FFCoeff,
                 double *E, int *Ktable, doft *dof, doft *cutoff, double *VT,
                      double *Kweights, int use_chebyshev);
    void Moment2Moment(int n, double *r, double *Schild, double *SS, doft *dof, double *Cweights);
    void InteractionList(nodeT **A, int levels);
    
    void DistributeSources(nodeT **A, vector3 *source, int *sourcelist,
                           int levels);
    void DistributeFieldPoints(nodeT **A, vector3 *field, int *fieldlist,
                               int levels);
    
    void FrequencyProduct(int N, double *Afre, double *xfre, double *res);
    void FreeNode(nodeT *A);
    ~H2_3D_Compute();
};


template <typename T>
H2_3D_Compute<T>::H2_3D_Compute(T * FMMTree,vector3 * field, vector3 *source, int Ns, int Nf, double* charge, int m, double* stress) {
    if (FMMTree->computed == true) {
        if (FMMTree->tree != NULL) {
            FMMTree->FreeNode(FMMTree->tree);
            FMMTree->tree = NULL;
        }
        vector3 center;
        center.x = 0;
        center.y = 0;
        center.z = 0;
        FMMTree->NewNode(&FMMTree->tree,center,FMMTree->L,FMMTree->n);
        FMMTree->BuildFMMHierarchy(&FMMTree->tree,FMMTree->level,FMMTree->n,&FMMTree->cutoff,FMMTree->dof);

    }
    this->FMMTree   =   FMMTree;
    this-> field    =   field;
    this-> source   =   source;
    this-> Ns       =   Ns;
    this-> Nf       =   Nf;
    this->charge    =   charge;
    this->m         =   m;
    FMMDistribute(&FMMTree->tree,field,source,Nf,Ns, FMMTree->level);
    for (int i = 0; i < m; i++) {
        FMMCompute(&FMMTree->tree,field,source,&charge[i*Ns*FMMTree->dof->s],FMMTree->K,FMMTree->U,FMMTree->VT,FMMTree->Tkz,FMMTree->Ktable,FMMTree->Kweights,FMMTree->Cweights,
                   FMMTree->homogen,&FMMTree->cutoff,FMMTree->n,FMMTree->dof,&stress[i*Nf*FMMTree->dof->f], FMMTree->use_chebyshev,FMMTree->p_r2c, FMMTree->p_c2r);
    }
    FMMTree->computed = true;

}

/*
 * Function: FMMDistribute
 * ------------------------------------------------------------------
 * Distribute the field points and sources to the appropriate location
 * in the FMM hierarchy and sets up the interaction list.
 */
template <typename T>
void H2_3D_Compute<T>::FMMDistribute(nodeT **A, vector3 *field, vector3 *source, int Nf,int Ns, int level) {
	int i;
	vector3 cshift;
	int *fieldlist, *sourcelist;
	fieldlist = (int *)malloc(Nf * sizeof(int));
	sourcelist = (int *)malloc(Ns * sizeof(int));
	
    // Initialize the point distribution for the root node
	for (i=0;i<Nf;i++)
		fieldlist[i] = i;
	for (i=0;i<Ns;i++)
		sourcelist[i] = i;
	(*A)->Nf = Nf;
	(*A)->Ns = Ns;
	
    // Distribute all of the sources and field points to the appropriate cell
	DistributeFieldPoints(A,field,fieldlist,level);
	DistributeSources(A,source,sourcelist,level);
	
    // Construct the interaction list for the entire octree
	(*A)->neighbors[0] = *A;
	(*A)->ineigh = 1;
	cshift.x = 0;
	cshift.y = 0;
	cshift.z = 0;
	(*A)->cshiftneigh[0] = cshift;
	InteractionList(A,level);  // build the interaction that need to be computed
	
    
	free(fieldlist);
    fieldlist = NULL;
	free(sourcelist);
    fieldlist = NULL;
}

template <typename T>
void H2_3D_Compute<T>::FMMCompute(nodeT **A, vector3 *field, vector3 *source, double *charge, double *K, double *U, double *VT, double *Tkz, int *Ktable, double *Kweights, double *Cweights, double homogen, doft *cutoff, int n, doft *dof, double *stress, int use_chebyshev, rfftw_plan p_r2c, rfftw_plan p_c2r) {
	UpwardPass(A,source,charge,Cweights,Tkz,VT,Kweights,cutoff,n,dof, use_chebyshev, p_r2c);
    // Computes all of the cell interactions
	FMMInteraction(A,K,Ktable,U,VT,Kweights,n,dof,cutoff,homogen,0, use_chebyshev, p_c2r);
    DownwardPass(A,field,source,Cweights,Tkz,charge,cutoff,n,dof,stress, use_chebyshev);
}

template <typename T>
void H2_3D_Compute<T>::UpwardPass(nodeT **A, vector3 *source, double *q, double *Cweights, double *Tkz, double *VT, double *Kweights, doft *cutoff, int n, doft *dof, int use_chebyshev, rfftw_plan p_r2c) {
    
    int Ns = (*A)->Ns;                  // Number of source points
    int i, l;
	int dofn3 = dof->s * n*n*n;
	
    double prefac  = 2./n;
    if(!use_chebyshev)
        prefac = 1.;

    double prefac3 = prefac*prefac*prefac;
    
	/* First gather the sources for all children cells and
	 * then gather for the parent cell - otherwise map all
	 * sources to Chebyshev nodes */
	if ((*A)->leaves[0] != NULL) {  // Going up the tree M2M
        
        // Allocate memory for the source values and initialize
		(*A)->sourceval = (double *)malloc(dofn3*sizeof(double));
		for (l=0;l<dofn3;l++) {
            (*A)->sourceval[l] = 0;
        }
		
        // Determine which children cells contain sources
		for (i=0;i<8;i++) {
		    if ((*A)->leaves[i]->Ns != 0) {
                UpwardPass(&((*A)->leaves[i]),source,q,Cweights,Tkz,VT,
						   Kweights,cutoff,n,dof,use_chebyshev,p_r2c);
                double *Schild;
                Schild = (*A)->leaves[i]->sourceval;
                
                // Manually set the vector from child to parent
                int  idx[3];
                double r[3];
                idx[2] = i%2, idx[1] = i/2%2, idx[0] = i/4%2;
                for (l=0; l<3; l++)
                    if (idx[l] == 0)
                        r[l] = -1;
                    else
                        r[l] =  1;
                
		        double SS[dofn3];
                Moment2Moment(n, r, Schild, SS, dof, Cweights);
                
                if (! use_chebyshev) {
                    if (Schild != NULL) {
                        free(Schild), Schild = NULL;
                    }
                    (*A)->leaves[i]->sourceval = NULL;

                }
                
                for (l=0; l<dofn3; l++)
                    (*A)->sourceval[l] += SS[l];
            }
        }
        
        for (i=0; i<dofn3; i++)
            // Prefactor of Chebyshev interpolation coefficients 'Sn': see FMMTree->ComputeSn
            (*A)->sourceval[i] *= prefac3;
        
	}
    
    
    
    else {
        
        int j, k, l1, l2, l3, l4, *sourcelist = (*A)->sourcelist;
        double sum, ihalfL = 2./(*A)->length;
        vector3 sourcet[Ns], Ss[n*Ns], scenter = (*A)->center;
        
        // Map all of the source points to the box ([-1 1])^3
        for (j=0;j<Ns;j++) {
            k = sourcelist[j];
            sourcet[j].x = ihalfL*(source[k].x - scenter.x);
            sourcet[j].y = ihalfL*(source[k].y - scenter.y);
            sourcet[j].z = ihalfL*(source[k].z - scenter.z);
        }
        // Compute Ss, the mapping function for the sources
        FMMTree->ComputeSn(sourcet,Tkz,n,Ns,Ss, use_chebyshev);
        
        // Compute the source values
        (*A)->sourceval = (double *)malloc(dofn3*sizeof(double));
        double *S = (*A)->sourceval;
        l = 0;
        for (l1=0;l1<n;l1++) {
            for (l2=0;l2<n;l2++) {
                for (l3=0;l3<n;l3++) {
                    for (l4=0;l4<dof->s;l4++) {
                        sum = 0;
                        for (j=0;j<Ns;j++) {
                            k = dof->s * sourcelist[j] + l4;
                            sum += q[k]*Ss[l1*Ns+j].x*Ss[l2*Ns+j].y*Ss[l3*Ns+j].z;
                        }
                        S[l] = prefac3*sum;
                        l++;
                    }
                }
            }
        }
	}
    
    
    
    if (!use_chebyshev) {
        
        // pad the source value
        int padSize = (int)round(pow(2*n-1, 3));
        int halfSize = padSize/2;
        int l3Size = n, l1Size = (int)round(pow(l3Size, 2));
        int l2Pad = 2*n-1, l1Pad = (int)round(pow(l2Pad, 2));
        int shift, count;
        int s, l1, l2, l3, n3 = n*n*n;
        
        double *x = (*A)->sourceval;
        double *padx = (double*)calloc(padSize *dof->s, sizeof(double));
        (*A)->sourcefre = (double*)malloc(padSize  *dof->s *sizeof(double));
        
        for (i=0, count=0; i<n3; i++) {
            l3 = i % l3Size;
            l1 = i / l1Size;
            l2 = i / l3Size % l3Size;
            shift = halfSize+l1*l1Pad+l2*l2Pad+l3;
            
            for (s=0; s<dof->s; s++) {
                padx[shift] = x[count++];
                shift += padSize;
            }
        }
        
        for (s=0, shift=0; s<dof->s; s++) {
            rfftw_one(p_r2c, padx + shift, (*A)->sourcefre + shift);
            shift += padSize;
        }
        
        x=NULL;
        free(padx), padx=NULL;
    }
}


template <typename T>
void H2_3D_Compute<T>::FMMInteraction(nodeT **A, double *E, int *Ktable, double *U,
                                      double *VT, double *Kweights, int n, doft *dof,
                                      doft *cutoff, double homogen, int curTreeLevel, int use_chebyshev,
                                      rfftw_plan p_c2r) {
    
    int shift = (curTreeLevel>=2) ? curTreeLevel-2: 0;
    shift *= !homogen;
        
    int i, j, l;
    int cutoff_f = cutoff->f;
    double tmp;
    vector3 fcenter, scenter, cshift;
    double *F;
    nodeT *B;
    
    int n3 = n*n*n;                      // n3 = n^3
    int dofn3_f  = dof->f * n3, dofn3_s = dof->s * n3;;
    double L     = (*A)->length;             // Length of cell
    double iL    = 1.0/L;                   // Inverse-length
    double scale = pow(iL,homogen);      // Scaling factor for SVD
    
    double beta=0;
    int incr=1;
    char trans[] = "n";
    
    int Usize = dofn3_f * cutoff->f;
    int Vsize = dofn3_s * cutoff->s;
    int Ksize = 316 * cutoff->f * cutoff->s;
    
    assert((*A)->Nf > 0); /* Cell cannot be empty. */
    // Allocate memory for field values
    (*A)->fieldval = (double *)malloc(dofn3_f*sizeof(double));
    
    // Initialize pointer to the field values of A
    F = (*A)->fieldval;
    
    int matSizeDof;
	if (use_chebyshev)
	    matSizeDof = cutoff_f;
	else
	    matSizeDof = (int)round(pow(2*n-1,3)) *dof->f;
    
    
    double *productfre;
    if (! use_chebyshev) {
        productfre = (double*)calloc(matSizeDof, sizeof(double));
    }
    
    // Obtain the center of the field cell
    fcenter = (*A)->center;
    // Obtain the number of interaction cells
    int ninter = (*A)->iinter;
    
    double *FFCoeff = (double*)malloc(matSizeDof *sizeof(double));
    double* Pf = (double*)malloc(cutoff_f *sizeof(double));
    for (i=0;i<cutoff_f;i++)
		Pf[i] = 0;

    // Compute the field values due to all members of the interaction list
    for (i=0;i<ninter;i++) {
		
        B = (*A)->interaction[i];
        
        // Obtain the center of the source cell
        scenter = B->center;
        cshift = (*A)->cshiftinter[i];
        scenter.x += cshift.x;
        scenter.y += cshift.y;
        scenter.z += cshift.z;
		
        // Note this is the normalized vector
        double R[3];
        R[0] = iL*(scenter.x-fcenter.x);
        R[1] = iL*(scenter.y-fcenter.y);
        R[2] = iL*(scenter.z-fcenter.z);
        
        // Determine the proxy values by pre-multiplying S by V^T
        
        for (j=0; j<matSizeDof; j++) {
		    FFCoeff[j] = 0;
		}
        
        if (!use_chebyshev) {
            Moment2Local(n, R, B->sourcefre, FFCoeff, E +Ksize*shift, Ktable,
                         dof, cutoff, VT +Vsize*shift, Kweights, use_chebyshev);
            
            // Add up the contribution
            for (j=0; j<matSizeDof; j++)
                productfre[j] += FFCoeff[j];
        }else {
        
        Moment2Local(n, R, B->sourceval, FFCoeff, E +Ksize*shift, Ktable,
                     dof, cutoff, VT +Vsize*shift, Kweights, use_chebyshev);
            for (j=0; j<cutoff_f; j++)
                Pf[j] += FFCoeff[j];

        }
    
    }
    free(FFCoeff);
    FFCoeff=NULL;
    
    if (! use_chebyshev) {
        int padSize = round(pow(2*n-1, 3));
        int l3Size = n, l1Size = (int)round(pow(l3Size, 2));
        int l2Pad = 2*n-1, l1Pad = (int)round(pow(l2Pad, 2));
        int count;
        int f, l1, l2, l3;
        
        double res[padSize *dof->f];
        for (f=0; f<dof->f; f++)
            rfftw_one(p_c2r, productfre + f*padSize, res + f*padSize);
        
        free(productfre);
        productfre=NULL;
        
        for (i=count=0; i<n3; i++) {
            l3 = i % l3Size;
            l1 = i / l1Size;
            l2 = i / l3Size % l3Size;
            shift = l1*l1Pad+l2*l2Pad+l3;
            
            for (f=0; f<dof->f; f++, shift+=padSize) {
                Pf[count++] = res[shift]/padSize;
                //printf("Pf[%d] = %f\t", count-1, Pf[count-1]);
            }
        }
    }



    // Compute the field values at the field Chebyshev nodes M2L
    if( use_chebyshev ) {
	    dgemv_(trans,&dofn3_f,&cutoff_f,&scale,U+Usize*shift,&dofn3_f,Pf,&incr,&beta,F,&incr);  // 3c
        
	    // Adjust the field values by the appropriate weight
	    l = 0;
	    for (i=0;i<n3;i++) {
            tmp = Kweights[i];
            for (j=0;j<dof->f;j++) {
                F[l] *= tmp;
                l++;
            }
	    }
	}
	else {
	    for (j=0; j<cutoff_f; j++)
            F[j] = scale *Pf[j];
	} 
	
	free(Pf);Pf = NULL;
    // Recursively compute the kernel interactions for all children cells
    
    // Go to the next level
    //fmmLevel += 1;
    if ((*A)->leaves[0] != NULL) {
        for (i=0;i<8;i++) {
            if ((*A)->leaves[i]->Nf != 0)
                FMMInteraction(&((*A)->leaves[i]), E, Ktable, U, VT,
                               Kweights, n, dof, cutoff, homogen,
                               curTreeLevel+1, use_chebyshev,p_c2r);
            
        }
    }
}

template <typename T>
void H2_3D_Compute<T>::DownwardPass(nodeT **A, vector3 *field, vector3 *source,
                  double *Cweights, double *Tkz, double *q, doft *cutoff,
                  int n, doft *dof, double *phi, int use_chebyshev) {
    /* Add the contributions from the parent cell to each child cell -
	 * otherwise compute all direct interactions and then interpolate
	 * to the field points */
	if ((*A)->leaves[0] != NULL) {
        
	    double *F = (*A)->fieldval;              // Field values for cell
        
	    // Determine which children cells contain field points
        int i, l;
	    for (i=0;i<8;i++) {
            if ((*A)->leaves[i]->Nf != 0) {
                
                double *Fchild = (*A)->leaves[i]->fieldval;
                
                // Manually set the vector from child to parent
                int  idx[3];
                double r[3];
                idx[2] = i%2, idx[1] = i/2%2, idx[0] = i/4%2;
                for (l=0; l<3; l++)
                    if (idx[l] == 0)
                        r[l] = -1;
                    else
                        r[l] =  1;
                
                Local2Local(n, r, F, Fchild, dof, Cweights, use_chebyshev);
                DownwardPass(&((*A)->leaves[i]),field,source,Cweights,Tkz,q,
							 cutoff,n,dof,phi, use_chebyshev);
            }
		}
		
	} else {
        int Nf = (*A)->Nf, Ns, i, j, k, m, *fieldlist = (*A)->fieldlist, nneigh = (*A)->ineigh,
        *sourcelist;
        int l, l1, l2, l3, l4;
        double sum, prefac  = 2.0/(double)n;          // Prefactor for Sn
        if(!use_chebyshev)
            prefac = 1.;
        double prefac3 = prefac*prefac*prefac;   // prefac3 = prefac^3
        double *F = (*A)->fieldval;
        
        double ihalfL = 2./(*A)->length, tmp1, tmp2, *qsource, fieldval[dof->f * Nf];
        vector3 fieldt[Nf], fieldpos[Nf], Sf[n*Nf], fcenter = (*A)->center, cshift,
        *sourcepos;
        nodeT *B;
        
        // Map all field points to the box ([-1 1])^3
        for (i=0;i<Nf;i++) {
            k = fieldlist[i];
            fieldt[i].x = ihalfL*(field[k].x - fcenter.x);
            fieldt[i].y = ihalfL*(field[k].y - fcenter.y);
            fieldt[i].z = ihalfL*(field[k].z - fcenter.z);
        }
        
        // Compute Sf, the mapping function for the field points
        FMMTree->ComputeSn(fieldt,Tkz,n,Nf,Sf, use_chebyshev);
        
        // Obtain the positions of the field points
        for (i=0;i<Nf;i++) {
            k = fieldlist[i];
            fieldpos[i].x = field[k].x;
            fieldpos[i].y = field[k].y;
            fieldpos[i].z = field[k].z;
        }
        
        // Compute the values at the field points
        for (i=0;i<Nf;i++) {
            k = dof->f * fieldlist[i];
            for (l4=0;l4<dof->f;l4++) {
                // Due to far field interactions
                sum = 0;
                l = l4;
                for (l1=0;l1<n;l1++) {
                    tmp1 = Sf[l1*Nf+i].x;
                    for (l2=0;l2<n;l2++) {
                        tmp2 = tmp1*Sf[l2*Nf+i].y;
                        for (l3=0;l3<n;l3++) {
                            sum += F[l]*tmp2*Sf[l3*Nf+i].z;
                            l += dof->f;
                        }
                    }
                }
                phi[k] = prefac3*sum;
                k++;
            }
        }
        
        // Due to near field interactions
        for (m=0;m<nneigh;m++) {
            B = (*A)->neighbors[m];
            sourcelist = B->sourcelist;
            Ns = B->Ns;
            cshift = (*A)->cshiftneigh[m];
            
            sourcepos = (vector3 *)malloc(Ns*sizeof(vector3));
            qsource = (double *)malloc(dof->s*Ns*sizeof(double));
            for (j=0;j<Ns;j++) {
                l = sourcelist[j];
                sourcepos[j].x = source[l].x + cshift.x;
                sourcepos[j].y = source[l].y + cshift.y;
                sourcepos[j].z = source[l].z + cshift.z;
                for (k=0;k<dof->s;k++)
                    qsource[dof->s*j+k] = q[dof->s*l+k];
            }
            
            EvaluateField(fieldpos,sourcepos,qsource,Nf,Ns,dof,fieldval);
            
            for (i=0;i<Nf;i++) {
                j = dof->f * fieldlist[i];
                l = dof->f * i;
                for (k=0;k<dof->f;k++) {
                    phi[j+k] += fieldval[l+k];
                }
            }
            
            free(sourcepos);
            sourcepos = NULL;
            free(qsource);
            qsource = NULL;
        }
	}
}

/*
 * Function: EvaluateField
 * -------------------------------------------------------------------
 * Evaluates the field due to interactions between a pair of cells.
 * P2P kernel.
 */
template <typename T>
void H2_3D_Compute<T>::EvaluateField(vector3 *field, vector3 *source, double *q, int Nf,int Ns, doft *dof, double *fieldval) {
	/* Compute K*q
	 * K is (dof->f*Nf) x (dof->s*Ns)
	 * q is a vector of length dof->s*Ns
	 * q stores q[Ns][dof->s] in row major
	 * so the result fieldval[Nf][dof->f] is also in row major
	 */
    
	int i, j, k, l, count, count_kernel;
	int dof2 = dof->f * dof->s;
	int LDA  = dof->f * Nf;
	int N    = dof->s * Ns;
    
	// Kij is the point to point kernel
	// Kcell stores all the Kij's
	
	double *Kij, *Kcell;
	Kij   = (double*) malloc(dof2 * sizeof(double));
	Kcell = (double*) malloc(LDA * N * sizeof(double));
    
	for (j=0;j<Ns;j++) {
		for (i=0;i<Nf;i++) {
			if (source[j].x != field[i].x || source[j].y != field[i].y || source[j].z != field[i].z) {
                FMMTree->EvaluateKernel(field[i],source[j],Kij,dof);
            }
			else { // The source point and the filed overlape
                FMMTree->EvaluateKernel(field[i],source[j],Kij,dof);
                for (k=0;k<dof2;k++) {
                    if (isinf(Kij[k])) {
                        Kij[k] = 0;
                    }
                }
            }
  
			count_kernel = dof->f * i + LDA * dof->s * j;
			count = 0;
			for (k=0;k<dof->s;k++)
				for (l=0;l<dof->f;l++, count++)
					// Kcell and Kij are column-major storage
					Kcell[count_kernel + k * LDA + l] = Kij[count];
		}
	}
	
	free(Kij);
    Kij = NULL;
	// Do matrix vector product
	double alpha = 1, beta = 0;
    int incr = 1;
	char trans[] = "n";
	dgemv_(trans, &LDA, &N, &alpha, Kcell, &LDA, q, &incr, &beta, fieldval, &incr);
    
	free(Kcell);
    Kcell = NULL;
}


/* Function: L2L
 * Input:
 *       n: Number of Chebyshev nodes
 *       r: (r[0], r[1], r[2]) is the vector from the child box to the parent box
 *       F: Field values of the parent
 *
 * Output:
 *  Fchild: Field values of the child
 *
 */
template <typename T>
void H2_3D_Compute<T>::Local2Local(int n, double *r, double *F, double *Fchild, doft *dof, double *Cweights, int use_chebyshev){
    
    int l, l1, l2, l3, l4, count1, count2, count3, wstart;
	
	int n2 = n*n;                       // n2 = n^2
	int n3 = n*n*n;                     // n3 = n^3
	int dofn  = dof->f * n;
	int dofn2 = dof->f * n2;
	int dofn3 = dof->f * n3;
	
	double prefac  = 2.0/(double)n;          // Prefactor for Sn
    if(!use_chebyshev)
        prefac = 1.;
	double prefac3 = prefac*prefac*prefac;   // prefac3 = prefac^3
	
	double Fx[dofn3], Fy[dofn3];
    
    
    // Recover the index of the child box
    int idx[3] = {0};
    if (  r[0] > 0)
        idx[0] = 1;
    if (  r[1] > 0)
        idx[1] = 1;
    if (  r[2] > 0)
        idx[2] = 1;
    
	// Interpolate the parent field along the x-component
    l = 0;
	if (idx[0] == 0)
	    wstart = 0;
	else
	    wstart = n2;
    
    for (l1=0;l1<n;l1++) {
		count2 = wstart + l1;
		for (l2=0;l2<n;l2++) {
            for (l3=0;l3<n;l3++) {
                count1 = l2*n + l3;
                for (l4=0;l4<dof->f;l4++) {
		            count3 = dof->f*count1 + l4;
                    Fx[l]  = ddot_(&n,&F[count3],&dofn2,&Cweights[count2],&n);
                    l++;
                }
		    }
		}
    }
    
    // Interpolate the parent field along the y-component
    l = 0;
    if (idx[1] == 0)
        wstart = 0;
    else
        wstart = n2;
    
    for (l1=0;l1<n;l1++) {
        for (l2=0;l2<n;l2++) {
            count2 = wstart + l2;
            for (l3=0;l3<n;l3++) {
                count1 = l1*n2 + l3;
                for (l4=0;l4<dof->f;l4++) {
                    count3 = dof->f*count1 + l4;
                    Fy[l] = ddot_(&n,&Fx[count3],&dofn,&Cweights[count2],&n);
                    l++;
                }
            }
        }
    }
    
    /* Interpolate the parent field along the z-component and add
     to child field */
    l = 0;
    if (idx[2] == 0)
        wstart = 0;
    else
        wstart = n2;
    
    for (l1=0;l1<n2;l1++) {
        count1 = l1*n;
        for (l3=0;l3<n;l3++) {
            count2 = wstart + l3;
            for (l4=0;l4<dof->f;l4++) {
                count3 = dof->f*count1 + l4;
				Fchild[l] += prefac3*ddot_(&n,&Fy[count3],&(dof->f), &Cweights[count2],&n);
				l++;
            }
        }
    }
    
}


template <typename T>
void H2_3D_Compute<T>::Moment2Local(int n, double *R, double *cell_mpCoeff, double *FFCoeff,
                                    double *E, int *Ktable, doft *dof, doft *cutoff, double *VT,
                                    double *Kweights, int use_chebyshev) {
    
	int n3 = n*n*n, cutoff_f = cutoff->f, cutoff_s = cutoff->s;
    int cutoff2  = cutoff_f * cutoff_s;
	int dofn3    = dof->s * n*n*n;
    
    // Final multipole expansion (Omega w) =  Sw ; v * Sw = Wl
    double Sw[dofn3];
    
    if(use_chebyshev) {
        int l = 0, l1, l2;
        double tmp;
        for (l1=0;l1<n3;l1++) {
            tmp = Kweights[l1];
            for (l2=0;l2<dof->s;l2++) {
                Sw[l] = cell_mpCoeff[l] * tmp;
                l++;
            }
        }
    }
    
    // 3a in page 8718
    int incr = 1;
    double  alpha = 1, beta = 0;
    char trans = 'n';
    
    /* TODO: The current code computes 'CompCoeff' each times for the same 'cell_mpCoeff' when
     *       the cell is in the interaction list of multiple cells. And one optimization is to
     *       overwrite 'cell_mpCoeff' with 'CompCoeff' and set the flag, so the next time
     *       'CompCoeff' can be used directly
     */
    double CompCoeff[cutoff_s];
    if(use_chebyshev) {
        dgemv_(&trans,&cutoff_s,&dofn3,&alpha,VT,&cutoff_s,Sw,&incr,&beta,CompCoeff,&incr);
    }

    
    // Determine the corresponding index in the lookup table
    int k1 = (int)round(R[0]) + 3;
    int k2 = (int)round(R[1]) + 3;
    int k3 = (int)round(R[2]) + 3;
    int ninteract = Ktable[49*k1+7*k2+k3];
    int count = ninteract*cutoff2;
    
    if(!use_chebyshev)
        count = ninteract*(2*n-1)*(2*n-1)*(2*n-1)*dof->s*dof->f;
    
    // Compute the field proxy values: Matrix Vector multiplication in BLAS: dgemv - 3b
    // cutoff: cutoff on the number of SVD values,
    // Ecell : M2L operator - E[count] : precalculated, P input, M expansion Output : Pf L expansion
    // CHANGE0731: no 'Ecell' used
    if(use_chebyshev) {
        dgemv_(&trans, &cutoff_f, &cutoff_s, &alpha, E+count, &cutoff_f,
               CompCoeff, &incr, &beta,FFCoeff,&incr); // 3b  Ecell is the kernel
    }
    else {
        
        // entry-wise product of cell_mpCoeff
        int N = (int)round(pow(2*n-1, 3));
        int f, s, shift1, shift2=0, shift3=0;
        for (f=0; f<dof->f; f++, shift2+=N)
            for (s=shift1=0; s<dof->s; s++) {
                FrequencyProduct(N, &E[count+shift3], &cell_mpCoeff[shift1],
                                 &FFCoeff[shift2]);
                shift1 += N;
                shift3 += N;
            }
    }
    
}

/*
 * Compute the entry-wise product of two frequencies from rfftw
 * Note 'res' has been initialized
 */
template <typename T>
void H2_3D_Compute<T>::FrequencyProduct(int N, double *Afre, double *xfre, double *res) {
    
    int i;
    res[0] += Afre[0]*xfre[0];
    for (i=1; i<N; i++) {
        if (i<(N+1)/2)
            res[i] += Afre[i]*xfre[i] + Afre[N-i]*xfre[N-i];
        else
            res[i] += Afre[N-i]*xfre[i] - Afre[i]*xfre[N-i];
    }
    if (N%2 == 0)
        res[N/2] += Afre[N/2]*xfre[N/2];

}


/* Function: M2M
 * Input:
 *       n: Number of Chebyshev nodes
 *       r: (r[0], r[1], r[2]) is the vector from the child box to the parent box
 *  Schild: Source values of the child
 *
 * Output:
 *  SS: Source values of the parent
 *
 */
template <typename T>
void H2_3D_Compute<T>::Moment2Moment(int n, double *r, double *Schild, double *SS, doft *dof, double *Cweights) {
    
	int l, l1, l2, l3, l4;
	int count1, count2, count3, wstart;
	
	int incr = 1;
	int n2 = n*n;                       // n2 = n^2
	int n3 = n*n*n;                     // n3 = n^3
	
	int dofn  = dof->s * n;
	int dofn2 = dof->s * n2;
	int dofn3 = dof->s * n3;
	
	double Sy[dofn3], Sz[dofn3];
    
    
    // Recover the index of the child box
    int idx[3] = {0};
    if (  r[0] > 0)
        idx[0] = 1;
    if (  r[1] > 0)
        idx[1] = 1;
    if (  r[2] > 0)
        idx[2] = 1;
    
    // Gather the children source along the z-component
    if (idx[2] == 0)
	    wstart =  0;
    else
	    wstart =  n2;
    
    l = 0;
    for (l1=0;l1<n2;l1++) {
	    count1 = l1*n;
	    for (l3=0;l3<n;l3++) {
            count2 = wstart + l3*n;
	        for (l4=0;l4<dof->s;l4++) {
                count3 = dof->s*count1 + l4;
                Sz[l]  = ddot_(&n,&Schild[count3],&(dof->s),&Cweights[count2],&incr);
                l++;
	        }
	    }
    }
    
    // Gather the children sources along the y-component
    if (idx[1] == 0)
        wstart =  0;
    else
        wstart =  n2;
    
    l = 0;
    for (l1=0;l1<n;l1++) {
        for (l2=0;l2<n;l2++) {
            count2 = wstart + l2*n;
            for (l3=0;l3<n;l3++) {
                count1 = l1*n2 + l3;
                for (l4=0;l4<dof->s;l4++) {
                    count3 = dof->s*count1 + l4;
                    Sy[l]  = ddot_(&n,&Sz[count3],&dofn,&Cweights[count2],&incr);
                    l++;
                }
            }
        }
    }
    
    /* Gather the children sources along the x-component and determine
     the parent sources */
    l = 0;
    if (idx[0] == 0)
        wstart =  0;
    else
        wstart =  n2;
    
    for (l1=0;l1<n;l1++) {
        count2 = wstart + l1*n;
        for (l2=0;l2<n;l2++) {
            for (l3=0;l3<n;l3++) {
                count1 = l2*n + l3;
                for (l4=0;l4<dof->s;l4++) {
                    count3  = dof->s*count1 + l4;
                    SS[l] = ddot_(&n,&Sy[count3],&dofn2,&Cweights[count2],&incr);
                    l++;
                }
            }
        }
    }
}

/*
 * Function: InteractionList
 * -------------------------------------------------------------------
 * For each node at each level of the octree, the interaction list
 * is determined.
 */
template <typename T>
void H2_3D_Compute<T>::InteractionList(nodeT **A, int levels) {
	int i, j, k, ineigh, iinter, nneigh;
	nodeT *B, *C;
	vector3 center1, center2, cshift, diff;
	double dist, cutoff;
	
	assert((*A)->Nf > 0);
	
	if (levels > 0) {
		nneigh = (*A)->ineigh;
		
		/* Sets the cutoff between near and far to be L (this is equivalent
		 * to a one cell buffer) */
		cutoff = (*A)->length;
		
		/*
		 * Finds all neighbors that are too close for the far field
		 * approximation and stores them in the neighbors array -
		 * Also finds the neighboring cells that are sufficiently
		 * far away and stores them in interaction array
		 */
		for (i=0;i<nneigh;i++) {
			B = (*A)->neighbors[i]; // All the neighbors of A
			cshift = (*A)->cshiftneigh[i];
			for (j=0;j<8;j++) {
				if (B->leaves[j]->Ns != 0) { // All the children j of the neighbor cluster A
					/* Skip empty neighbor clusters. */
					
					center1 = B->leaves[j]->center;
					center1.x += cshift.x;
					center1.y += cshift.y;
					center1.z += cshift.z;
					
					for (k=0;k<8;k++) {
						C = (*A)->leaves[k]; // Child k of cluster A
						
						if (C->Nf != 0) {
							/* Skip empty children clusters. */
							
							ineigh = C->ineigh;
							iinter = C->iinter;
							center2 = C->center;
							
							diff.x = center1.x - center2.x;
							diff.y = center1.y - center2.y;
							diff.z = center1.z - center2.z;
							dist = sqrt(diff.x*diff.x+diff.y*diff.y+diff.z*diff.z);
							
							if (dist < cutoff) {
								C->neighbors[ineigh] = B->leaves[j]; // This is a neighbor
								C->cshiftneigh[ineigh] = cshift;
								C->ineigh++;
							} else {
								C->interaction[iinter] = B->leaves[j]; // This is a well-separated cluster
								C->cshiftinter[iinter] = cshift;
								C->iinter++;
							}
						}
						
					}
				}
			}
		}
		
        // recursively build the interaction lists
		if ((*A)->leaves[0] != NULL) {
			for (k=0;k<8;k++) {
				if ((*A)->leaves[k]->Nf != 0) { /* Skip empty clusters */
					InteractionList(&((*A)->leaves[k]),levels-1);
				}
			}
		}
	}
}


/*
 * Function: DistributeSources
 * -------------------------------------------------------------------
 * Distributes all of the sources to the appropriate cell at the
 * bottom of the octree.
 */
template <typename T>
void H2_3D_Compute<T>::DistributeSources(nodeT **A, vector3 *source, int *sourcelist,
                                         int levels) {
	int i, j, k, m, z;
	int Ns = (*A)->Ns;
	vector3 point, center;
	int *sourcecell[8], *S;
	for (z = 0; z < 8; z++)
		sourcecell[z] = (int *)malloc(Ns * sizeof(int));
	
	if (levels > 0) {
        // Obtain the center of the cell
		center = (*A)->center;
		
        // Distribute the field points
		if (Ns != 0) {
            // Determine which child cell each field and source point belongs to
			for (i=0;i<Ns;i++) {
				k = sourcelist[i];
				point = source[k];
				
                // Determine which cell the point is in
				if (point.x < center.x) {
					if (point.y < center.y) {
						if (point.z < center.z)
							j = 0;
						else
							j = 1;
					} else {
						if (point.z < center.z)
							j = 2;
						else
							j = 3;
					}
				} else {
					if (point.y < center.y) {
						if (point.z < center.z)
							j = 4;
						else
							j = 5;
					} else {
						if (point.z < center.z)
							j = 6;
						else
							j = 7;
					}
				}
				
                // Add the field point to the list for child cell j
				m = (*A)->leaves[j]->Ns;
				sourcecell[j][m] = k;
				(*A)->leaves[j]->Ns++;
			}
			
            // Recursively distribute the points
			for (j=0;j<8;j++) {
				S = sourcecell[j];
				DistributeSources(&((*A)->leaves[j]),source,S,levels-1);
			}
		}
	} else if (levels == 0) {
		Ns = (*A)->Ns;
		(*A)->sourcelist = (int *)malloc(Ns*sizeof(int));
		S = (*A)->sourcelist;
		for (i=0;i<Ns;i++)
			S[i] = sourcelist[i];
	}
	for (z = 0; z < 8; z++) {
		free(sourcecell[z]);
        sourcecell[z] = NULL;
    }
}


/*
 * Function: DistributeFieldPoints
 * -------------------------------------------------------------------
 * Distributes all of the field points to the appropriate cell at the
 * bottom of the octree.
 */
template <typename T>
void H2_3D_Compute<T>::DistributeFieldPoints(nodeT **A, vector3 *field, int *fieldlist,
                                             int levels) {
	int i, j, k, m, z;
	int Nf = (*A)->Nf;
	vector3 point, center;
	int *fieldcell[8], *F;
	for (z = 0; z < 8; z++)
		fieldcell[z] = (int *)malloc(Nf * sizeof(int));
	
	if (levels > 0) {
        // Obtain the center of the cell
		center = (*A)->center;
		
        // Distribute the field points
		if (Nf != 0) {
            // Determine which child cell each field and source point belongs to
			for (i=0;i<Nf;i++) {
				k = fieldlist[i];
				point = field[k];
				
                // Determine which cell the point is in
				if (point.x < center.x) {
					if (point.y < center.y) {
						if (point.z < center.z)
							j = 0;
						else
							j = 1;
					} else {
						if (point.z < center.z)
							j = 2;
						else
							j = 3;
					}
				} else {
					if (point.y < center.y) {
						if (point.z < center.z)
							j = 4;
						else
							j = 5;
					} else {
						if (point.z < center.z)
							j = 6;
						else
							j = 7;
					}
				}
				
                // Add the field point to the list for child cell j
				m = (*A)->leaves[j]->Nf;
				fieldcell[j][m] = k;
				(*A)->leaves[j]->Nf++;
			}
			
            // Recursively distribute the points
			for (j=0;j<8;j++) {
				F = fieldcell[j];
				DistributeFieldPoints(&((*A)->leaves[j]),field,F,levels-1);
			}
		}
	} else if (levels == 0) {
		Nf = (*A)->Nf;
		(*A)->fieldlist = (int *)malloc(Nf*sizeof(int));
		F = (*A)->fieldlist;
		for (i=0;i<Nf;i++)
			F[i] = fieldlist[i];
	}
	for (z = 0; z < 8; z++) {
		free(fieldcell[z]);
        fieldcell[z] = NULL;
    }
}



template <typename T>
H2_3D_Compute<T>::~H2_3D_Compute() {
    if (FMMTree->tree != NULL) {
        FMMTree->FreeNode(FMMTree->tree);
        FMMTree->tree = NULL;
    }
    if (FMMTree->Kweights!= NULL) {
        free(FMMTree->Kweights);
        FMMTree->Kweights = NULL;
    }
    if (FMMTree->Cweights!= NULL) {
        free(FMMTree->Cweights);
        FMMTree->Cweights = NULL;
    }
    if (FMMTree->Tkz!= NULL) {
        free(FMMTree->Tkz);
        FMMTree->Tkz = NULL;
    }
    if (FMMTree->K!= NULL) {
        free(FMMTree->K);
        FMMTree->K = NULL;
    }
    if (FMMTree->U!= NULL) {
        free(FMMTree->U);
        FMMTree->U = NULL;
    }
    if (FMMTree->VT!= NULL) {
        free(FMMTree->VT);
        FMMTree->VT = NULL;
    }
}



#endif //(__compute_hpp__)
