"""
    Randomized algorithm for Hermitian eigenvalue problems using onepass method
    Computs k largest eigenvalues computed using the randomized algorithm

    Parameters:
    k	= int, 
            number of eigenvalues/vectors to be estimated
    p   = int,
            oversampling parameter which can improve accuracy of resulting solution
    FMM parameters.. 
    
    Returns:

    w	= double, k
            eigenvalues
    u 	= n x k 
            eigenvectors - commented 

    References:
    -----------
    .. [1] Halko, Nathan, Per-Gunnar Martinsson, and Joel A. Tropp. "Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions." SIAM review 53.2 (2011): 217-288.
    ---------
"""    
import matplotlib
matplotlib.use("Agg")
from FMMTree import *
from FMMCompute import *
import numpy as np
import timeit
from scipy.linalg import qr, inv, pinv,svd, pinv2, eigh
import matplotlib.pyplot as plt	
import random
import math
from scipy.spatial import distance_matrix
#npts = np.array([10,22,46,100,171]) # 1,000, 9,261, 97,336, 1000,000, 5,000,000


def check_accuracy(myTree, num_rows, nCols, Ns, Nf, target, source, q, stress):
    stress_exact = [0] * num_rows * nCols;
    for k in range(nCols):
        for i in range(num_rows):
            for j in range(Ns): 
                val = myTree.EvaluateKernel(target[i], source[j]);
                if math.isnan(val) or math.isinf(val):
                    val = 0
                stress_exact[k*num_rows+i] += val * q[k*Nf+j];
    diff, sum_Ax = 0, 0
    for k in range(nCols):
        for i in range(num_rows): 
            diff += (stress[k*Nf+i] - stress_exact[k*num_rows+i])*(stress[k*Nf+i] - stress_exact[k*num_rows+i]); 
            sum_Ax += stress_exact[k*num_rows+i] * stress_exact[k*num_rows+i]
    print("diff of first 10 values = %f"% (math.sqrt(diff) / math.sqrt(sum_Ax)))



npts = np.array([10000, 80000])
tree_levels = [3,4]

nCols = 100 # number of eigenvalues
nOverSample =  20 # oversample for randomized method
nTotalCols = nCols + nOverSample
# set FMM parameters
L, interpolation_order, eps, use_chebyshev = 1, 4, 1e-5, 1
print("********* Onepass Randomized method starts *********\n")
print("- computes %d largest eigenvalues of a covariance matrix with oversample %d" % (nCols,nOverSample))

random.seed(1)

for i, Ns  in enumerate(npts):
    data = np.random.random([Ns, 3])
    tree_level = tree_levels[i]
    Nf = Ns
    tStart = timeit.default_timer()
    tBegin = tStart 

    print("(1) define target and source")    


    print("L                            : %3.3f"%L)
    print("interpolation order          : %d"%interpolation_order)
    print("Ns                           : %d"%Ns)
    print("nTotalCols(#cols in weights) : %d"% nTotalCols)
    print("tree tree_level              : %d"%tree_level)
    print("eps                          : %3.3e"%eps)
    print("\n")


    print("(2) Generate gaussian random matrix and copy it to q - need to be updated")

    t_fmm = 0
    t_rest = 0

    # Generate gaussian random matrix 
    Omega = np.random.randn(Ns,nTotalCols)

    # Convert python to C
    target_c = vecOfvec3()
    source_c = vecOfvec3()
    weight_c = vecOfdouble()
    output_c = vecOfdouble()
    convert_to_vecOfvec3(data, source_c)
    convert_to_vecOfvec3(data, target_c)
    convert_to_vecOfdouble(Omega.T.reshape(Ns*nTotalCols), weight_c)

    output_c[:] = np.zeros(Ns*nTotalCols)

    # Build FMM Tree
    print("(3) Building FMM Tree")
    tStart = timeit.default_timer()


    myTree = myKernel(L, tree_level, interpolation_order, eps, use_chebyshev)
    myTree.buildFMMTree()

    Compute(myTree, target_c, source_c, weight_c, nTotalCols, output_c)

    t_fmm = t_fmm + timeit.default_timer() - tStart
	
    tStart = timeit.default_timer()

    Y = np.empty([Ns*nTotalCols], dtype=np.float64)
    convert_to_numpy(output_c, Y)
    Y = Y.reshape(nTotalCols, Ns).T


    Q,_ = qr(Y, mode = 'economic')

    convert_to_vecOfdouble(Q.T.reshape(Ns*nTotalCols), weight_c)
    output_c[:] = np.zeros(Ns*nTotalCols)

    t_rest = t_rest + timeit.default_timer() - tStart


    tStart = timeit.default_timer()

    myTree = myKernel(L, tree_level, interpolation_order, eps, use_chebyshev)
    myTree.buildFMMTree()
    Compute(myTree, target_c, source_c, weight_c, nTotalCols, output_c)
    t_fmm = t_fmm + timeit.default_timer() - tStart

    tStart = timeit.default_timer()

    Y1 = np.empty([Ns*nTotalCols], dtype=np.float64)
    convert_to_numpy(output_c, Y1)
    Y1 = Y1.reshape(nTotalCols, Ns).T

    B  = np.dot(Q.T, Y1)
  
	#Eigen subproblem
    w, v = eigh(B)
    
    w = w[::-1]
    w = w[:nCols]
    t_rest = t_rest + timeit.default_timer() - tStart

    tEnd = timeit.default_timer()
    print "t_fmm", t_fmm, "t_rset", t_rest, "t_tot", t_fmm + t_rest
    print("-----------------------------------------------------------")
    print("** For %d by %d cov. matrix, total elapsed time is %3.3f sec" % (Ns,Ns,tEnd - tBegin))



    # Exact 
    K = distance_matrix(data, data)
    K = np.exp(-K, out=K) # need to change 

    t_mat = 0
    t_res = 0

    tStart = timeit.default_timer()

    Y = np.dot(K, Omega)
    t_mat = t_mat + timeit.default_timer() - tStart

    tStart = timeit.default_timer()
    q,_ = qr(Y, mode = 'economic')
    t_res = t_res + timeit.default_timer() - tStart
    tStart = timeit.default_timer()

    Y = np.dot(K, q)
    t_mat = t_mat + timeit.default_timer() - tStart
    tStart = timeit.default_timer()

    B  = np.dot(q.T, Y) 
    w_ext, v = eigh(B)
    w_ext = w_ext[::-1]
    w_ext = w_ext[:nCols]
    t_res = t_res + timeit.default_timer() - tStart
    print "t_mat", t_mat, "t_res", t_res, "t_tot", t_mat + t_res
 

    print np.linalg.norm(w - w_ext) / np.linalg.norm(w_ext)


    # plot eigen decay
    fig = plt.figure()
    plt.semilogy(range(nCols), w, 'k.-')
    plt.title("eigen values for a %d by %d covariance matrix" % (Ns,Ns))
    plt.axis('tight')
    fname = "eigv_%d.png" % (Ns)
    fig.savefig(fname, dpi=fig.dpi)
    plt.close(fig)

    print("** eig decay plot saved in %s" % (fname))
    print("-----------------------------------------------------------\n\n")

if __name__ == '__main__':
    exit(main())
    