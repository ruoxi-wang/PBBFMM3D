#!/usr/bin/env python2.7
import math
import random
from FMMTree import *
from FMMCompute import *
import numpy as np
import timeit

def load_data(file_xcoord, file_ycoord, file_zcoord):
    xcoord = np.loadtxt(file_xcoord)
    ycoord = np.loadtxt(file_ycoord)
    zcoord = np.loadtxt(file_zcoord)
    nx = xcoord.shape[0]
    ny = ycoord.shape[0]
    nz = zcoord.shape[0]
    print("nx = %d, ny = %d, nz = %d"% (nx, ny, nz))
    Ns = Ny = nx*ny*nz
   
    source = np.empty([Ns, 3], dtype=np.float64)
    target = np.empty([Ny, 3], dtype=np.float64)
    index = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                source[index,:] = [xcoord[i], ycoord[j], zcoord[k]]
                target[index,:] = [xcoord[i], ycoord[j], zcoord[k]]
                index = index + 1

    return source, target



def check_accuracy(myTree, num_rows, nCols, Ns, Nf, target, source, weight, output):
    output_exact = [0] * num_rows * nCols;
    for k in range(nCols):
        for i in range(num_rows):
            for j in range(Ns): 
                val = myTree.EvaluateKernel(target[i], source[j]);
                output_exact[k*num_rows+i] += val * weight[k*Nf+j];
    diff, sum_Ax = 0, 0
    for k in range(nCols):
        for i in range(num_rows): 
            diff += (output[k*Nf+i] - output_exact[k*num_rows+i])*(output[k*Nf+i] - output_exact[k*num_rows+i]); 
            sum_Ax += output_exact[k*num_rows+i] * output_exact[k*num_rows+i]
    print("diff of first 10 values = %f"% (math.sqrt(diff) / math.sqrt(sum_Ax)))


def main():

    print("Loading data...")
    source, target = load_data("../input/xcoord.txt", "../input/ycoord.txt", "../input/zcoord.txt")
    Ns, Nf, nCols = len(source), len(target), 1
    weight = np.ones(Ns*nCols)

    print("Ns = %d, Nf = %d"% (Ns, Nf))

    # set FMM parameters
    L, tree_level, interpolation_order, eps, use_chebyshev = 90, 6, 5, 1e-5, 0

    print("L                                     : %3.3f"%L)
    print("interpolation_order (# chebyshev)     : %d"%interpolation_order)
    print("Ns (=Nf)                              : %d"%Ns)
    print("Nf (=Nf)                              : %d"%Nf)
    print("nCols(# of cols in weight)            : %d"%nCols)
    print("tree_level                            : %d"%tree_level)
    print("eps                                   : %3.3e"%eps)


    # Convert python to c
    target_c = vecOfvec3()
    source_c = vecOfvec3()
    weight_c = vecOfdouble()
    convert_to_vecOfdouble(weight, weight_c)
    convert_to_vecOfvec3(source, source_c)
    convert_to_vecOfvec3(source, target_c)

    # Build FMM Tree
    print("Building FMM Tree...")
    start = timeit.default_timer()
    myTree = myKernel(L, tree_level, interpolation_order, eps, use_chebyshev)
    myTree.buildFMMTree()
    tPre = timeit.default_timer() - start


    output_c = vecOfdouble()
    output_c[:] = np.zeros(Ns)
    # FMM compute
    print("Computing...")
    start = timeit.default_timer()
    Compute(myTree, target_c, source_c, weight_c, nCols, output_c)
    tFMM = timeit.default_timer() - start

    print("Pre-computation time: %3.3f"%tPre)
    print("FMM computing time:   %3.3f"%tFMM)
    print("FMM total time:       %3.3f"%(tPre+tFMM))


    # Convert c to python
    output = np.empty(Ns, dtype=np.float64)
    convert_to_numpy(output_c, output)


    # Checking accuracy
    num_rows = 10
    check_accuracy(myTree, num_rows, nCols, Ns, Nf, target_c, source_c, weight_c, output_c)




if __name__ == '__main__':
    exit(main())



