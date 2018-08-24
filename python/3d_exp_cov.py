#!/usr/bin/env python2.7
import math
import random
from FMMTree import *
from FMMCompute import *
import numpy as np
import timeit

def load_data(file_xcoord, file_ycoord, file_zcoord, source, field, data_size, nCols, q):
    xcoord = np.loadtxt(file_xcoord)
    ycoord = np.loadtxt(file_ycoord)
    zcoord = np.loadtxt(file_zcoord)
    nx = xcoord.shape[0]
    ny = ycoord.shape[0]
    nz = zcoord.shape[0]
    print("nx = %d, ny = %d, nz = %d"% (nx, ny, nz))
    data_size[0] = data_size[1] = nx*ny*nz
    Ns = Ny = nx*ny*nz
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                source.append(vector3(xcoord[i], ycoord[j], zcoord[k]))
                field.append(vector3(xcoord[i], ycoord[j], zcoord[k]))
 
    for l in range(nCols):
        for i in range(Ns):
            q.append(1)

def check_accuracy(myTree, num_rows, nCols, Ns, Nf, field, source, q, stress):
    stress_exact = [0] * num_rows * nCols;
    for k in range(nCols):
        for i in range(num_rows):
            for j in range(Ns): 
                val = myTree.EvaluateKernel(field[i], source[j]);
                stress_exact[k*num_rows+i] += val * q[k*Nf+j];
    diff, sum_Ax = 0, 0
    for k in range(nCols):
        for i in range(num_rows): 
            diff += (stress[k*Nf+i] - stress_exact[k*num_rows+i])*(stress[k*Nf+i] - stress_exact[k*num_rows+i]); 
            sum_Ax += stress_exact[k*num_rows+i] * stress_exact[k*num_rows+i]
    print("diff of first 10 values = %f"% (math.sqrt(diff) / math.sqrt(sum_Ax)))


def main():
    field = vector_vector3()
    source = vector_vector3()
    q = vector_double()
    data_size = [None, None]
    nCols = 1
    print("Loading data...")
    load_data("../input/xcoord.txt", "../input/ycoord.txt", "../input/zcoord.txt", source, field, data_size, nCols, q)
    Ns, Nf = data_size

    print("Ns = %d, Nf = %d"% (Ns, Nf))

    # set FMM parameters
    L, level, nChebNode, eps, use_chebyshev = 90, 6, 5, 1e-5, 0

    print("L                    : %3.3f"%L)
    print("n (# chebyshev)      : %d"%nChebNode)
    print("Ns (=Nf)             : %d"%Ns)
    print("Nf (=Nf)             : %d"%Nf)
    print("nCols(# of cols in q): %d"%nCols)
    print("level                : %d"%level)
    print("eps                  : %3.3e"%eps)

    # Build FMM Tree
    print("Building FMM Tree...")
    start = timeit.default_timer()
    myTree = myKernel(L, level, nChebNode, eps, use_chebyshev)
    myTree.buildFMMTree()
    tPre = timeit.default_timer() - start

    # FMM compute
    stress = vector_double()
    for i in range(Ns):
        stress.append(0)

    print("Computing...")
    start = timeit.default_timer()
    Compute(myTree, field, source, Ns, Nf, q, nCols, stress)
    tFMM = timeit.default_timer() - start

    # Checking accuracy
    num_rows = 10
    check_accuracy(myTree, num_rows, nCols, Ns, Nf, field, source, q, stress)


    print("Pre-computation time: %3.3f"%tPre)
    print("FMM computing time:   %3.3f"%tFMM)
    print("FMM total time:       %3.3f"%(tPre+tFMM))

if __name__ == '__main__':
    exit(main())



