import numpy as np
from mpi4py import MPI
from sympy import fwht
import time
from numpy.random import seed, random, choice
from math import sqrt


def sketch_left(A: np.array, l: int, comm: MPI.Comm, seedLocal: int = 166297, seedGlobal: int = 1263):
    '''
    Sketch matrix A from left with block hrhrt, computes Omega*A, A is row-block wise
    distributed of size mxn, each block of rows is of size rxn, Omega is of size lxm.
    :param A: matrix to sketch of size mxn
    :param l: sketching size
    :param comm: mpi communicator to allow parallelism
    :param seedLocal: seed to use when getting the local diagonal Rademacher
    :param seedGlobal: seed to use to get the random columns, shared among processors
    '''
    m,n = A.shape()
    p = comm.Get_size()
    rank = comm.Get_rank()
    r = int(m/p)
    t0 = time.perf_counter()
    seedGlobal = comm.bcast(seedGlobal, root = 0)
    

