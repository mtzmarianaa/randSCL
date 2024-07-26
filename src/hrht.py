import numpy as np
from mpi4py import MPI
from sympy import fwht
import time
from numpy.random import seed, random, choice
from math import sqrt
import pdb
import scipy
from numpy.linalg import svd

def hrht(Ai: np.array,l: int, r:int, s: int, seedLocal: int = 166297, seedGlobal: int = 1263):
    '''
    Sketches Ai using a hashed randomized Hadamard transform,
    can be block wise. Returns B = Omega x Ai,
    Omega is lxr and B is lxn, Ai is rxn.
    :param Ai: numpy array, (sub)matrix to sketch
    :param l: sketch dimension
    :param r: number of rows in Ai
    :param s: parameter for s hashing
    :param seedLocal: seed to use when getting the local diagonal Rademacher
    :param seedGlobal: seed to use to get the random columns, shared among processors
    '''
    m, n = Ai.shape
    t0 = time.perf_counter()
    B = np.zeros((l, n))
    # Beging sketching
    seed(seedLocal)
    Di = np.array([1 if random() < 0.5 else -1 for i in range(r)])
    Ai = (1/sqrt(r))*np.array(fwht(np.diag(Di)@Ai), dtype=np.float64) # scaled randomized Hadamard transform
    seed(seedGlobal) # Seed to generate the same hashing matrix
    S = np.array( [ np.sort(choice(range(l), s, replace = False)) for i in range(r)   ] ).transpose()
    t1 = time.perf_counter()
    for j in range(r):
        seed(seedLocal + j)
        randN = [1 if random() < 0.5 else -1]
        B[ S[:, j], :] += (1/sqrt(s))*randN[0]*Ai[j, :]
    t2 = time.perf_counter()
    return B, t1-t0, t2-t0

def blockHRHT_serial(A: np.array, l: int, s: int, p: int, seedGlobal: int = 1263):
    '''
    Serial implementation of block HRHT.
    :param A: numpy array, mxn
    :param l: sketch dimension
    :param s: parameter for s hashing
    :param p: number of (fake) processors to use
    :param seedGlobal: seed to use to get the random columns, shared among processors
    '''
    m, n = A.shape
    assert( m % p == 0 ) # Has to be divisible to distribute this
    B = np.zeros((l, n))
    # Do block partitions
    r = int(m/p)
    tF = 0
    tL = 0
    for i in range(p):
        Bi, tf, tl = hrht(A[ r*i:(r*(i+1)), : ], l = l, r = r, s = s, seedGlobal = seedGlobal, seedLocal = 16*i)
        tF += tf
        tL += tl
        B += Bi
    return (1/sqrt(p))*B, tF, tL


def rand_nystrom_svd(A, rank, Omega = None, l = None, s = None, p = None, seedGlobal = 1263):
    """
    Compute the randomized Nyström rank k approximation given the sketching
    matrix Omega (uses SVD instead of Cholesky). The method relies on the
    Nyström approximation which incorporates the 'Q' factor of QR decomposition.
    """
    if (Omega is None):
        assert( l is not None and s is not None and p is not None )
        C = blockHRHT_serial(A.T, l = l, s = s, p = p, seedGlobal = seedGlobal)
        B = blockHRHT_serial(C.T, l = l, s = s, p = p, seedGlobal = seedGlobal)
        B = B.T
    else:
        C = A @ Omega
        B = Omega.T @ C
    Ub, Sb, Vb = svd(B)
    L = Ub@np.diag(np.sqrt(Sb))@Vb
    Z = np.linalg.solve(L, C.T).T
    Q, R = np.linalg.qr(Z)
    if min(R.shape) == rank:
        U_t, Sigma, _ = svd(R)
    else:
        U_t, Sigma, _ = scipy.sparse.linalg.svds(R, k=rank)
    U = Q @ U_t
    return U, np.diag(Sigma**2)



