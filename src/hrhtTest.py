import numpy as np
from scipy.linalg import norm
from math import sqrt
from numpy.linalg import svd, cond
from scipy.linalg import hadamard as sp_hadamard
import matplotlib.pyplot as plt
from math import ceil, log
import pandas as pd
from pathlib import Path

plt.ion()


def generate_Sh(s, l, n):
    '''
    Generate s-hashing matrix of dimensions l x n
    This method might be slow but we use this to 
    generate experiments and deduce properties 
    from the s-hashing matrix
    '''
    S = np.zeros((l, n))
    for i in range(n):
        d = np.array([1 if np.random.random() < 0.5 else -1 for i in range(s)])
        choice = np.random.choice(range(l), s, replace=False)
        S[choice, i] = d
    S = 1/sqrt(s)*S
    return S


# def sym(A):
#     l, r = A.shape
#     return np.block( [ [np.zeros((l,l)), A],
#                        [np.transpose(A) , np.zeros((r,r)) ] ]  )


def testOrthogonality(V, p = 2**5, plot = False):
    '''
    Test if the block randomized Hadamard matrix 
    actually preserves orthogonality
    '''
    m, n = V.shape
    r = int(m/p)
    ds = np.array([1 if np.random.random() < 0.5 else -1 for i in range(m)])
    h = sp_hadamard(r, dtype="d") / np.sqrt(r)
    h = np.eye(r)
    blockH = np.block(  [  [h for i in range(p)]       ]  )* ds.reshape([1, -1])
    T = (1/sqrt(p))*blockH@V
    S = T.T@T
    _, sig, _ = svd(T)
    if plot == True:
        # plt.figure()
        # plt.imshow(np.abs(S))
        # plt.colorbar()
        # plt.title("Imshow of (HV)^T(HV), r = " + str(r) + " m = " + str(m) + " n =" + str(n) )
        # plt.figure()
        # plt.imshow(np.abs(S) - np.eye(n))
        # plt.title("Imshow of (HV)^T(HV) - I, r = " + str(r) + " m = " + str(m) + " n =" + str(n) )
        # plt.colorbar()
        plt.figure()
        plt.plot( np.arange(min(n, r)), sig )
        plt.scatter( np.arange(min(n, r)), sig, s = 0.01, alpha = 0.8 )
        plt.title("Singular values of sqrt(m/r)HV, r = " + str(r) + " m = " + str(m) + " n =" + str(n) )
        plt.suptitle("Without scaling")
        plt.figure()
        sig = sqrt(r/m)*sig
        plt.plot( np.arange(min(n, r)), sig )
        plt.scatter( np.arange(min(n, r)), sig, s = 0.01, alpha = 0.8 )
        plt.title("Singular values of sqrt(m/r)HV, r = " + str(r) + " m = " + str(m) + " n =" + str(n) )
        plt.suptitle("With scaling")
    return norm( (S - np.eye(n)), 'fro' )/n, norm( S - np.eye(n), 2), sig

def testEV(V, p = 2**3, N = 25):
    m, n = V.shape
    r = int(m/p)
    meanSig = np.zeros((n, ))
    for k in range(N):
        ds = np.array([1 if np.random.random() < 0.5 else -1 for i in range(m)])
        h = sp_hadamard(r, dtype="d") / np.sqrt(r)
        #h = np.eye(r)
        blockH = np.block(  [  [h for i in range(p)]       ]  )* ds.reshape([1, -1])
        T = blockH@V
        _, thisS, _ = svd(T)
        meanSig += thisS
    return meanSig/N

def testSingularValues(V, p = 2**3, plot = False, l = None, s = None):
    '''
    Test how close the singular values of block HRHT * orthogonal behave
    '''
    _, sigmaV, _ = svd(V) # The true singular values of V, should all be 1
    # Compute the sketching
    m, n = V.shape
    if l is None:
        l = ceil(max( n, log(4/0.5) )/0.5**2)
    if s is None:
        s = ceil( min( l/3, ceil( 0.1*(log(n/0.5)**4)/0.5**6 ) ) )
    r = int(m/p)
    ds = np.array([1 if np.random.random() < 0.5 else -1 for i in range(m)])
    h = sp_hadamard(r, dtype="d") / np.sqrt(r)
    blockH = np.block(  [  [h for i in range(p)]       ]  )* ds.reshape([1, -1])
    T = (1/sqrt(p))*blockH@V
    S = np.zeros( (l, r) )
    for i in range(r):
        choice = np.random.choice(range(l), s, replace=False)
        S[choice, i] = np.array([1 if np.random.random() < 0.5 else -1 for i in range(s)])
    sketchM = S@T
    _, sigmaSketch, _ = svd(sketchM)
    if plot == True:
        ind = np.arange(n)
        indSketched = np.arange(len(sigmaSketch))
        plt.figure()
        plt.plot( ind, sigmaV,  label = 'True singular values', c = "#111441" )
        plt.scatter( ind, sigmaV, color = "#111441" )
        plt.plot( indSketched, sigmaSketch,  label = 'Sketched', c = "#FF5050" )
        plt.scatter( indSketched, sigmaSketch, color = "#FF5050" )
        plt.colorbar()
        plt.title("Imshow of (HV)^T(HV), r = " + str(r) + " m = " + str(m) + " n =" + str(n) )
    return sigmaV, sigmaSketch, cond(sketchM)

def testOrthogonality_2D(V, p = 2**4, height = 2):
    '''
    Test if the block randomized Hadamard matrix 
    actually preserves orthogonality. We use a 2D block
    randomized Hadamard transform. 
    '''
    m, n = V.shape
    prow = p/height
    r = int(m/prow)
    ds = np.array([1 if np.random.random() < 0.5 else -1 for i in range(m)])
    h = sp_hadamard(r, dtype="d") / np.sqrt(r)
    blockH = np.block(  [  [h for i in range(prow)]       ]  )* ds.reshape([1, -1])
    T = blockH@V
    S = T.T@T
    plt.figure()
    plt.imshow(np.abs(S))
    plt.colorbar()
    plt.title("Imshow of (HV)^T(HV), r = " + str(r) + " m = " + str(m) + " n =" + str(n) )
    plt.figure()
    plt.imshow(np.abs(S) - np.eye(n))
    plt.title("Imshow of (HV)^T(HV) - I, r = " + str(r) + " m = " + str(m) + " n =" + str(n) )
    plt.colorbar()
    plt.figure()
    _, sig, _ = svd(T)
    plt.plot( np.arange(min(n, r)), sig )
    plt.title("Singular values of HV, r = " + str(r) + " m = " + str(m) + " n =" + str(n) )
    
def createTestMatrices(m, n, type = 1, param = 0.2):
    '''
    Creates different type of tall skinny tests matrices.
    Type 1: ORTHOGONAL V factor from SVD of a gaussian matrix
    Type 2: ORTHOGONAL first columns of a hadamard matrix
    Type 3: ORTHOGONAL first rows are Hadamard, rest are zeros
    Type 4: ORTHOGONAL rows of identity matrix scatter along
    Type 5: diagonal with polynomial decay
    Type 6: diagonal with exponential decay
    '''
    if type == 1:
        V = np.random.normal(size = (m, n))
        V, _, _ = svd(V, full_matrices=False)
    elif type == 2:
        V = sp_hadamard(m, dtype="d") / np.sqrt(m)
        V = V[:, 0:n]
    elif type == 3:
        V = np.zeros((m,n))
        V[0:n, :] = sp_hadamard(n, dtype="d") / np.sqrt(n)
    elif type == 4:
        I = np.eye(n)
        V = np.zeros((m,n))
        np.random.seed(166297)
        ind = np.random.choice(range(m), n, replace=False)
        V[ind, :] = I
    elif type == 5:
        d = [1.0 for _ in range(n)] + [(2.0 + o) ** (-param) for o in range(n)]
        V = np.zeros((m, n))
        np.fill_diagonal(V, d)
    elif type == 6:
        d = [1.0 for _ in range(n)] + [(10.0) ** (-(o + 1)) for o in range(n)]
        V = np.zeros((m, n))
        np.fill_diagonal(V, d)
    elif type == 7:
        # MNIST DATA SET
        matrixFile = Path("./data/A.csv")
        dataFile = Path("./data/denseDataMNIST.csv")
        if matrixFile.is_file():
            print("Matrix pre saved, just reading it")
            V = np.genfromtxt("sample_data.csv", delimiter=",", dtype=float)
        elif dataFile.is_file():
            print("Data formated, building matrix")
            labels = np.genfromtxt('./data/labelsMNIST.csv', delimiter=",", dtype=float)
            V = build_A(filename = "./data/denseDataMNIST.csv", labels = labels, save = True)
        else:
            print("Data unprocessed, formatting and building")
            V, _ = build_A(filename = "./data/mnist.scale", save = True)
    return V


def read_data(filename, size=784, save=False):
    '''
    Read MNIST sparse data from filename and transforms this into a dense
    matrix, each line representing an entry of the database (i.e. a "flattened"
    image)
    '''
    dataR = pd.read_csv(filename, sep=',', header=None)
    n = len(dataR)
    data = np.zeros((n, size))
    labels = np.zeros((n, 1))
    # Format accordingly
    for i in range(n):
        l = dataR.iloc[i, 0]
        # We know that the first digit is the label
        labels[i] = int(l[0])
        l = l[2:]
        indices_values = [tuple(map(float, pair.split(':')))
                          for pair in l.split()]
        # Separate indices and values
        indices, values = zip(*indices_values)
        indices = [int(i) for i in indices]
        # Fill in the values at the specified indices
        data[i, indices] = values
    if save:
        data.tofile('./data/denseDataMNIST.csv', sep=',', format='%10.f')
        labels.tofile('./data/labelsMNIST.csv', sep=',', format='%10.f')
    return data, labels

def build_A(filename = None, data = None, labels = None, c=10**4, save=False):
    '''
    Function to build A out of a data base using the RBF exp( -|| x i - x j || /
    c). The function solely based on NumPy broadcasting and is therefore
    superior to methods involving Python loops.
    '''
    assert( filename is not None or data is not None)
    if(data is None and labels is None):
        data, labels = read_data(filename = filename, save = True)
    elif(labels is not None):
        data = np.genfromtxt(filename, delimiter=",", dtype=float)
    data_norm = np.sum(data ** 2, axis=-1)
    A = np.exp(-(1/c) * (data_norm[:, None] +
               data_norm[None, :] - 2 * np.dot(data, data.T)))
    if save:
        A.tofile('./data/A.csv', sep=',', format='%10.f')
    return A, labels


