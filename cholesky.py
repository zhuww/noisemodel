import numpy as np
import ctypes
import scipy.linalg as sl
import time

_chofac = np.ctypeslib.load_library('chofactor', '.')
#_chofac.chofactor.argtypes = [ctypes.c_int, ctypes.c_float_p, ctypes.c_float_p]
#_chofac.chofactor.restype = ctypes.c_void

def _chofactor(N, A, LDA):
    _chofac.chofactor(N, ctypes.byref(A), LDA)
    return A
def _chosolve(N, Nrhs, C, ldc, b, ldb):
    _chofac.chosolve(N, Nrhs, ctypes.byref(C), ldc, ctypes.byref(b), ldb)
    return b

    
def Gchosolve(c, b):
    m,n = c.shape
    if b.ndim == 2:
        l,h = b.shape
    elif b.ndim == 1:
        l, = b.shape
        h = 1
    else:
        raise "Ill dimention of b ", b.shape 

    if not n == l:
        raise "Error: shape mismatch", c.shape, b.shape
    mbynmat = ctypes.c_float * (m*n)
    C = mbynmat(*np.asfarray(c.T.flatten(), dtype=np.float32))
    lbyhmat = ctypes.c_float * (l*h)
    B = lbyhmat(*np.asfarray(b.T.flatten(), dtype=np.float32))
    B = _chosolve(m, h, C, n, B, l)
    return np.reshape(np.asfarray(B), (l,h), order='F')

def Gcholesky(a):
    LDA, N = a.shape
    NbyNmat = ctypes.c_float * (LDA*N)
    A = NbyNmat(*np.asfarray(a.T.flatten(), dtype=np.float32))
    _chofactor(N, A, LDA)
    Out = np.asfarray(A, dtype=np.float32)
    Out = np.reshape(Out, (LDA,N), order='F')
    return Out

if __name__ == '__main__':
    N = 50
    a = np.zeros((N, N)); o = np.zeros((N, N))
    np.fill_diagonal(a, np.ones(N))
    for i in range(N):
        for j in range(i+1,N):
            r = np.random.rand() * 0.001
            a[i][j] = r
            a[j][i] = r
    
    tstart = time.time()
    Out1 = Gcholesky(a)
    now = time.time()
    print "{0}x{1} cholesky factorization took GPU {2} s".format(N, N, now - tstart)
    Out2 = sl.cholesky(a)
    #sl.cho_factor(a)
    print "{0}x{1} cholesky factorization took CPU {2} s".format(N, N, time.time() - now)
    print np.allclose(Out1, Out2)

    M = 10
    b = np.random.randn(N,M)
    now = time.time()
    x1 = Gchosolve(Out1, b)
    print "{0}x{1} cholesky solve took GPU {2} s".format(N, M, time.time() - now)
    now = time.time()
    x2 = sl.cho_solve((Out2, False), b)
    print "{0}x{1} cholesky solve took CPU {2} s".format(N, M, time.time() - now)
    print np.allclose(x1, x2)


