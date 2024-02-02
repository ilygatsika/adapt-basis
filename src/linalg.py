import scipy
import numpy as np
from pyscf.lib import scipy_helper

'''
    Linear algebra routines
'''

# Linear dependency tolerance set to machine precision
LIN_DEP_TOL = 1E-14

def PCD(A, tol=LIN_DEP_TOL):
    '''
    Wrapper for PySCF pivoted Cholesky wrapper for scipy

    L   :: lower triangular matrix
    piv :: array of selected indices (order is important)
    k   :: approximate rank of A for tolerance

    '''
    return scipy_helper.pivoted_cholesky(A, tol=tol, lower=True)

def permutation(n,k,idx):
    '''
    Return permutation matrix of dimension n by k, k<n
    '''

    P = np.zeros((n,k), dtype=int)
    P[idx[:k], np.arange(k)] = 1
    return P

def PCD_solver(A,b,r=None,tol=LIN_DEP_TOL):
    '''
    Solve Ax=b using pivoted Cholesky decomposition
    for low-rank approximation equal to r
    '''

    L,piv,k = PCD(A, tol)
    n = A.shape[0]
    Pinit = permutation(n, k, piv) # select r-first pivots
    if (r is None): r = k
    P = permutation(n, r, piv) # select r-first pivots
    Ak = P.T @ A @ P
    bk = np.dot(P.T, b)
    xk = scipy.linalg.solve(Ak, bk)
    x = np.dot(P, xk)

    #print(tol, np.square(np.diag(L[:k,:k])) )
    #print("cond PCD", np.linalg.cond(Ak))

    #return x
    Lk = L[:r,:r]
    Ainit = Pinit.T @ A @ Pinit
    B = np.zeros(Ainit.shape)
    B[:r,:r] = Lk @ Lk.T
    err = np.linalg.norm(Ainit - B)/np.linalg.norm(Ainit)
    return x,err

def PPCD_solver(A,b,r=None,tol=LIN_DEP_TOL):
    '''
    Solve Ax=b using preconditioner and 
    pivoted Cholesky decomposition
    for low-rank approximation equal to r
    '''

    L,piv,k = PCD(A, tol)
    n = A.shape[0]
    if (r is None): r = k
    P = permutation(n, r, piv) # select r-first pivots
    Ak = P.T @ A @ P
    bk = np.dot(P.T, b)

    D = np.power(np.diag(Ak), -0.5)
    H = np.diag(D) @ Ak @ np.diag(D)
    c = np.diag(D) @ bk

    xk = scipy.linalg.solve(H, c)
    x = P @ np.diag(D) @ xk
    return x


def orth_solver(A,b,r=None,tol=LIN_DEP_TOL):
    '''
    Solve Ax=b by modified Gram-Schmidt orthonormalization 
    based on pivoted Cholesky decomposition
    Remark:
    1/ yk.T @ Linv @ P.T @ A @ P @ Linv.T @ yk - yk.T @ yk
    2/ yk.T @ yk - x.T @ A @ x
    should be small
    '''

    L,piv,k = PCD(A, tol)
    n = A.shape[0]
    if (r is None): r = k
    P = permutation(n, r, piv) # select r-first pivots
    Lk = L[:r,:r]
    bk = np.dot(P.T, b)
    yk = scipy.linalg.solve_triangular(Lk, bk, lower=True)
    Linv,_ = dtrtri(Lk, lower=True)
    C = np.dot(P, Linv.T) # ortho coeff
    #x = np.dot(C, yk)
    # Lk @ Lk.T @ x = bk
    # Lk @ y        = bk, y = Lk.T @ x
    x = scipy.linalg.solve_triangular(Lk.T, yk, lower=False)
    err = np.linalg.norm(np.dot(Linv,Lk) - np.eye(r),2)
    
    #return (x,yk,err)
    return (yk,Linv,P,x)

def TSVD_solver(A,b,r):
    '''
    Solve Ax=b using truncated singular value decomposition
    for low-rank approximation equal to r
    ''' 
    
    U,s,Vh = scipy.linalg.svd(A, full_matrices=False)
    n = A.shape[0]
    #print(np.allclose(A, (U * s) @ Vh))
    #print(np.linalg.norm((U * s) @ Vh - A))

    #print("cond tsvd", np.linalg.cond((U[:r,:r]*s[0:r])@Vh[:r,:r]))

    sinv = np.zeros(n, dtype=float)
    sinv[0:r] = 1./s[0:r] # select r-first singular values
    c = np.dot(U.T,b)
    w = np.dot(np.diag(sinv),c)
    x = np.dot(Vh.T,w)
    #return x
    sr = np.zeros(n,dtype=float)
    sr[0:r] = s[0:r]
    err = np.linalg.norm((U * sr) @ Vh - A)/np.linalg.norm(A)
    return x,err


