from pyscf import scf, df, lib
import src.sparse as sparse
import src.linalg as linalg
import scipy

'''
    Density fitting routines 
    for applications to KS-DFT and GEM calculations
'''

def cderi(mol, auxmol):
    '''
    Cholesky decomposed ERI tensor
    Adapted from pyscf/examples/df/10-access_df_integrals.py
    '''

    # ints_3c is the 3-center integral tensor (ij|P), where i and j are the
    # indices of AO basis and P is the auxiliary basis
    ints_3c2e = df.incore.aux_e2(mol, auxmol, intor='int3c2e')
    ints_2c2e = auxmol.intor('int2c2e')

    nao = mol.nao
    naux = auxmol.nao

    # Compute the DF coefficients (df_coef) and the DF 2-electron (df_eri)
    df_coef = scipy.linalg.solve(ints_2c2e, ints_3c2e.reshape(nao*nao, naux).T)
    df_coef = df_coef.reshape(naux, nao, nao)

    df_eri = lib.einsum('ijP,Pkl->ijkl', ints_3c2e, df_coef)

    return df_eri

class DensityFitting:
    '''
    Wrapper for different methods of computing density fitting coefficients
    Input is matrix A and vector rhs of least-squares system Ax=rhs
    Different methods are specified as options in kernel routine
    '''

    def __init__(self, A, rhs, lin_dep_tol = 1e-14, sum_trun_tol = 1e-10,
            rdc_bas=None, rdc_coef=None):

        self.A = A
        self.rhs = rhs
        self.ldtol = lin_dep_tol
        self.sttol = sum_truc_tol
        self.rdc_bas = rdc_bas
        self.rdc_coef = rdc_coef

    def kernel(self, orthonrml=False, sparsify=False, interpolate=False):

        # Generate auxiliary basis if necessary
        rank = basis.format_rank(aux_bas)
        if rank:
            aux_bas = basis.generate_auxbas(monomer, ao_bas, metric, rank, unit=unit)

        if interpolate:
            assert (not orthonorml) and (not sparsify)

        ldtol = self.ldtol
        if orthonorml:
            x,y = linalg.orth_solver(A,b,r,ldtol,verbose=1)
            w = y
        else:
            x = linalg.PCD_solver(A,b,ldtol)
            w = x

        if sparsify:
 
            self.sttol = sttol
            W = np.abs(np.outer(w,w))
            mask_orth = sparse.make_anti_triu(W, sttol)
            # post-process the mask in order to apply to A
            mask = np.dot(mask_orth, np.outer())

        elif interpolate:

            # interpolate coefficients on reduced basis
            a = 1

        else:
            
            a = 1
            # perform PCD with lin dep tol
            # solve system LLTx = b

        # Compute Gram matrix of monomer

        # Compute right-hand side
    
        # Solve the least-squares
        df_coef = drv(A, rhs)

        return df_coef



def pcd(A, rhs):

    # etc
    return 1

