import pymp
import numpy as np
from scipy.special import hyp1f1
from pyscf import gto
import ctypes

"""
    Integral routines for GTOs
"""

_cint = np.ctypeslib.load_library('libcint', '/home/io/Documents/rec/code/libcint/build')

# Number of cores for parallel computation
NC = 8

MAPM2SHELL = ['x','y','z']

ANG_OF    = 1
NCTR_OF   = 3

def to_triplet(xyz):
    """
    input    string: xyz, yyz
    output   tuple: 210, 021
    """
    shell = [0,0,0]
    if not len(xyz.strip()): return shell
    for s in xyz:
        shell[MAPM2SHELL.index(s)] += 1
    return shell

def boys(j,T):
    """
    Boys function
    """
    return hyp1f1(j+0.5,j+1.5,-T) * pow(2*j+1,-1)

def R(t,u,v,n,s,SDx,SDy,SDz,RSD):
    """
    Auxiliary function evaluated recursively
    """
    T = s*RSD*RSD
    val = 0.0
    if t == u == v == 0:
        val += np.power(-2*s,n)*boys(n,T)
    elif t == u == 0:
        if v > 1:
            val += (v-1)*R(t,u,v-2,n+1,s,SDx,SDy,SDz,RSD)
        val += SDz*R(t,u,v-1,n+1,s,SDx,SDy,SDz,RSD)
    elif t == 0:
        if u > 1:
            val += (u-1)*R(t,u-2,v,n+1,s,SDx,SDy,SDz,RSD)
        val += SDy*R(t,u-1,v,n+1,s,SDx,SDy,SDz,RSD)
    else:
        if t > 1:
            val += (t-1)*R(t-2,u,v,n+1,s,SDx,SDy,SDz,RSD)
        val += SDx*R(t-1,u,v,n+1,s,SDx,SDy,SDz,RSD)
    return val

def EH(i,t,A):
    """
    Hermite coefficients for one-center Gaussian
    """
    
    if (t < 0) or (t > i):
        return 0.0
    elif i == t == 0:
        # base case
        return 1
    elif i > 0:
        # reduce i
        e1 = EH(i-1,t-1,A)
        e3 = EH(i-1,t+1,A)
        return 0.5*e1 + (t+1)*e3

def int1c1e_nuc(shell,a,A,X):
    """
    1e Coulomb at center X
    """

    AX = A - X
    RAX = np.linalg.norm(AX)
    l,m,n = np.asarray(shell, dtype=int)
    
    grid = np.asarray([(t,u,v) \
            for t in range(l+1) \
            for u in range(m+1) \
            for v in range(n+1) \
                ])
    ngrid = grid.shape[0]

    val = 0.0
    for z in range(ngrid):

        t,u,v = grid[z]

        # Hermite coefficients
        Et = EH(l,t,A[0])
        Eu = EH(m,u,A[1])
        Ev = EH(n,v,A[2])
        
        val += Et * Eu * Ev * \
                R(t,u,v,0,a,AX[0],AX[1],AX[2],RAX)
    
    val *= 2*np.pi/a
    
    return val

def get_int1c_by_shell(mol, i, X):
    """
    i   index of a shell
    """
    
    # only works for primitive basis sets
    if not ( len(mol._libcint_ctr_coeff(i)) == 1 ):
        raise ValueError("basis set should be uncontracted")

    # spherical basis function
    ia = mol.bas_atom(i)
    A = mol.atom_coords()[ia]
    l = mol.bas_angular(i)
    a = mol.bas_exp(i)[0]
    coef = mol._libcint_ctr_coeff(i).ravel()
    u = gto.mole.cart2sph(l)
    ncart = ( (l + 1) * (l + 2) // 2 )

    vcart = np.empty(ncart, dtype=float)
    
    for ic in range(ncart):

        # cartesian basis function
        xyz = mol.cart_labels(ia)[ic][-1]
        shell = to_triplet(xyz)
    
        val = int1c1e_nuc(shell,a,A,X)

        vcart[ic] = val

    vcart *= coef

    # Transform cartesian to spherical
    vsph = u.T @ vcart

    return vsph

def int1e_nuc(mol, X):
    
    off = gto.ao_loc_nr(mol)
    norb = mol.nao
    npair = norb*(norb+1)//2
    res = np.empty(norb, dtype=float)

    for i in range(mol.nbas):

        val = get_int1c_by_shell(mol, i, X)
        i_start, i_end = off[i], off[i+1]
        res[i_start:i_end] = val
    
    return res


def int4c1e(mol):
    """
    Adapted from github.com/sunqm/libcint/testsuite/test_cint4c1e.py
    Credits: Qiming Sun
    """

    natm = mol.natm
    nbas = mol.nbas
    atm = mol._atm
    bas = mol._bas
    env = mol._env

    natm = ctypes.c_int(natm)
    nbas = ctypes.c_int(nbas)
    c_atm = atm.ctypes.data_as(ctypes.c_void_p)
    c_bas = bas.ctypes.data_as(ctypes.c_void_p)
    c_env = env.ctypes.data_as(ctypes.c_void_p)
    
    fnpp1 = _cint.cint2e_ipip1_sph
    fnp1p = _cint.cint2e_ipvip1_sph
    nullptr = ctypes.POINTER(ctypes.c_void_p)()

    def by_pp(shls, shape):
        buf = np.empty(shape+(9,), order='F')
        fnpp1(buf.ctypes.data_as(ctypes.c_void_p), (ctypes.c_int*4)(*shls),
              c_atm, natm, c_bas, nbas, c_env, nullptr)
        ref = buf[:,:,:,:,0] + buf[:,:,:,:,4] + buf[:,:,:,:,8]
        fnp1p(buf.ctypes.data_as(ctypes.c_void_p), (ctypes.c_int*4)(*shls),
              c_atm, natm, c_bas, nbas, c_env, nullptr)
        ref+=(buf[:,:,:,:,0] + buf[:,:,:,:,4] + buf[:,:,:,:,8])*2
        shls = (shls[1], shls[0]) + shls[2:]
        shape = (shape[1], shape[0]) + shape[2:] + (9,)
        buf = np.empty(shape, order='F')
        fnpp1(buf.ctypes.data_as(ctypes.c_void_p), (ctypes.c_int*4)(*shls),
              c_atm, natm, c_bas, nbas, c_env, nullptr)
        ref+= (buf[:,:,:,:,0] + buf[:,:,:,:,4] + buf[:,:,:,:,8]).transpose(1,0,2,3)
        return ref * (-.25/np.pi)

    # Offset of every shell
    ncomp = np.array([(bas[i,ANG_OF] * 2 + 1) * bas[i,NCTR_OF] for i in
        range(mol.nbas)], dtype=int)
    off = gto.ao_loc_nr(mol)

    # precompute indices using 4-fold permutation symmetry
    idx = np.array([ (i,j,k,l)                  \
            for l in range(nbas.value)          \
                for k in range(nbas.value)             \
                    for j in range(nbas.value)  \
                        for i in range(nbas.value)], dtype=int)
    
    nidx = idx.shape[0]
    norb = mol.nao
    npair = norb*(norb+1)//2
    res = pymp.shared.array((norb,norb,norb,norb), dtype=float)
    with pymp.Parallel(NC) as p:
        for ip in p.range(nidx):
            i,j,k,l = idx[ip]
            shls = (i, j, k, l)
            di = ncomp[i]
            dj = ncomp[j]
            dk = ncomp[k]
            dl = ncomp[l]
            ref = by_pp(shls, (di,dj,dk,dl))
            i_start, i_end = off[i], off[i+1]
            j_start, j_end = off[j], off[j+1]
            k_start, k_end = off[k], off[k+1]
            l_start, l_end = off[l], off[l+1]
            res[i_start:i_end,j_start:j_end,k_start:k_end,l_start:l_end] = ref

    return res


