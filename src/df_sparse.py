from pyscf import gto, scf, lib
import src.linalg as LA

PTR_LIN_DEP_TOL = 0
PTR_SUM_TRU_TOL = 1

ANGMAP = ['S','P','D','F','G','H','I']
ATMSYM = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O']

def df_coefficients(mol, auxmol, dm):
    # (should also work with XR metric)

    int3c = df.incore.aux_e2(mol, auxmol, 'int3c2e', aosym='s1', comp=1)
    int2c = auxmol.intor('int2c2e', aosym='s1', comp=1)
    rhs = np.einsum('ijp,ji->p', int3c, dm)
    caux = scipy.linalg.solve(int2c, rhs)

    return caux

def gen_metric_solver(int2c):
    try:
        j2c = scipy.linalg.cho_factor(int2c, lower=True)
        j2c_solver = lambda v: scipy.linalg.cho_solve(j2c, v, overwrite_b=True)
    except (numpy.linalg.LinAlgError, scipy.linalg.LinAlgError):
        w, v = scipy.linalg.eigh(int2c)
        mask = w > LINEAR_DEP_THRESHOLD
        #logger.debug(mf_grad, 'int2c2e cond = %.4g, drop %d bfns',
        #             w[-1]/w[0], w.size-numpy.count_nonzero(mask))
        v1 = v[:,mask]
        j2c = lib.dot(v1/w[mask], v1.conj().T)
        def j2c_solver(v):
            if v.ndim == 2:
                return lib.dot(j2c, v)
            else:
                return j2c.dot(v)
    return j2c_solver

def ctr_nao_to_unctr_nbas(mol, pmol):
    """
    Regroup aos by uncontracted shells
    """ 
    pcsh = np.zeros((mol.nao, pmol.nbas), dtype=int)
    j, k = 0, 0
    for i in range(mol.nbas):

        l = mol.bas_angular(i)
        nc = mol.bas_nctr(i)
        if mol.cart:
            nsph = (l+1)*(l+2)//2 * nc
        else:
            nsph = (2 * l + 1) * nc
        pcsh[j:j+nsph, k] = 1
    
        k += nc
        j += nsph

    return pcsh

def restrict_atom_products(pmol):

    idx_atom = []
    count = -1
    for i in range(pmol.nbas):

        ai = pmol.bas_atom(i)
        li = pmol.bas_angular(i)
    
        for j in range(pmol.nbas):

            count += 1
        
            aj = pmol.bas_atom(j)
            lj = pmol.bas_angular(j)

            if (ai == aj) and (li + lj) <= MAXANGM:
                idx_atom.append(count)

    nprod = len(idx_atom)
    idx_atom = np.asarray(idx_atom,dtype=int)
    R = np.zeros((pmol.nbas**2, nprod),dtype=int)
    R[idx_atom, np.arange(nprod)] = 1

    return R

def gto_product_Gram(mol, pmol, dm):
    """
    Create the Gram matrix of GTO products
    """

    # first multiply int2e on naos of mol by dm
    eri1 = lib.einsum('pqrs, pq, rs->pqrs', mol.intor('int2e_sph'), dm, dm)
    # second regroup uncontracted naos by shell
    pcsh = ctr_nao_to_unctr_nbas(mol, pmol)
    eri3 = lib.einsum('pqrs, pi, qj, rk, sl->ijkl', eri1, pcsh, pcsh, pcsh, pcsh)
    G = eri3.reshape(pmol.nbas**2, pmol.nbas**2)

    # Initialize Gram matrix of auxbas
    R = restrict_atom_products(pmol)
    Gatom = R.T @ G @ R

    return (Gatom, R)

def extract_pcd3_products(pmol, pivk):
    """
    Store auxiliary basis set atom-wise
    """

    pcd3_params = {iatm: [] for iatm in range(pmol.natm)}
    ip = 0
    for i in range(pmol.nbas):

        li = pmol.bas_angular(i)
        ei = pmol.bas_exp(i)[0] # uncontr
        ai = pmol.bas_atom(i)

        assert( pmol.bas_nctr(i) == 1 ) 

        for j in range(pmol.nbas):
       
            lj = pmol.bas_angular(j)
            ej = pmol.bas_exp(j)[0] # uncontr
            aj = pmol.bas_atom(j)

            assert( pmol.bas_nctr(j) == 1 ) 

            if (ai != aj): 
                # discard non-atom products
                continue

            if (ip in pivk):
                # select products of PCD pivots
                pcd3_params[ai].append([li+lj, ei+ej])

            ip += 1

    return pcd3_params

def dump_pcd3(pcd3_params, pmol, out_format, angmax):
    """
    Take as input a parameter object and return a string
    out_format 0 is GEM_fit, 1 is PySCF
    """

    prev_elem = []
    count = {i: -1 for i in range(pmol.natm)}
    bases = {}
    for k in pcd3_params:

        elem = pmol._atom[k][0]

        # elem exists already, then append
        if elem in prev_elem:
            elem_bas = bases[elem]
        else: 
            elem_bas = ""
            count[elem] = 0

        prev_elem.append(elem)

        params = pcd3_params[k]
        for (ell,expo) in params:

            # discard functions bigger than MAXANG
            if ell > angmax: continue

            # first line 
            if out_format == 0:
                elem_bas += "   1  %i\n" %ell
            elif out_format == 1:
                elem_bas += "%s   %s\n" %(elem,ANGMAP[int(ell)])
        
            # expo and coeff
            elem_bas += "    %.5f\t1.0\n" %expo

            count[elem] += 1

        if out_format == 0:
            elem_bas = "   %i\n%s" %(count[elem],elem_bas)

        bases[elem] = elem_bas


    return bases

def make_pcd3(mol, angmax, cart, tol, prefactor, out_format=1):
    """
    PCD3 is the Pivoted Cholesky selected basis
    with three centers (for water)

    out_format 0 is GEM_fit, 1 is PySCF

    The target nao size is prefactor*mol.nbas
    """
    
    mol.cart = cart
    mol.build()
    pmol, ctr_coeff = mol.decontract_basis()

    mf = scf.RKS(mol)
    mf.xc = 'b3lyp'
    mf.kernel()
    dm = mf.make_rdm1()

    # Initialize the Gram matrix
    Gatom, R = gto_product_Gram(mol, pmol, dm)

    # Select indices with pivoted Cholesky
    L,piv,_ = LA.PCF(Gatom,eps0)
    Maux = int(prefactor*mol.nbas)
    pivk = piv[:Maux]
    Pk = LA.prmt_mat(Gatom.shape[0], Maux, piv)

    # Prepare extraction as string
    params = extract_pcd3_params(pmol, pivk)
    bases = dump_pcd3(pcd3_params, pmol, out_format, angmax)

    return (bases, Pk, R, Gatom)

def quality(auxbasis, caux):
    """
    Compute general properties of the basis such as
        * condition number
        * vector space span error (orthogonal projection error)
        * modified Gram Schmidt error
        * non-regroup fit error
        * atomic charges
    """

    # Compute the Gram matrix of tha pcd3 basis
    G = auxmol.intor('int1e_ovlp_sph') # show that this yields worst fit res
    # regroup to uncontracted shells
    csh = nao_to_nbas(auxmol)
    Ggrp = lib.einsum('pq, pi, qj-> ij', G, csh, csh) # this regroup be better
    
    Pk = LA.prmt_mat(nprod, Maux, piv)
    print("Condition number", np.linalg.cond(Pk.T @ Gatom @ Pk))

    # check the quality of the span of the space of GTO products
    # by fitting to the shell-wise basis
    b = Pk.T @ R.T @ G @ np.ones(pmol.nbas**2)
    x = scipy.linalg.solve(Pk.T @ Gatom @ Pk, b)
    e_ref = np.sum(G)
    err = np.sqrt(e_ref - x.T @ b)
    print("Span error ", err)

    # Compute coefficients on the orthonormal ao basis
    Lk = L[:Maux,:Maux]
    x_ortho = solve_triangular(Lk, b, lower=True)
    res = np.linalg.norm(Lk @ x_ortho - b,2)
    print('Residual solve_triangular: %.2E' %res)

    # Orthonormalisation coefficients and inversion error
    m = np.size(x_ortho)
    C,_ = dtrtri(Lk, lower=True)
    print('Residual dtrtri: %.2E' %np.linalg.norm(C @ Lk - np.eye(Maux),2))

    # Compute error when auxbasis is not regrouped
    auxmol = df.addons.make_auxmol(mol, auxb)
    # (ij|P)
    int3c = df.incore.aux_e2(mol, auxmol, 'int3c2e', aosym='s1', comp=1)
    # (P|Q)
    int2c = auxmol.intor('int2c2e', aosym='s1', comp=1)
    b = np.einsum('ijP,ji->P', int3c, dm)
    # fit using pivoted Cholesky
    # notice that this fit yields a bigger error than the one we computed
    # shell-wise
    eps = 1e-10
    L,piv,k = LA.PCF(int2c,eps)
    P = LA.prmt_mat(auxmol.nao, k, piv)
    x = scipy.linalg.solve(P.T @ int2c @ P, P.T @ b)
    err = np.sqrt(e_ref - x.T @ P.T @ b)

    return 0

def nao_to_nbas(mol)
    """
    Regroup ao by shell
    """
    csh = np.zeros((mol.nbas, mol.nao), dtype=int)
    j = 0
    for i in range(mol.nbas):

        l = mol.bas_angular(i)
        nc = mol.bas_nctr(i)
        if mol.cart:
            nsph = (l+1)*(l+2)//2 * nc
        else:
            nsph = (2 * l + 1) * nc

        csh[i,j:j+nsph] = 1
        j += nsph

    return mol

def pcd3_coefficients(mol, auxmol, dm, tol, bases, Pk, R, Gatom, is_ortho):

    # we basically want to keep the shell-wise grouping 
    # as it yields lower errors (the basis has been optimized for this)

    # Compute the Gram matrix of tha pcd3 basis
    G = auxmol.intor('int1e_ovlp_sph') # show that this yields worst fit res
    # regroup to uncontracted shells
    csh = nao_to_nbas(auxmol)
    Ggrp = lib.einsum('pq, pi, qj-> ij', G, csh, csh) # this regroup be better

    if is_ortho:
        # Compute coefficients on the orthonormal ao basis
        Lk = L[:Maux,:Maux]
        x_ortho = solve_triangular(Lk, b, lower=True)
        res = np.linalg.norm(Lk @ x_ortho - b,2)
        print('Residual solve_triangular: %.2E' %res)
    else:
        b = Pk.T @ R.T @ G @ np.ones(pmol.nbas**2)
        x = scipy.linalg.solve(Pk.T @ Gatom @ Pk, b)

    # Truncate domain 
    eps = err*2
    D = np.abs(np.outer(x_ortho,x_ortho)) 
    is_node, _, _ = graph.minimal_nodes(D, eps, time_lim=1e3)
    print(np.sum(is_node), Maux**2)
    
    return 1

def write_to_file(filename, bases):

    with open(filename, "w") as file:
        for sym in ATMSYM:

            if sym not in bases: 
                continue
            elif not len(bases[sym]):
                raise IndexError("%s atom auxbasis is empty" %sym)
            
            file.write("#BASIS SET %s\n" %sym)
            file.write(bases[sym])


