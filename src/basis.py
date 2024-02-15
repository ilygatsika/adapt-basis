from pyscf import gto, lib
import scipy
import numpy as np
from os.path import realpath, join

"""
    Routines for generating adaptive basis sets
"""

# PySCF configuration
from pyscf import __config__
setattr(__config__, 'B3LYP_WITH_VWN5', True)

# Global parameters
ANGSYM = ["S","P","D","F","G","H","I"]
ATMSYM = ["H","He","Li","Be","B","C","N","O"]
MAXL = 2 # maximum orbital type (2 for d-type)

def intermediary_basis(mol, option=0, debug=False):
    """
    Create intermediary basis of s,p,d type from given AO basis

    option 0 :: no constraint
    option 1 :: discard higher than d-type (ABS-1)
    option 2 :: add s,p,d orbitals (ABS-2)
    """
    
    # Initialize empty atomic orbital basis
    newbas = {mol._atom[i][0]: [] for i in range(mol.natm)}

    # loop over given orbitals
    for i in range(mol.nbas):
        
        # Get orbital type
        l = mol.bas_angular(i)
        
        # Get atom
        atm = mol.bas_atom(i)
        elem = mol._atom[atm][0]
        
        # Get orbital exponent
        expo = mol.bas_exp(i)[0]
        num_expo = round(expo,9)
        
        # Check orbital type constraint
        if (option != 0):

            if (l > MAXL):

                # If ABS-2: add new s,p,d orbitals
                if (option == 2):
                    
                    for li in range(MAXL+1):
                        newbas[elem].append((li, num_expo))
           
                # If ABS-1: discard directly
                continue
        
        # Add existing orbital
        newbas[elem].append((l, num_expo))
   
    # Remove duplicated orbitals
    count = 0
    for elem in newbas: 
        newbas[elem] = set(newbas[elem])
        count += len(newbas[elem])

    if (debug):
        # Print number of discarded orbitals
        diff = mol.nbas - count 
        print("Number of discarded orbitals: %i" %diff)

    return newbas

def get_4c_overlap(mol, metric):
    """
    Return four-center matrix of dimension mol.nao by mol.nao
    """
    
    if metric.lower() == 'j':
        A = mol.intor("int2e")
    elif metric.lower() == 's':
        A = mol.intor("int4c1e", comp=1)

    return A

def gto_by_angular(mol):
    """
    Return block-diagonal matrix of dimension mol.nbas by mol.nao,
    that regroups mol.nao GTOs by angular momenta
    """

    # list of number of angular momenta components
    p_angm = []
    for i in range(mol.nbas):

        l = mol.bas_angular(i)
        nc = mol.bas_nctr(i)
        if mol.cart: ncomp = ((l+1) * (l+2) // 2) * nc
        else: ncomp = (2 * l + 1) * nc
        p_angm.append(np.ones(ncomp))

    P = scipy.linalg.block_diag(*p_angm)

    return P

def get_4c_gram(mol, dm, metric):
    '''
    Construct 4-index Gram matrix using
    given metric (j for Coulomb or s for overlap)
    '''

    if metric.lower() not in ['j','s']:
        raise ValueError('invalid metric')

    # Gram matrix of uncontracted basis
    G = get_4c_overlap(mol, metric)
    
    # Weight by density matrix
    G_w = lib.einsum('ij, kl, ijkl->ijkl', dm, dm, G)

    # Block matrix of angular momenta
    P = gto_by_angular(mol)
    G_w = lib.einsum('pi, qj, pqrs, kr, ls->ijkl', P.T, P.T, G_w, P, P)
    norb = mol.nbas**2
    
    gram = G_w.reshape(norb, norb)

    return gram

def get_gram(mol):
    """
    Contract matrix blocks per orbital type
    """
    
    # Overlap (Gram) matrix L2 norm
    # is (mol.nao, mol.nao)
    G = mol.intor("int1e_ovlp")

    # Output is (mol.nbas,mol.nbas)
    P = gto_by_angular(mol)
    B = lib.einsum('qj, qr, kr->jk', P.T, G, P)

    return B

def prod_gto_atom_mask(mol):
    '''
    Return indices of GTO products centered on atoms
    '''
    
    atm = np.array([mol.bas_atom(i) for i in range(mol.nbas)], dtype=int)
    atm2c = np.dstack(np.meshgrid(atm,atm)).reshape(-1, 2)
    idx_atom = np.where( atm2c[:,0] == atm2c[:,1] )[0]

    return idx_atom

def restrict_atoms(mol):
    '''
    Return restriction matrix 
    of primitive GTO products to atomic-centered ones
    '''

    # Restriction matrix 
    atom_mask = prod_gto_atom_mask(mol)
    natp = np.size(atom_mask)
    norb = mol.nbas
    R = np.zeros((norb**2, natp), dtype=int)
    R[atom_mask, np.arange(natp)] = 1

    return R


def extract_basis(mol, indices):
    """
    Extract basis from given indices in internal format
    """
    
    # Initialize empty atomic orbital basis
    newbas = {mol._atom[i][0]: [] for i in range(mol.natm)}

    nidx = len(indices)
    
    # loop over selected orbitals
    for i in range(nidx):
        
        si = indices[i]

        # Get orbital type
        l = mol.bas_angular(si)
        
        # Get atom
        atm = mol.bas_atom(si)
        elem = mol._atom[atm][0]
        
        # Get orbital exponent
        expo = mol.bas_exp(si)[0]
        num_expo = round(expo,9)
        
        # Add existing orbital
        newbas[elem].append((l, num_expo))
    
    # Remove duplicated orbitals
    for elem in newbas: 
        newbas[elem] = set(newbas[elem])

    return newbas

def format(mol):
    """
    Convert AO basis of PySCF Mole object
    in internal form. Assumes primitive
    """
    # Initialize empty atomic orbital basis
    newbas = {mol._atom[i][0]: [] for i in range(mol.natm)}

    # loop over orbitals
    for i in range(mol.nbas):
        
        # Get orbital type
        l = mol.bas_angular(i)
        
        # Get atom
        atm = mol.bas_atom(i)
        elem = mol._atom[atm][0]
        
        # Get orbital exponent
        expo = mol.bas_exp(i)[0]
        num_expo = round(expo,9)
        
        # Add existing orbital
        newbas[elem].append((l, num_expo))
    
    # Remove duplicated orbitals
    for elem in newbas: 
        newbas[elem] = set(newbas[elem])

    return newbas

def init_adapt(mol, option=0, debug=False):
    """
    Return intermediary basis for maximal 
    orbital constraint input option
    as a PySCF Mole object
    """

    # Uncontract atomic orbital basis set
    pmol, ctr_coeff = mol.decontract_basis()

    #print("Primitive atomic orbitals\t%i" %pmol.nbas)

    # Create helper basis
    newbas = intermediary_basis(pmol, option=option, debug=debug)
    help_bas = parse_nwchem(newbas)

    # Return helper basis as PySCF object
    help_mol = gto.M(atom=mol.atom, basis=help_bas)

    return help_mol

def low_rank_adapt(mol, piv, M):
    """
    Generate adaptive basis of rank M based on Cholesky pivots
    from intermediary basis stored in mol
    """
    
    indices = piv[:M]
    sub_bas = extract_basis(mol, indices)
    adapt_bas = parse_nwchem(sub_bas)
    
    # Return adapt basis as PySCF object
    # adapt_mol = gto.M(atom=mol.atom, basis=adapt_bas)

    #return adapt_mol
    return adapt_bas

def parse_nwchem(basis_dic):
    '''
    Convert basis from internal dictionary to NWChem dictionary
    All coefficients are set to 1
    '''

    nw_bases = {}
    for elem in basis_dic:
        elem_bas = ""
        for (ell,expo) in basis_dic[elem]:
            elem_bas += "%s   %s\n" %(elem,ANGSYM[int(ell)])
            elem_bas += "    %.9f\t1.0\n" %expo
        nw_bases[elem] = elem_bas
    
    return nw_bases

def dump(filename, bases):
    '''
    Write string bases to file per atom
    '''
    with open(filename, "w") as file:

        file.write('BASIS "ao basis" PRINT\n')
        for sym in ATMSYM:

            if sym not in bases: 
                continue
            elif not len(bases[sym]):
                raise IndexError("%s atom auxbasis is empty" %sym)
            
            file.write("#BASIS SET %s\n" %sym)
            file.write(bases[sym])

        file.write('END')

def format_basis(basis):
    '''
    Format basis name from local file in basis directory
    '''

    if gto.basis._format_basis_name(basis) in gto.basis.ALIAS:
        return basis
    else:
        dirnow = realpath(join(__file__,'../../basis/'))
        user_basis = join(dirnow, basis+'.dat')
        return user_basis

def PCD_atom_wise(M, Nrbas):
    """
    Greedy pivoted Cholesky (PCF ref Schneider et al.)
    Return subindices of reduced basis

    M       K x Naop x Naop Gram matrices
    Nrbas   size of reduced basis
    """

    K, Naop = M.shape[0], M.shape[1]
    Irbas = np.empty(Nrbas, dtype=int)

    # initialize diagonal, index, lower triangular
    d = np.empty((K, Naop), dtype=float)
    for a in range(K):
        d[a] = np.diag(M[a])
    p = np.arange(Naop)
    L = np.zeros((K, Naop, Naop), dtype=float)

    for m in range(Nrbas):
       
        # choose next reduced basis index
        dm = np.zeros((K, Naop), dtype=float)
        for a in range(K):
            dm[a,m:] = d[a,p[m:]]
        i = np.argmax(np.min(dm, axis=0)) # size Naop

        # Find the same AO at different centers
        aos_i = ao_same_center[i]
        nb_aos = np.size(aos_i)

        # Use the AO at all centers as pivots
        for k in range(nb_aos):

            # project all lines to the set of the AOs
            test = 1


        #print(m, p[i], np.min(d, axis=0)[p[i]])
        Irbas[m] = p[i]
        
        # swap indices
        p[i], p[m] = p[m], p[i]
        assert(i >= m)
    
        L[:,m,p[m]] = np.sqrt(d[:,p[m]])

        for i in range(m+1, Naop):
        
            for a in range(K):
                
                # update lower triangular
                s = 0
                for j in range(m):
                    s += L[a,j,p[m]] * L[a,j,p[i]]
                nrm = L[a,m,p[m]]
                L[a,m,p[i]] = (M[a,p[m],p[i]] - s)/nrm
       
            # update diagonal
            d[:,p[i]] -= L[:,m,p[i]]**2
    
    return Irbas
