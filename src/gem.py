import src.geometry as geometry
import src.df_helper as df_helper
import src.linalg as LA
import src.basis as basis
import src.utils as utils
from pyscf import gto, lib, df
import numpy as np

'''
    Gaussian electrostatic model (GEM) for
    intermolecular interaction energies computation
    implemented in PySCF

    Support for
        - rigid fragments of frozen geometry (TIP3P)
        - polarizable fragments of flexible geometry (AMOEBA) TODO
'''

create_mol = lambda atom, basis, unit: gto.M(atom=atom, basis=basis, 
        symmetry=False, cart=False, unit=unit)

class RigidFragment:
    '''
    A rigid fragment is an electronic density on a basis (Mole object)
    multiplied by coefficients (array). It moves around space preserving its 
    internal geometrical configuration.
    '''

    def __init__(self, mol, coef):

        self.mol = mol
        self.coef = coef
        self.basis = mol.basis
        self.unit = mol.unit
        self.coord = mol.atom_coords(unit=self.unit)
    
    def make_frame(self, atom_global):
        '''
        atom_global the atoms in global frame as string
        Return a fake RigidFragment object on the new frame
        '''

        unit = self.unit
        bas = self.basis
        mol_global = create_mol(atom_global, bas, unit)
        coord_global = mol_global.atom_coords(unit=unit)
        coef_global = self.rotate_coef(coord_global)
        den_global = RigidFragment(mol_global, coef_global)

        return den_global
        
    def rotate_coef(self, coord_global): 
        '''
        Rotate coefficients on basis to global coordinates
        '''

        Q = geometry.rotate_water(self.coord, coord_global)
        u = gto.ao_rotation_matrix(self.mol, Q)
        
        coef_loc = self.coef
        if (coef_loc.ndim == 1):
            coef_global = np.dot(coef_loc, u)
        elif (coef_loc.ndim == 2):
            coef_global = np.dot(u.T, np.dot(coef_loc, u))

        return coef_global

class PolarizableFragment:
    '''
    A polarizable fragment has a flexible density whose centers and coefficients
    change as a function of atomic relative distances
    '''
    
    def __init__(self):

        raise NotImplementedError()

    def make_frame(self, atom_global):
        '''
        The coefficients are no longer direct to compute.
        Either we run a SCF calculation at the new global frame or we go back to
        the local frame, perform the interpolation step and map to the global
        frame.
        '''

        if tag == "qm":
        
            # compute density matrix running a SCF-KS
            a = 1

        elif tag == "df":

            # compute df_coef by solving the least-squares fit
            a = 1

        elif tag == "rdc":
            
            # go back to the local reference H2O frame by isometry (matrix u)
            # compute the perturbation wrt the reference frame
            # for this perturbation compute the interpolation coefficients
            # map to the global frame (matrix inverse u)
            a = 1

        return 1

def interact(intor_name, den1, den2, mask=None): 
    '''
    A PySCF Moleintor object wrapper for computing 
    pairwise interaction between two Fragment objects
    Return intermolecular interaction energy scalar value
    '''
        
    mol1 = den1.mol
    mol2 = den2.mol
    coef1 = den1.coef
    coef2 = den2.coef
    nbas1 = mol1.nbas
    nbas2 = mol2.nbas

    shls = None
    if intor_name.startswith('int2e') or intor_name.startswith('int4c'):
        # (ij|kl)
        shls = (0, nbas1, 0, nbas1, nbas1, nbas1 + nbas2, nbas1, nbas1 + nbas2)
        drv = get_energy_4c
    elif intor_name.startswith('int3c'):
        # (ij|P)
        shls = (0, nbas1, 0, nbas1, nbas1, nbas1 + nbas2)
        drv = get_energy_3c
    elif intor_name.startswith('int2c') or intor_name.startswith('int1e_ovlp'):
        # (P|Q)
        shls = (0, nbas1, nbas1, nbas1 + nbas2)
        drv = get_energy_2c
    elif intor_name.startswith('int1e_rinv'):
        # (ij|rinv)
        shls = (0, nbas1, 0, nbas1) 
        drv = get_energy_nuc
    elif intor_name.startswith('int1c'):
        # (P|rinv)
        shls = (0, nbas1)
        drv = get_energy_nuc

    res = drv(intor_name, mol1, mol2, coef1, coef2, shls, mask)

    return res

def get_energy_4c(intor_name, mol1, mol2, coef1, coef2, shls, mask=None): 
    
    mol = mol1 + mol2
    out = mol.intor(intor_name, comp=1, aosym = 's1', shls_slice = shls)
    val = lib.einsum('pqrs, pq,rs->', out, coef1, coef2)
    return val

def get_energy_3c(intor_name, mol1, mol2, coef1, coef2, shls, maks=None):
        
    if intor_name.startswith('int3c2e'):
        drv = df.incore.aux_e2
        out = drv(mol1, mol2, intor_name, aosym='s1', comp=1)
    elif intor_name.startswith('int3c1e'):
        mol = mol1 + mol2
        out = mol.intor(intor_name, aosym='s1', comp=1, shls_slice = shls)

    val = lib.einsum('pqr, pq,r->', out, coef1, coef2)
    return val

def get_energy_2c(intor_name, mol1, mol2, coef1, coef2, shls, mask=None):

    mol = mol1 + mol2
    out = mol.intor(intor_name, aosym='s1', comp=1, hermi=1, 
            shls_slice = shls, mask=mask)
    val = lib.einsum('pq, p,q->', out, coef1, coef2)
    return val

def get_energy_nuc(self, intor_name, mol1, mol2, coef1, coef2, shls, mask=None):

    mol = mol1 + mol2
    if intor_name.startswith('int1c'):
        drv = int1c_nuc(mol) # our integrator
    else:
        drv = mol.intor
    nao = mol1.nao
    natm = mol2.natm
    offset = mol1.natm
    out = np.empty((nao, nao), dtype=float)
    z = - np.array([mol2.atom_charge(c) for c in range(natm)], dtype=int)
    for c in range(natm):
        with mol.with_rinv_as_nucleus(offset + c):  
            # set center of 1/r to the coordinate of atom c belonging to mol1 
            out += z[c] * mol.intor(intor_name, shls_slice = shls)
    val = lib.einsum('pq, pq->', out, coef1)
    return val

def GEM(*args):
    '''
    This is a shortcut to build up GaussianElectrostaticModel object.
    Args: Same to :func:`GaussianElectrostaticModel.__init__`
    '''
    return GaussianElectrostaticModel(*args)

class GaussianElectrostaticModel:
    '''
    A GEM object for computation of energies with different methods: 
    exact one using qm densities, or density fitting or robust density fitting.
    Input is coordinates of monomer on local frame and basis names. 
    '''

    def __init__(self, atom, aobas, auxbas, metric, unit='Angstrom'):

        self.atom = atom
        self.unit = unit

        self.rank = basis.format_rank(auxbas)
        self.metric = metric
        self.bas_df = auxbas
        self.bas_qm = aobas

        # QM densities by default
        self.with_qm = True
        self.with_df = False
        self.with_sparse = False
    
    def density_fit(self):
        
        self.with_qm = False
        self.with_df = True
        self.with_sparse = False
        self.build()
        return self

    def robust_density_fit(self):
        
        self.with_qm = True
        self.with_df = True
        self.with_sparse = False
        self.build()
        return self

    def sparse_density_fit(self):

        self.with_qm = False
        self.with_df = False
        self.with_sparse = True
        self.build()
        return self

    def build(self):
    
        atom = self.atom
        unit = self.unit
        rank = self.rank
        aobas = self.bas_qm
        metric = self.metric
        with_qm = self.with_qm
        with_df = self.with_df
        with_sparse = self.with_sparse

        # QM coefficients
        mol = create_mol(atom, aobas, unit)
        dm = basis.get_density_matrix(mol)

        # AUX basis
        if (rank and not isinstance(self.bas_df, dict)):
            bas_df = basis.generate_auxbas(atom, aobas, metric, 
                     rank, unit=unit)
            self.bas_df = bas_df
        
        # DF coefficients
        df_coef = None
        if (with_df or with_sparse):
        
            A,b = basis.least_squares(self.bas_df, mol, dm, metric)
            if (with_df):
                solve = LA.PCD_solver
            elif (with_sparse):
                solve = LA.orth_solver
            df_coef = solve(A,b)
            
        # Mask
        mask = None
        if (with_sparse):

            # generate the mask
            mask = 1

        self.mask = mask
        self.coef_df = df_coef
        self.coef_qm = dm

    def kernel(self, polymer, comp_name):
        '''
        Run energy decomposition computation
        for the energy component named comp_name
        '''
 
        # Parse polymer geometry as string
        fragments = geometry.parse_xyz(polymer)
        
        atom = self.atom
        unit = self.unit
        with_qm = self.with_qm
        with_df = self.with_df
        with_sparse = self.with_sparse
        mask = self.mask

        if (with_df and with_qm):
            mol_qm = create_mol(atom, self.bas_qm, unit)
            mol_df = create_mol(atom, self.bas_df, unit)
            den_qm_local = RigidFragment(mol_qm, self.coef_qm)
            den_df_local = RigidFragment(mol_df, self.coef_df)
            den_local = (den_qm_local, den_df_local)
            drv = GEM_robust_drv
        else: # xor
            if (with_qm):
                bas = self.bas_qm
                coef = self.coef_qm
            elif (with_df):
                bas = self.bas_df
                coef = self.coef_df
            mol = create_mol(atom, bas, unit)
            den_local = RigidFragment(mol, coef)
            drv = GEM_pair_drv
        
        # Associate integral names to energy component
        if comp_name == 'EE':
            if with_qm and with_df:
                intor_name = ('int3c2e', 'int2c2e')
            elif with_df:
                intor_name = 'int2c2e'
            elif with_qm:
                intor_name = 'int2e'
        elif comp_name == 'XR':
            if with_qm and with_df:
                intor_name = ('int3c1e', 'int1e_ovlp')
            elif with_df:
                intor_name = 'int1e_ovlp'
            elif with_qm: 
                intor_name = 'int4c1e'
        elif comp_name == 'EN':
            if (with_df == with_qm or with_sparse):
                return 0
            if with_df:
                intor_name = 'int1c_rinv'
            elif with_qm:
                intor_name = 'int1e_rinv'
        else:
            raise ValueError('unknown energy component')

        out = drv(intor_name, fragments, den_local, mask)
        
        return out

def GEM_pair_drv(intor_name, fragments, den_local, mask=None):
    '''
    Return res array of interaction energies
    '''

    n = len(fragments)
    npair = n * (n-1) // 2
    res = np.empty(npair, dtype=float)

    # Loop over pairs of fragments
    for i in range(n):

        frag1 = fragments[i]
        den1 = den_local.make_frame(frag1)

        for j in range(i):

            frag2 = fragments[j]
            den2 = den_local.make_frame(frag2)
    
            k = j + ((i-1) * i) // 2
            # (ij|kl) or (P|Q)
            res[k] = interact(intor_name, den1, den2, mask) 

    return res

def GEM_robust_drv(intor_name, fragments, den_local, mask=None):
    '''
    Return res array of interaction energies
    including robust density fitting correction terms
    '''

    assert len(den_local) == 2
    assert len(intor_name) == 2

    den_qm_local, den_df_local = den_local
    intor_name_3c, intor_name_2c = intor_name

    n = len(fragments)
    npair = n * (n-1) // 2
    res = np.empty(npair, dtype=float)

    # Loop over pairs of fragments
    for i in range(n):

        frag1 = fragments[i]
        den_qm1 = den_qm_local.make_frame(frag1)
        den_df1 = den_df_local.make_frame(frag1)

        for j in range(i):

            frag2 = fragments[j]
            den_qm2 = den_qm_local.make_frame(frag2)
            den_df2 = den_df_local.make_frame(frag2)
    
            k = j + ((i-1) * i) // 2
            # (ij|P) + (kl|Q) - (P|Q)
            res[k] = interact(intor_name_3c, den_qm1, den_df2, mask) +   \
                    interact(intor_name_3c, den_qm2, den_df1, mask)  -   \
                    interact(intor_name_2c, den_df1, den_df2, mask)

    return res

# In the flexible geometry case a density fitting step should be added inside
# the loop of i and j. Precisely we should add a density fitting step for
# density fitted densities and an interpolation step for reduced basis
# densities. The function make_frame should be replaced by a generic function
# that tunes the coefficients accordingly.

# Some post-processing functions
def root_mean_square_error(vec_ref, vec_app):

    MSE = np.square(np.subtract(vec_ref, vec_app)).mean()
    res = np.sqrt(MSE)

    return res

def parse_key(options):

    monomer = options.monomer
    cluster = options.cluster
    aobas = options.aobas
    auxbas = options.auxbas
    metric = options.metric
    ecomp = options.ecomp
    sttol = options.sttol
    
    key = (monomer, cluster, aobas, auxbas, metric, ecomp, sttol)
    return key

def store_to_file(outfile, key, vref, vapp1, vapp2, vapp3, verb=0):
    '''
    First four entries are total energies
    next three are RMSE
    Last array are pairwise interaction energies
    '''
    res1 = root_mean_square_error(vref, vapp1)
    res2 = root_mean_square_error(vref, vapp2)
    res3 = root_mean_square_error(vref, vapp3)
    t0 = np.sum(vref)
    t1 = np.sum(vapp1)
    t2 = np.sum(vapp2)
    t3 = np.sum(vapp3)

    res = list()
    if (verb == 0):
        # store errors only
        res = np.array([t0, t1, t2, t3, 
            res1, res2, res3], dtype=float)
    elif (verb == 1):
        # store reference values as well
        n = vref.shape[0]
        res = np.empty(4+3+n, dtype=float)
        res[0:4] = [t0, t1, t2, t3]
        res[4:7] = [res1, res2, res3]
        res[7:] = vref

    utils.store_to_file(outfile, key, res)


