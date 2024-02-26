from pyscf import gto
from src import utils

bases = ['sto-3g', '6-31g', '6-31g_st', '6-311g', 'cc-pvdz', 'cc-pvtz', 'cc-pvqz']
         #'cc-pv5z', 'cc-pv6z'] 
molecules = ['h2', 'h2o']

for molecule in molecules:
    print('%s =========== nb of MOs ========== Hartree-Fock' %molecule)
    coord = 'system/%s.xyz' %molecule
    for basis in bases:
        file_nw = 'dat/%s.nw' %basis
        input_basis = utils.create_nw_basis(coord, file_nw)
        mol = gto.M(atom=coord, basis=input_basis, verbose=0)
        myhf = mol.RHF() # Hartree-Fock
        E = myhf.kernel()
        print('%s\t\t%i\t\t\t%.6f' %(basis, mol.nao, E))

