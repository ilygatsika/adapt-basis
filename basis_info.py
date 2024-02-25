from pyscf import gto
from src import utils

bases = ['sto-3g', '6-31g', 'cc-pvdz', 'cc-pvtz', 'cc-pvqz', 'cc-pv5z',
         'cc-pv6z'] 
molecules = ['h2', 'h2o']

for molecule in molecules:
    print('%s =========== nb of MOs' %molecule)
    coord = 'system/%s.xyz' %molecule
    for basis in bases:
        file_nw = 'dat/%s.nw' %basis
        input_basis = utils.create_nw_basis(coord, file_nw)
        mol = gto.M(atom=coord, basis=input_basis)
        print('%s\t\t%i' %(basis, mol.nao))

