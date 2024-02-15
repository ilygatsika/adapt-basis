from src.utils import count_orbitals as test

#test('HHO', 'out/test_3/')

from pyscf import gto

coord = 'system/h2o.xyz'

print('6-31g')
input_basis = {
        'H': gto.basis.load('dat/6-31g.nw', 'H'),
        'O': gto.basis.load('dat/6-31g.nw', 'O')
        }
mol = gto.M(atom=coord, basis=input_basis)
print(mol.nbas)
mf = mol.RHF()
mf.kernel()

ao_nw = 'dat/test.nw'
print(ao_nw)
input_basis = {
        'H': gto.basis.load(ao_nw, 'H'),
        'O': gto.basis.load(ao_nw, 'O')
        }
mol = gto.M(atom=coord, basis=input_basis)
print(mol.nbas)
mf = mol.RHF()
mf.kernel()


