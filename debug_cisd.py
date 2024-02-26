from pyscf import gto, ci
from src import utils
import numpy as np

#gcisd, ucisd, qcisd, fullci

coord = 'system/h2.xyz'
#ao_nw = 'dat/cc-pvdz.nw'
#input_basis = utils.create_nw_basis(coord, ao_nw)
input_basis = '6-31g'

mol = gto.M(atom=coord, basis=input_basis)
myhf = mol.RHF().run() # Hartree-Fock
mf = ci.CISD(myhf).run() # Single, double excitation
dm_mo = mf.make_rdm1()
c = mf.mo_coeff
dm = c @ dm_mo @ c.T
nelec_val = np.trace(mol.intor("int1e_ovlp") @ dm)
nelec = np.sum(mol.nelec)
print("Number of electrons (should be %i) = %f" %(nelec, nelec_val))

# Hartree-Fock
dm_hf = myhf.make_rdm1()
print(dm_hf, dm)
