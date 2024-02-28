#TODO
'''
- remove rows and columns associated to discarded components
- remove eri 4-index entries equaly
- set init density matrix guess
- dump to file
'''

# PySCF example adapted from mcscf/40-customizing_hamiltonian.py

# 1D anti-PBC Hubbard model at half filling
n, u = 12, 2.
mol_hub = gto.M()
mol_hub.nelectron = n // 2
# Setting incore_anyway=True to ensure the customized Hamiltonian (the _eri
# attribute) to be used in the post-HF calculations.  Without this parameter,
# some post-HF method (particularly in the MO integral transformation) may
# ignore the customized Hamiltonian if memory is not enough.
mol_hub.incore_anyway = True
h1 = np.zeros([n] * 2, dtype=np.float64)
for i in range(n-1):
    h1[i, i+1] = h1[i+1, i] = -1.
h1[n-1, 0] = h1[0, n-1] = -1.
eri = np.zeros([n] * 4, dtype=np.float64)
for i in range(n):
    eri[i, i, i, i] = u
rhf_hub = scf.RHF(mol_hub)
rhf_hub.get_hcore = lambda *args: h1
rhf_hub.get_ovlp = lambda *args: np.eye(n)
rhf_hub._eri = ao2mo.restore(8, eri, n) # 8-fold symmetry
rhf_hub.init_guess = '1e'
rhf_hub.kernel()


# Notes : init guess can be the true density matrix with missing lines and
# columns

# What is the incore_anyway option: A second choice is to set incore_anyway in cell which forces the program to generate and hold _eri in the mean-field object.

"""
pyscf.tools.fcidump.from_scf(mf, filename, tol=1e-15, float_format=' %.16g', molpro_orbsym=False)[source]¶
Use the given SCF object to transform the 1-electron and 2-electron integrals then dump them to FCIDUMP.
"""

# Second example based on pyscf/examples/cc/40-ccsd_custom_hamiltonian.py
'''
Six-site 1D U/t=2 Hubbard-like model system with PBC at half filling.
The model is gapped at the mean-field level
'''

import numpy
from pyscf import gto, scf, ao2mo, cc

mol = gto.M(verbose=4)
n = 6
mol.nelectron = n
# Setting incore_anyway=True to ensure the customized Hamiltonian (the _eri
# attribute) to be used in the post-HF calculations.  Without this parameter,
# some post-HF method (particularly in the MO integral transformation) may
# ignore the customized Hamiltonian if memory is not enough.
mol.incore_anyway = True

h1 = numpy.zeros((n,n))
for i in range(n-1):
    h1[i,i+1] = h1[i+1,i] = -1.0
h1[n-1,0] = h1[0,n-1] = -1.0
eri = numpy.zeros((n,n,n,n))
for i in range(n):
    eri[i,i,i,i] = 2.0

mf = scf.RHF(mol)
mf.get_hcore = lambda *args: h1
mf.get_ovlp = lambda *args: numpy.eye(n)
mf._eri = ao2mo.restore(8, eri, n)
mf.kernel()


# In PySCF, the customized Hamiltonian needs to be created once in mf object.
# The Hamiltonian will be used everywhere whenever possible.  Here, the model
# Hamiltonian is passed to CCSD object via the mf object.

mycc = cc.RCCSD(mf)
mycc.kernel()
e,v = mycc.ipccsd(nroots=3)
print(e)

