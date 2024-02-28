import src.basis as basis
import src.linalg as LA
import src.utils as utils
from pyscf import gto, lib
import numpy as np

# input
coord = 'system/h2.xyz'
ao_nw = 'dat/gto/cc-pvtz.nw'
'''
{0, 1, 2, 3, 6, 7, 8, 9}
{0, 1, 2, 5, 14, 15, 16, 19}
Pendant Hartree-Fock: use
mask = 1 1 1 0 0 1 1 1 1 0 0 1
After Hartree-Fock: discard orbitals 
=====
0  s
1  s
2  s
3  p
4  p
5  d
6  s
7  s
8  s
9  p
10 p
11 d
=====
0  s
1  s
2  s
3  px
4  py
5  pz
6  px 
7  py
8  pz
9  xy
10 yz
11 z2
12 xz
13 x2-y2
14 s
15 s
16 s
17 px
18 py
19 pz
20 px 
21 py
22 pz
23 xy
24 yz
25 z2
26 xz
27 x2-y2
'''
M_target = 9

input_basis = utils.create_nw_basis(coord, ao_nw)
mol = gto.M(atom=coord, basis=input_basis)
nelec = np.sum(mol.nelec)

print("Running Hartree-Fock to get density:")
myhf = mol.RHF().run() # Hartree-Fock
dm = myhf.make_rdm1()
nelec_val = np.trace(mol.intor("int1e_ovlp") @ dm) 
print("Number of electrons (should be %i) = %f" %(nelec, nelec_val))

metric = "j"
G = basis.get_4c_gram(mol, dm, metric)
R = basis.restrict_atoms(mol)
#A = R.T @ G @ R; mask = True
A = G; mask = False

# Select products
L,piv,k = LA.PCD(A)
pivk = piv[:M_target]
idx = basis.orbitals_in_products(mol, pivk, mask=mask)
print(idx)

# Now select among separate components
G = basis.get_4c_overlap(mol, metric)
G_w = lib.einsum('ij, kl, ijkl->ijkl', dm, dm, G)
norb = mol.nao**2
A = G_w.reshape(norb, norb)

# Select products
L,piv,k = LA.PCD(A)
pivk = piv[:M_target]
# orbitals in products but for nao
aoi = np.arange(mol.nao)
prd_idx = np.dstack(np.meshgrid(aoi,aoi)).reshape(-1, 2)
idx_i, idx_j = [], []
for ip in pivk:
    i,j = prd_idx[ip]
    idx_i.append(i)
    idx_j.append(j)
idx = set(idx_i + idx_j)
print(idx)
print(mol.nao)
print(mol.nbas)

