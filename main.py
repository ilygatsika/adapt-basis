import src.basis as basis
import src.linalg as LA
import src.utils as utils
from pyscf import gto
import numpy as np
import argparse
import sys

"""
    Routine that reduces the size of a given atomic orbital basis set, 
    given a molecular geometry and target size. The result is stored in 
    out/<AO>-<method>_<M>.dat in NWChem format.

    Note: 
    All output basis sets are made up of *primitive* Gaussian-type orbitals. 
    We do not generate contracted orbitals.

    Description:
    Our generation method is named Adapted Basis Set (ABS).
    The metric used to evaluate the orbital overlap is L2.
    The selection is based on selecting the most linear independent orbitals
    that approximately span the full space of orbitals.
    We propose three different methods to handle 
    orbital type constraints imposed to the generated basis set:

        ABS-0       no constraint
        ABS-1       all output orbitals are lower than d-type and before the
                    selection procedure the f,g,h ones were discarded
        ABS-2       all output orbitals are lower than d-type and before the
                    selection procedure we replaced each orbital higher than 
                    d-type by a full spd channel

    ========================================================================
    Usage:

        python3 main.py [--coord] h2o [--dm] input_dm [--out_dir]
        out_dir [--AO] aug-cc-pvdz [--M] 20 [--option] 0

    mol       molecular geometry in Angstrom
    AO        atomic orbital basis set
    M         target number of selected atomic orbitals
    option    desired orbital type constraints for method ABS-<option>
"""

# Input parameters
parser = argparse.ArgumentParser()
parser.add_argument("--coord", 
                    help="molecular geometry (in Angstrom) stored in \
                            system/<coord>.xyz", 
                    type=str)
parser.add_argument("--dm", 
                    help="density matrix of converged CBS solution stored in \
                            dat/<dm>.txt (optional)", 
                    type=str, default=None)
parser.add_argument("--out_dir", 
                    help="directory for placing output results", 
                    type=str)
parser.add_argument("--AO", 
                    help="atomic orbital basis set", 
                    type=str)
parser.add_argument("--M_target", 
                    help="target number of selected atomic orbitals",
                    type=int)
parser.add_argument("--option", 
                    help="orbital type constraint for ABS-<option> method",
                    type=int, default=0)
args = parser.parse_args()

# Read user input
coord = args.coord
dm_in = args.dm
out_dir = args.out_dir
ao_basis = args.AO
M_target = args.M_target
option = args.option

# Input filenames
coord = 'system/%s.xyz' %coord
dm_file = 'dat/%s.txt' %dm_in
ao_nw = 'dat/%s.nw' %ao_basis # NwChem
ao_bas = 'dat/%s.bas' %ao_basis # GAMESS US
out_file = out_dir+'/%s-abs_%i_%i.bas' %(ao_basis, M_target, option)

print("Reading coord from %s file." %coord)

# Create molecule in PySCF
input_basis = {
        'H': gto.basis.load(ao_nw, 'H'),
        'O': gto.basis.load(ao_nw, 'O')
        }
mol = gto.M(atom=coord, basis=input_basis, cart=True)
nelec = np.sum(mol.nelec)

# AO basis on fixed geometry
print("Number of atomic orbitals (M must be smaller)\t\t", mol.nbas)
print("Number of orbital components\t\t\t\t", mol.nao)

if (mol.nbas < M_target):
    raise ValueError("target M is too large")

# Get density matrix
if (dm_in is not None):

    # Read density matrix
    # Important: multiply by 2 for occupation
    print("\nImporting CBS density")
    dm = 2 * utils.read_matrix_from_file(dm_file)

    print("\nImporting overlap")
    S = utils.read_matrix_from_file("dat/H2O_ao_overlap.txt")

    #print(mol.ao_labels())
    print("Reference", np.trace(S @ dm))

    # Should be equal to electron number
    nelec_val = np.trace(mol.intor("int1e_ovlp") @ dm) 
    print("Number of electrons (should be %i) = %f" %(nelec, nelec_val))

else: 

    # Get Hartree-Fock density matrix
    print("\nRunning Hartree-Fock to get density")
    mf = mol.RHF()
    mf.run()
    dm = mf.make_rdm1()
    nelec_val = np.trace(mol.intor("int1e_ovlp") @ dm) 
    print("Number of electrons (should be %i) = %f" %(nelec, nelec_val))


print("\nStarting reduction with target %i and %i constraint" %(M_target, option)) 

# Assembly Gram matrix of products
metric = "j"
G = basis.get_4c_gram(mol, dm, metric)
R = basis.restrict_atoms(mol)
A = R.T @ G @ R

# Select products
L,piv,k = LA.PCD(A)
pivk = piv[:M_target]

# Recover orbitals in products
aoi = np.arange(mol.nbas)
prd_idx = np.dstack(np.meshgrid(aoi,aoi)).reshape(-1, 2)
atm_idx = basis.prod_gto_atom_mask(mol)
idx_i, idx_j = [], []
for ip in pivk:

    i,j = prd_idx[atm_idx[ip]]
    idx_i.append(i)
    idx_j.append(j)

print("First index")
print(set(idx_i))
print("Second index")
print(set(idx_j))

print("Final choice")
idx = set(idx_i + idx_j)
print(idx)

# Prepare to write to out_file
basis = utils.read_orbitals(coord, ao_bas)
utils.extract_orbitals(coord, basis, idx, out_file)

print("Results written in %s." %out_file)

exit()

# Initialize intermediary basis respecting the option
mol0 = basis.init_adapt(mol, option=option)

# Gram matrix of intermediary AO basis
G0 = basis.get_gram(mol0)

# Perform PCD with very small tolerance
# to obtain rank-k approximation of G
L, piv, k = LA.PCD(G0)

# Select M orbitals
# Note that orbitals to be selected are primitives
basis_M = basis.low_rank_adapt(mol0, piv, M_target)

mol_M = gto.M(atom=coord, basis=basis_M)
print("Number of selected orbitals\t\t\t", mol_M.nbas)

# Print to file
basis.dump(filename, basis_M)

print("Write adapted basis to %s done." %filename)

