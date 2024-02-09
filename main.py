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

        python3 main.py [--coord] h2o [--AO] aug-cc-pvdz [--M] 20 [--option] 0

    mol       molecular geometry in Angstrom
    AO        atomic orbital basis set
    M         target number of selected atomic orbitals
    option    desired orbital type constraints in ABS-<option>
"""

# Input parameters
parser = argparse.ArgumentParser()
parser.add_argument("--coord", 
                    help="coordinates of molecular geometry (in Angstrom)", 
                    type=str)
parser.add_argument("--dm", 
                    help="density matrix of CBS solution", 
                    type=str)
parser.add_argument("--AO", 
                    help="atomic orbital basis set", 
                    type=str)
parser.add_argument("--M_target", 
                    help="target number of selected atomic orbitals",
                    type=int)
parser.add_argument("--option", 
                    help="desired orbital type constraint",
                    type=int)
args = parser.parse_args()

# Read input
coord = args.coord
dm_file = args.dm
ao_basis = args.AO
M_target = args.M_target
option = args.option

# Input filenames
coord_file = "system/%s.xyz" %coord
dm_file = "dat/%s.txt" %dm_file

# Output filename
filename = "out/%s-abs-%i_%i.dat" %(ao_basis,option,M_target)


print("Reading coord from %s file" %coord_file)

# Read density matrix
# Important: multiply by 2 for occupation
print("\nImporting CBS density")
dm = 2 * utils.read_matrix_from_file(dm_file)

basis = {'H': gto.basis.load('dat/cc-pvdz.nw', 'H')}
mol = gto.M(atom=coord_file, basis=basis)

# Should be equal to electron number
nelec_val = np.trace(mol.intor("int1e_ovlp") @ dm) 
nelec = np.sum(mol.nelec)
print("Number of electrons (should be %i) = %f" %(nelec, nelec_val))

# Get Hartree-Fock density matrix
print("\nRunning Hartree-Fock")
mf = mol.RHF()
mf.run()
dm = mf.make_rdm1()
nelec_val = np.trace(mol.intor("int1e_ovlp") @ dm) 
print("Number of electrons (should be %i) = %f" %(nelec, nelec_val))

exit()


# AO basis on fixed geometry
mol = gto.M(atom=coord_file, basis=ao_basis)
print("Number of contracted atomic orbitals\t\t", mol.nbas)

# Decontract the basis
#mol, _ = mol.decontract_basis() 
print("Number of decontracted atomic orbitals\t\t", mol.nbas)

if (mol.nbas < M_target):
    raise ValueError("target M is too large")

print("** Starting reduction using ABS-%i method with target %i"
      %(option,M_target)) 

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

mol_M = gto.M(atom=coord_file, basis=basis_M)
print("Number of selected orbitals\t\t\t", mol_M.nbas)

# Print to file
basis.dump(filename, basis_M)

print("Write adapted basis to %s done." %filename)

