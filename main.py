import src.basis as basis
import src.linalg as LA
import src.utils as utils
from pyscf import gto, ci
import numpy as np
import argparse
import sys

"""
    Routine that reduces the size of a given atomic orbital basis set, 
    given a molecular geometry and target size. The result is stored in 
    out/<AO>-<method>_<M>.dat in GAMESS US format.

    Note:
    Our code supports contracted AO bases given as input.

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
    Example usage:

        python3 main.py --coord h2 --dm H2_alpha_rdm --AO ccpvdz --M_target 4 
"""

# Input parameters
parser = argparse.ArgumentParser()
parser.add_argument("--coord", 
                    help="molecular geometry (in Angstrom) stored in \
                            system/<coord>.xyz", 
                    type=str)
parser.add_argument("--dm", 
                    help="density matrix used as reference: \
                            dat/<dm>.txt, HF, CISD", 
                    type=str)
parser.add_argument("--AO", 
                    help="atomic orbital basis is stored in dat/<AO>.nw \
                            (NWChem) and dat/<AO>.bas (GAMESS US)", 
                    type=str)
parser.add_argument("--M_target", 
                    help="target number of selected atomic orbitals, \
                            default is each value from minimal to nao",
                    type=int, default=0)
parser.add_argument("--gram_atom", 
                    help="use of atomic (1) or full (0) Gram matrix, \
                            default is 1",
                    type=int, default=1)
parser.add_argument("--option", 
                    help="orbital type constraint for ABS-<option> method, \
                            default is 0",
                    type=int, default=0)
parser.add_argument("--out_dir", 
                    help="output basis stored in \
                            <out_dir>/<AO>-<M_target>_g<is_atom>_<dm>.bas \
                            (in GAMESS US)", 
                    type=str)
args = parser.parse_args()

# Read user input
coord = args.coord
dm_in = args.dm
out_dir = args.out_dir
ao_basis = args.AO
M_target = args.M_target
option = args.option
gram_atom = args.gram_atom

# Filenames
coord = 'system/%s.xyz' %coord
dm_file = None
if (dm_in not in ['HF','CISD']): 
    dm_file = 'dat/ref/%s.txt' %dm_in
# Input basis files
ao_nw = 'dat/gto/%s.nw' %ao_basis # NWChem
ao_bas = 'dat/gto/%s.bas' %ao_basis # GAMESS US
# Load basis on different formats
res_gamess = utils.read_orbitals(coord, ao_bas, option='gamess')
res_nwchem = utils.read_orbitals(coord, ao_nw, option='nwchem')

print("\nReading coord from %s file." %coord)

# Create molecule in PySCF
input_basis = utils.create_nw_basis(coord, ao_nw)
mol = gto.M(atom=coord, basis=input_basis)
nelec = np.sum(mol.nelec)

# AO basis on fixed geometry
print("Number of atomic orbitals (M must be smaller)\t\t%i MOs %i" \
        %(mol.nbas, mol.nao))

# Get density matrix
if (dm_file is not None):

    # Read density matrix
    # Important: multiply by 2 for occupation
    print("Importing CBS density")
    dm = 2 * utils.read_matrix_from_file(dm_file)

    #print("Importing overlap")
    #S = utils.read_matrix_from_file("dat/H2O_ao_overlap.txt")

    #print(mol.ao_labels())
    #print("Reference", np.trace(S @ dm))

    # Should be equal to electron number
    nelec_val = np.trace(mol.intor("int1e_ovlp") @ dm) 
    print("Number of electrons (should be %i) = %f" %(nelec, nelec_val))

elif (dm_in == 'HF'): 

    # Get Hartree-Fock density matrix
    print("Running Hartree-Fock to get density:")
    myhf = mol.RHF().run() # Hartree-Fock
    dm = myhf.make_rdm1()
    nelec_val = np.trace(mol.intor("int1e_ovlp") @ dm) 
    print("Number of electrons (should be %i) = %f" %(nelec, nelec_val))

elif (dm_in == 'CISD'): 

    # Get Hartree-Fock density matrix
    print("Running Hartree-Fock+CISD to get density:")
    myhf = mol.RHF().run() # Hartree-Fock
    mf = ci.CISD(myhf).run() # Single, double excitation
    dm_mo = mf.make_rdm1() # this is dm on MOs
    c = mf.mo_coeff 
    dm = c @ dm_mo @ c.T # convert dm on MOs to AOs
    nelec_val = np.trace(mol.intor("int1e_ovlp") @ dm) 
    print("Number of electrons (should be %i) = %f" %(nelec, nelec_val))

# Assembly Gram matrix of products
metric = "j"
G = basis.get_4c_gram(mol, dm, metric)
R = basis.restrict_atoms(mol)

if (gram_atom):
    A = R.T @ G @ R; mask = True
else: 
    A = G; mask = False

# Select products
L,piv,k = LA.PCD(A)

if (M_target == 0):
    # Target size is minimal one
    minimal_mol = gto.M(atom=coord, basis='sto-3g')
    M_target = minimal_mol.nbas
    do_once = False
else: 
    do_once = True

if (mol.nbas < M_target):
    raise ValueError("target size is too large")

# Loop over target sizes
adapt_size = 0
while (adapt_size < mol.nbas):

    # Output basis files
    key = (ao_basis, M_target, gram_atom, dm_in.lower())
    out_file_nw = out_dir+'/%s-%i_g%i_%s.nw' %key # NWChem
    out_file    = out_dir+'/%s-%i_g%i_%s.bas' %key # GAMESS US

    pivk = piv[:M_target]
    idx = basis.orbitals_in_products(mol, pivk, mask=mask)

    # Here we try the multiple pivot strategy
    '''
    print('DEBUG')
    A = G; mask = False
    pivk = basis.PCD_2pivot(mol, A, M_target)
    idx = basis.orbitals_in_products(mol, pivk, mask=mask)
    '''
    
    # Write to out_file in GAMESS
    utils.extract_orbitals(coord, res_gamess, idx, out_file, option='gamess')

    # Compute Hartree-Fock energy with PySCF (needs NWChem)
    utils.extract_orbitals(coord, res_nwchem, idx, out_file_nw, option='nwchem')
    input_basis = utils.create_nw_basis(coord, out_file_nw)
    mol_adapt = gto.M(atom=coord, basis=input_basis, unit='A', verbose=0)
    mf_adapt = mol_adapt.RHF()
    adapt_size = mol_adapt.nbas # this is the true size (=/= target)
    adapt_energy = mf_adapt.kernel()
    
    # Print result
    print("Target = %3.i, energy = %.10f, True=\t\t%3.i MOs %3.i\t saved %s" \
            %(M_target, adapt_energy, adapt_size, mol_adapt.nao, out_file))
    
    # Increment target size
    M_target += 1

    if (do_once): break

# Below is the orbital type constraint, not useful for the moment
"""

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
"""

