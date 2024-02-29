import trexio
import h5py
import numpy as np
from pyscf import gto, scf, ao2mo

import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--system",
                     help="Name of the molecular system",
                     type=str)
parser.add_argument("--geom_file",
                     help="Geometry file, in angström.",
                     type=str)
parser.add_argument("--basis",
                     help="Basis set.",
                     type=str)
parser.add_argument("--charge",
                     help="Charge. Default=0.",
                     type=int,
		     default=0)
parser.add_argument("--spin",
                     help="Total spin = number of unpaired electrons. 2S. Default=0",
                     type=int,
		     default=0)

args = parser.parse_args()

system          = args.system
geom_file       = args.geom_file
basis_input     = args.basis
charge_input    = args.charge
spin_input      = args.spin

# Definition of the system
mol = gto.Mole()
mol.build(
        atom   = geom_file,
        basis  = basis_input,
        charge = charge_input,
        spin   = spin_input,
        )

mf = scf.RHF(mol)
mf.kernel()
dm1 = mf.make_rdm1()

hcore_ao = mol.intor_symmetric('int1e_kin') + mol.intor_symmetric('int1e_nuc')
hcore_mo = np.einsum('pi,pq,qj->ij', mf.mo_coeff, hcore_ao, mf.mo_coeff)

eri_4fold_ao = mol.intor('int2e_sph', aosym=4)
eri_4fold_mo = ao2mo.incore.full(eri_4fold_ao, mf.mo_coeff)
print(eri_4fold_mo.shape)

print(mf.mo_coeff)
print(mf.make_rdm2())

n_mo = mf.mo_coeff.shape[1]
buffsize=n_mo**4
buff_index=np.array([(i,j,k,l) for i in range(n_mo) for j in range(n_mo) for k in range(n_mo) for l in range(n_mo)])
print("Shape buff_index= ", buff_index.shape)
print("buffindex=", buff_index)

# Write infos in trexio file
with trexio.File(f"{system}.trexio", 'w', back_end=1) as f:	
	# Nucleus group
	trexio.write_nucleus_num(f, mol.natm)
	trexio.write_nucleus_coord(f, mol.atom_coords())
	# Basis group
	trexio.write_basis_type(f, "Gaussian")
	# ao group
	trexio.write_ao_num(f, int(mol.nao))
	# MO group
	trexio.write_mo_num(f, mf.mo_coeff.shape[1])
	trexio.write_mo_coefficient(f, mf.mo_coeff)
	# mo_1e_int group
	trexio.write_mo_1e_int_overlap(f, mol.intor("int1e_ovlp"))        
	trexio.write_mo_1e_int_kinetic(f, mol.intor("int1e_kin"))        
	trexio.write_mo_1e_int_potential_n_e(f, mol.intor("int1e_nuc"))        
	trexio.write_mo_1e_int_core_hamiltonian(f, hcore_mo)        
	# mo_2e_int group
	trexio.write_rdm_1e(f, mf.make_rdm1()) 
	trexio.write_rdm_2e(f, 0, buffsize, buff_index.flatten() , mf.make_rdm2().flatten())


