#-----------------------------------------------------------------------
QP_ROOT=/home/dtraore/programs/qp2   # <----------- ADRESSE DE TON QP2
source ${QP_ROOT}/quantum_package.rc
#-----------------------------------------------------------------------

system=H2 # Nom de la molecule
basis=cc-pvdz # nom de la base
#basis=cc-pvdz.bas # nom du fichier de la base en format GAMESS US

multiplicity=1 # spin multiplicity 2S+1

geom=${system}.xyz # geometry file

# Create de EZFIO database : (-o is optionnal, it is the name of the ezfio folder)
qp create_ezfio -b ${basis} -m ${multiplicity} -o ${system}_${basis}.ezfio ${geom}

# Run Hartree-Fock SCF
qp run scf | tee scf_${system}_${basis}.out

