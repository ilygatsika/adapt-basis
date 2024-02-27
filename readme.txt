Setup

        make setup

Help

        python3 main.py --help

Example
    
    Generate basis sets for the H2 molecule.
        
        ./script/h2.sh 

    Output is atomic orbital basis parameter file in GAMESS US format placed in 
    the out/h2 folder.

Notes

    The target M value should be smaller than the total number of molecular
    orbitals in the reference basis. 

    Molecular geometries for tests are placed in the system folder (coordinates in 
    Angstrom).


TODO list

    [OK] update H2O geometry from Diata file
    [OK] verify that basis from BSE or PySCF is the same as Diata file
    [OK] read density matrix data from file
    [OK] initialize atomic orbital basis from GAMESS US file
    [OK] deactivate the orbital decontraction before Cholesky
    [OK] perform Cholesky on four-index tensor
    [OK] recover set(idx_selection) from orbital products
    [OK] voir quel pourcentage des orbitales on a besoin pour créer les produits
    [OK] remove the GAMESS US lines that do not contain selected orbitals
    [OK] expand contractions from BSE data in NWChem format
    [OK] automatize different AO basis sets
    [OK] target M is not equal to actual basis size. Store true basis sizes
    [OK] fix -75 Hartree-Fock energy bug
    [OK] implement NWChem output of our adapted basis
    [OK] initialise report on tex
    [OK] load minimal basis sets
    [OK] on-off atom restriction in Gram matrix assembly
    [OK] debug the modified PCD algorithm with conventional pivot selection
    [OK] incorporate atom-wise forcing during the pivot selection
    [OK] scripts that automatize tests on basis generation
    [OK] fix CISD error on electrons
    [OK] add F element to all basis sets in dat folder
    [*] clean main program and create hidden routines for compact calls
    [*] remove abs-x constraints in options and source because not used

TODO simulations

    [OK] H2O (wrong geom) shows that 9 orbitals extracted from cc-pv5z
        give better HF energy (-75.8211) than 9 orbitals of 6-31g (-75.7947)
    [OK] H2O molecule : voir si toutes les bases adaptées à partir de cc-pvXz,
        X=D,T,Q,5, de taille cible égale à 9 sont meilleures en précision que 6-31Gs 
        de taille 9 en énergie DBCBS
    [OK] H2O molecule : réduire cc-pv5z et voir si on peut atteindre la taille de 
        cc-pvqz avec la même précision que cc-pv5z
    [OK] générer des bases de la taille minimal (STO-3G) jusqu'à 14-(15-16) MOs
    [OK] faire des tableaux de comparaison contre STO-3G et 6-31g
    [OK] molécules H2O et H2
    [OK] compare H2 density HF vs CISD
    [*] FH molecule
    [*] matrice de densité alternatives gcisd, ucisd, qcisd, full ci
    [*] matrice de densité CIPSI
    [*] H2 molecule generate size 6 adapted basis, add line in h2.pdf with *-6
    [*] FH et frozen core

TODO long term
    
    [*] test properties: dipoles
    [*] use other metrics such as L2
    [*] try orbitals mixing atoms and types, i/o compatiblity and typed bases


