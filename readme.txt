Usage
        
        make setup
        python3 main.py [--coord] h2o [--AO] aug-cc-pvdz [--M] 20 [--option] 0
 
    Output is atomic orbital basis parameter file in GAMESS US format placed in 
    the out folder.

Help
    
        python3 main.py --help

Example
    
    Reduce H2 molecule cc-pVDZ basis.
        
        ./script/01-h2_cc-pvdz.sh 

Notes

    Tuning of the M value should be decided by trial and error. Try different
    values and evaluate the error on a desired quantity. There is no general rule,
    as the tuning depends on the target desired quantity (energy, density, etc).

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

TODO simulations

    [OK] H2O (wrong geom) shows that 9 orbitals extracted from cc-pv5z
        give better HF energy (-75.8211) than 9 orbitals of 6-31g (-75.7947)
    [OK] H2O molecule : voir si toutes les bases adaptées à partir de cc-pvXz,
        X=D,T,Q,5, de taille cible égale à 9 sont meilleures en précision que 6-31Gs 
        de taille 9 en énergie DBCBS
    [OK] H2O molecule : réduire cc-pv5z et voir si on peut atteindre la taille de 
        cc-pvqz avec la même précision que cc-pv5z
    [*] produce ABS results for densities (HF or CISD), for atom restrictions, 
        for manual forcing and pivot forcing
    [*] matrice de densité CIPSI
    [*] générer des bases de la taille minimal (STO-3G) jusqu'à 14-(15-16) MOs
    [*] faire des tableaux de comparaison contre STO-3G et 6-31g
    [*] molécules H2O et H2
    [*] compare H2 density matrix true vs HF vs CISD

TODO long term
    
    [*] FH molecule
    [*] test properties: dipoles
    [*] use other metrics such as L2
    [*] optimize contraction coefficients using a linear dependence criterion on
        the density 


