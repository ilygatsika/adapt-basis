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
    [OK] copy the GAMESS US file to out folder and remove the lines that do not
        contain selected orbitals (for the moment it is manual)
    [*] Create two methods to force atom-wise consistency. 
        The second method does so during the PCD procedure, since it will change the
        pivot order. La proposition de Susi c'est de regarder si le pivot correspond à
        deux fonctions centrées sur différents atomes. Si oui, il faut choisir les
        deux fonctions comme pivot. Prendre le min ou la somme des deux comme pivot à
        ordonner. Dans ce cas deux fonctions vont entre choisi à l'itération actuelle.
        Sinon une fonction est choisi. Pour implémenter cela, il faut ajouter une
        condition qui vérifie si l'indice du AO correspond à deux centres ou pas. Dans
        ce cas il faut stocker pour chaque indice de AO, les indices des AO avec les
        mêmes paramètres centrés aux différents atomes.
    [OK] expand contractions from BSE data in NWChem format
    [OK] automatize different AO basis sets
    [OK] target M is not equal to actual basis size. Store true basis sizes

TODO simulations

    [OK] toy simulation shows that 9 orbitals extracted from cc-pv5z
        give better HF energy (-75.8211) than 9 orbitals of 6-31g (-75.7947)
    [*] is the -75 normal?
    [*] H2O molecule : voir si toutes les bases adaptées à partir de cc-pvXz,
        X=D,T,Q,5, de taille cible égale à 9 sont meilleures en précision que 6-31Gs 
        de taille 9 en énergie DBCBS
    [*] H2O molecule : réduire cc-pv5z et atteindre la taille de cc-pvqz avec la
        même précision que cc-pv5z


