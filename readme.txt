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
    [OK] verify that basis from BSE is the same as Diata file
    [*] read density matrix data from file
    [*] initialize atomic orbital basis from GAMESS US file
    [*] deactivate the decontraction of orbitals before Cholesky
    [*] recover set(idx_selection) from orbital products
    [*] copy the GAMESS US file to out folder and remove the lines that do not
      contain selected orbitals
    [*] force atom-wise consistency during the PCD procedure, since it will change the
      pivot order. La proposition de Susi c'est de regarder si le pivot correspond à
      deux fonctions centrées sur différents atomes. Si oui, il faut choisir les
      deux fonctions comme pivot. Prendre le min ou la somme des deux comme pivot à
      ordonner. Dans ce cas deux fonctions vont entre choisi à l'itération actuelle.
      Sinon une fonction est choisi. Pour implémenter cela, il faut ajouter une
      condition qui vérifie si l'indice du AO correspond à deux centres ou pas. Dans
      ce cas il faut stocker pour chaque indice de AO, les indices des AO avec les
      mêmes paramètres centrés aux différents atomes. 

Long term goals

    [*] voir quel pourcentage des orbitales ont a besoin pour créer les produits
    [*] H2O molecule : réduire cc-pv5z et atteindre la taille de cc-pvqz avec la
      même précision que cc-pv5z


