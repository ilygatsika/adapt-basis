Usage
        
        make setup
        python3 main.py [--coord] h2o [--AO] aug-cc-pvdz [--M] 20 [--option] 0

Help
    
        python3 main.py --help

Example
    
    Reduce H2 molecule cc-pVDZ basis.
        
        ./script/01-h2_cc-pvdz.sh 

Notes

    Tuning of the M value should be decided by trial and error. Try different
    values and evaluate the error on a desired quantity. There is no general rule,
    as the tuning depends on the target desired quantity (energy, density, etc).

    Molecular geometries for tests are in system folder (all data is in Angstrom).
    Molecules are at equilibrium geometry.


# TODO 8/2
==========
* GAMESS US
* desactiver les primitives
* H2O cc-pv5z on aimerait atteindre les tailles de la base cc-pvqz
* Angstrom tableau géométrie H2O
* add Gershgorin circle theorem for PCD initialization wrt maximal overlap
* force atom-wise consistency during the PCD procedure, since it will change the
  pivot order. La proposition de Susi c'est de regarder si le pivot correspond à
  deux fonctions centrées sur différents atomes. Si oui, il faut choisir les
  deux fonctions comme pivot. Prendre le min ou la somme des deux comme pivot à
  ordonner. Dans ce cas deux fonctions vont entre choisi à l'itération actuelle.
  Sinon une fonction est choisi. Pour implémenter cela, il faut ajouter une
  condition qui vérifie si l'indice du AO correspond à deux centres ou pas. Dans
  ce cas il faut stocker pour chaque indice de AO, les indices des AO avec les
  mêmes paramètres centrés aux différents atomes. 
* add a loop over contractions in order to take into account the contraction
  coefficients



