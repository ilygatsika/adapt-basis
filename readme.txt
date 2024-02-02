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

