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

