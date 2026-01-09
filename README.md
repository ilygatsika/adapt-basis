# adapt-basis

Python code for the automatic generation of atomic orbital Gaussian basis sets using pivoted Cholesky decomposition. This code is used in our paper cited as: Traore, D., Adjoua, O., Feniou, C. et al. _Shortcut to chemically accurate quantum computing via density-based basis-set correction_ . Commun Chem 7, 269 (2024). https://doi.org/10.1038/s42004-024-01348-3


# Setup

        make setup

# Help

        python3 main.py --help

# Example
    
    Generate basis sets for the H2 molecule.
        
        ./script/h2.sh 

    Output is atomic orbital basis parameter file in GAMESS US format placed in 
    the out/h2 folder.

# Notes

    The target M value should be smaller than the total number of molecular
    orbitals in the reference basis. 

    Molecular geometries for tests are placed in the system folder (coordinates in 
    Angstrom).
    
