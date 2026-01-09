# adapt-basis

Python code for the automatic generation of atomic orbital Gaussian basis sets using pivoted Cholesky decomposition. This code is used in our paper cited as: Traore, D., Adjoua, O., Feniou, C., Lygatsika, I.-M., Maday, Y., Posenitskiy, E., Hammernik, K., Peruzzo, A., Toulouse, J., Giner, E., & Piquemal, J.-P. (2024). _Shortcut to chemically accurate quantum computing via density-based basis-set correction._ **Communications Chemistry, 7**, Article 269. https://doi.org/10.1038/s42004-024-01348-3

# Setup

        make setup

# Help

        python3 main.py --help

# Example
    
Generate basis sets for the H2 molecule.
        
        ./script/h2.sh 

Output is atomic orbital basis parameter file in GAMESS US format placed in the out/h2 folder.

# Notes

The target M value should be smaller than the total number of molecular orbitals in the reference basis. 

Molecular geometries for tests are placed in the system folder (coordinates in Angstrom).
    
