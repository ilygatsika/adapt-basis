'''
Customizing Hamiltonian for SCF module.

Three steps to define Hamiltonian for SCF:
1. Specify the number of electrons. (Note mole object must be "built" before doing this step)
2. Overwrite three attributes of scf object
    .get_hcore
    .get_ovlp
    ._eri
3. Specify initial guess (to overwrite the default atomic density initial guess)

Note you will see warning message on the screen:

        Overwritten attributes  get_ovlp get_hcore  of <class 'pyscf.scf.hf.RHF'>

'''

from pyscf import gto, scf, ao2mo
from pyscf.tools import fcidump
import numpy as np

def select_components(coord, basis, idx, filename):
    
    mol = gto.M(atom=coord, basis=basis)
    mf = scf.RHF(mol).run()
    print('should be', mol.nelectron, np.trace(mf.get_ovlp() @ mf.make_rdm1()))
    print('original MOs', mol.nao)

    # SCF for the custom Hamiltonian
    mol_custom = gto.M(atom=coord)                                                                                                                                                                                                    
    mol_custom.incore_anyway = True                                                                                                                                                                                                   
    mol_custom.nelectron = mol.nelectron                                                                                                                                                                                              
    mf_custom = scf.RHF(mol_custom)                                                                                                                                                                                                   
                                                                                                                                                                                                                                      
    # selected indices                                                                                                                                                                                                                
    nao = len(idx)                                                                                                                                                                                                                    
                                                                                                                                                                                                                                      
    # Overwrite three attributes of SCF object                                                                                                                                                                                        
    mf_custom.get_hcore = lambda *args: mf.get_hcore()[np.ix_(idx, idx)]                                                                                                                                                              
    mf_custom.get_ovlp = lambda *args: mf.get_ovlp()[np.ix_(idx,idx)]                                                                                                                                                                 
                                                                                                                                                                                                                                      
    # Two equivalent ways to recover eri as four-index                                                                                                                                                                                
    eri = ao2mo.restore(1, mf._eri, mol.nao)                                                                                                                                                                                          
    mf_custom._eri = ao2mo.restore(8, eri[np.ix_(idx,idx,idx,idx)], nao)                                                                                                                                                              
                                                                                                                                                                                                                                      
    mf_custom.init_guess = mf.make_rdm1()[np.ix_(idx,idx)]                                                                                                                                                                            
    mf_custom.kernel()                                                                                                                                                                                                                
                                                                                                                                                                                                                                      
    print('should be', mol.nelectron, np.trace(mf_custom.get_ovlp() @ mf_custom.make_rdm1()))                                                                                                                                         
                                                                                                                                                                                                                                      
    # Call the second order SCF solver in case converging the DIIS-driven HF method                                                                                                                                                   
    # without a proper initial guess is difficult.                                                                                                                                                                                    
    #print(mf_custom.newton().run())                                                                                                                                                                                                  
    assert(mf_custom.converged)                                                                                                                                                                                                       
    #mf_custom.stability()                                                                                                                                                                                                            
                                                                                                                                                                                                                                      
    print('final MOs', mf_custom.make_rdm1().shape[0])                                                                                                                                                                                
                                                                                                                                                                                                                                      
    # Use the given SCF object to transform the 1-electron and 2-electron integrals then dump them to FCIDUMP.                                                                                                                        
    # Convert an SCF object to FCIDUMP                                                                                                                                                                                                
    fcidump.from_scf(mf_custom, filename)                                                                                                                                                                                             
    print('Result print to %s\n' %filename)                                                                                                                                                                                           
                                                                                                                                                                                                                                      

basis = {'O': gto.parse('''
O    S
      1.172000E+04           7.100000E-04
      1.759000E+03           5.470000E-03
      4.008000E+02           2.783700E-02
      1.137000E+02           1.048000E-01
      3.703000E+01           2.830620E-01
      1.327000E+01           4.487190E-01
      5.025000E+00           2.709520E-01
      1.013000E+00           1.545800E-02
      3.023000E-01          -2.585000E-03
O    S
      1.172000E+04          -1.600000E-04
      1.759000E+03          -1.263000E-03
      4.008000E+02          -6.267000E-03
      1.137000E+02          -2.571600E-02
      3.703000E+01          -7.092400E-02
      1.327000E+01          -1.654110E-01
      5.025000E+00          -1.169550E-01
      1.013000E+00           5.573680E-01
      3.023000E-01           5.727590E-01
O    S
      3.023000E-01           1.000000E+00
O    P
      1.770000E+01           4.301800E-02
      3.854000E+00           2.289130E-01
      1.046000E+00           5.087280E-01
      2.753000E-01           4.605310E-01
O    P
      2.753000E-01           1.000000E+00
O    D
      1.185000E+00           1.0000000'''), 
         'H': gto.parse('''
H    S
      1.301000E+01           1.968500E-02           
      1.962000E+00           1.379770E-01           
      4.446000E-01           4.781480E-01           
      1.220000E-01           5.012400E-01           
H    S
      1.220000E-01           1.000000E+00
H    P
      7.270000E-01           1.0000000''')
         }

# H2 molecule, planar on the z-axis
h2 = 'H 0 0 0; H 0 0 0.8'
# H2O molecule, planar on the yz axis
h2o = '''
O   0.000000    0.000000   0.1173000;
H   0.000000   -0.757200 -0.4692000; 
H   0.000000    0.757200  -0.4692000'''

# selected indices
h2_idx = [0,1,4,5,6,9]
h2o_idx = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23]

filename_h2 = 'dat/h2.example'
filename_h2o = 'dat/h2o.example'

select_components(h2, basis, h2_idx, filename_h2)
select_components(h2o, basis, h2o_idx, filename_h2o)

'''
cc-pvdz basis for Hydrogen
--------------------------
0 s
1 s
2 px
3 py
4 pz
5 s
6 s
7 px
8 py
9 pz
--------------------------
cc-pvdz basis for H2O
--------------------------
0 Os
1 Os
2 Os
3 Opx
4 Opy
5 Opz
6 Opx
7 Opy
8 Opz
9 Odxy
10 Odyz
11 Oz2
12 Oxz
13 Ox2-y2
14 Hs
15 Hs
16 Hpx
17 Hpy
18 Hpz
19 Hs
20 Hs
21 Hpx
22 Hpy
23 Hpz
'''
