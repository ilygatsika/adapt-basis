from src.utils import read_orbitals as test1
from src.utils import extract_orbitals as test2

geom = 'system/h2o.xyz'
file_in = 'dat/cc-pvdz.bas'
file_out = 'out/test_3/cc-pvdz.bas'
idx = [1,2,5,6,10,11]

basis = test1(geom, file_in)
#print(basis)

test2(geom, basis, idx, file_out)

