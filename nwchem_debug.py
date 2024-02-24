from src.utils import read_orbitals as test1 
from src.utils import extract_orbitals as test2

'''
Test to debug NWChem sleected orbital output
to file
'''

file_in = 'dat/6-31g.nw'
file_out = 'dat/6-31g-abs.nw'
geom = 'system/h2o.xyz'
idx = [0,1,2,3,4]

basis = test1(geom, file_in, option='nwchem')
print(basis)
#test2(geom, basis, idx, file_out)



