import numpy as np
import warnings
from scipy.spatial.transform import Rotation as R
import math

"""
    Molecular geometry processing
"""

def parse_xyz(filename, natom=3):
    """
    filename    system of identical molecules in XYZ format
    return list of molecular fragments in string format
    """

    with open(filename,'r') as file:
        data = file.read()

    lines = data.split('\n')
    comment_sep = '#'
    start = 0
    while comment_sep in lines[start]:
        # while the line is a comment
        start += 1
    lines = lines[start:]
    n = len(lines)//natom
    frag = [""]*n
    for i in range(n):
        i_start, i_end = i*natom, (i+1)*natom
        frag[i] = '\n '.join(lines[j].strip() for j in range(i_start,i_end))
    
    return frag

def rotate_water(water, water_other, verbose=0):
    """
    water         array coordinates in Angstrom ordered as OHH
    water_other   array coordinates
    compute R matrix and displ vector such that 

        water_other = water.dot(R.T) + displ

    """
    
    if not (water.shape[0] == 3):
        raise ValueError("only for 3-atom molecules")
   
    # origin is reference oxygen
    O_ref = water[0]
    O_A = water_other[0]

    A = water - O_ref
    B = water_other - O_A

    with warnings.catch_warnings():

        warnings.simplefilter("ignore")
        rot, rssd = R.align_vectors(B, A)

    Q = rot.as_matrix()
    
    if verbose:

        err = np.linalg.norm(B - rot.apply(A))
        print("rotation error ", err)
        print("rotation ortho", np.allclose(Q.T, np.linalg.inv(Q)))
        
        displ = O_A - rot.apply(O_ref)
        errt = np.linalg.norm(water_other - rot.apply(water) - displ)
        print("rotation+translation error ", errt)

    return Q

def linear_scan(fragment, dist):
    '''
    Input is a dimer
    Return a dimer with second monomer translated O-H bond by some distance
    '''
    return 1

def get_angle(a,b):
    '''
    calculate angle between two vectors
    '''
    na = a/np.linalg.norm(a)
    nb = b/np.linalg.norm(b)
    return np.arccos(np.clip(np.dot(na,nb),-1.0,1.0))

def get_coord(mol_str):
    '''
    input is H 0 0 0; H 1 0 0 (dimer)
    output is coordinates in array
    '''
    atom = mol_str.split("\n")
    natom = len(atom)
    coord = np.empty((natom,3), dtype=float)
    for i in range(natom):
        xyz = atom[i].split()
        coord[i,0] = float(xyz[1])
        coord[i,1] = float(xyz[2])
        coord[i,2] = float(xyz[3])
   
    return coord

def mol_coord(mol_str):
    '''
    input is H 0 0 0; H 1 0 0 (molecule)
    output is coordinates in array
    '''
    atom = mol_str.split(";")
    natom = len(atom)
    coord = np.empty((natom,3), dtype=float)
    for i in range(natom):
        xyz = atom[i].split()
        coord[i,0] = float(xyz[1])
        coord[i,1] = float(xyz[2])
        coord[i,2] = float(xyz[3])
   
    return coord

def rotate_point(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy

def create_water(l1,l2,theta):
    '''
    Return coordinates in local frame 
    with length l1 between OH1
    l2 between OH2, O at zero and H1OH2 angle theta
    return water as string
    '''
    omega = theta*0.5
    origin = [0,0]
    point1 = [l1,0]
    point2 = [l2,0]
    x1,y1 = rotate_point(origin, point1, omega)
    x2,y2 = rotate_point(origin, point2, omega)
    water = 'O 0 0 0; H %.9f %.9f 0; H %.9f %.9f 0' \
            % (x1,y1,x2,-y2)
    
    return water

def local_frame(coord, verbose=0):
    '''
    return coordinates in local frame of water
    local frame is O to zero and midpoint of H1H2 on x axis
    assumes order of coordinates is OHH
    
        coord = coord_ocal.dot(R.T) + displ
    '''
    
    # local water
    l1,l2,theta = get_inner_params(coord)
    water_loc = create_water(l1,l2,theta)
    coord0 = mol_coord(water_loc)
    coord -= coord[0]

    with warnings.catch_warnings():

        warnings.simplefilter("ignore")
        rot, rssd = R.align_vectors(coord, coord0)

    Q = rot.as_matrix()
    
    if verbose:

        err = np.linalg.norm(coord - rot.apply(coord0))
        print("rotation error ", err, rssd)
    
    return (Q, water_loc)

def get_inner_params(coord):

    O = coord[0]
    H1 = coord[1] - O
    H2 = coord[2] - O
    l1 = np.linalg.norm(H1)
    l2 = np.linalg.norm(H2)
    theta = get_angle(H1, H2)

    return (l1,l2,theta)

def stat(arr):

    v1 = np.min(arr)
    v2 = np.max(arr)
    v3 = np.mean(arr)
    v4 = np.std(arr)

    return (v1,v2,v3,v4)

def analyze(filename):
    '''
    analysis of inner parameters of cluster
    return mix, max, mean, std
    '''
    cluster = parse_xyz(filename)
    nfrag = len(cluster)
    vec_l1 = np.empty(nfrag, dtype=float)
    vec_l2 = np.empty(nfrag, dtype=float)
    vec_angle = np.empty(nfrag, dtype=float)

    for i in range(nfrag):

        frag = cluster[i]
        coord = get_coord(frag)
        l1,l2,theta = get_inner_params(coord)

        vec_l1[i] = l1
        vec_l2[i] = l2
        vec_angle[i] = theta

    stat_l1 = stat(vec_l1)
    stat_l2 = stat(vec_l2)
    stat_angle = stat(vec_angle)
    
    return (stat_l1, stat_l2, stat_angle)




