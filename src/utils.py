import matplotlib as mpl
import src.basis as basis
from pyscf import gto, scf
import argparse
import _pickle as cpickle
import numpy as np

'''
    Various utility functions
    (parsers, option setters, test templates)
'''

# Map first letter to elements
elem_map = {'H': 'HYDROGEN', 
            'O': 'OXYGEN'
            }

def read_matrix_from_file(filename):
    '''
    Read a matrix stored as a column in a text file
    '''

    with open(filename, 'r') as file:
    
        # Read data word by word as a list
        data = file.read().split()

        istart = int(data[0]) + 1
        ncol, nrow = int(data[1]), int(data[2])
        mat = np.array([float(line) for line in data[istart:]], dtype=float)
        mat = mat.reshape(ncol, nrow)

    return mat

def what_elements(geom):
    '''
    Map element first letters to full names
    for every atom in a given molecule in xyz format
    respecting the order of atoms in file
    '''
    
    with open(geom, 'r') as file:

        data = file.read().split()
        elem = []
        for word in data:
            if word.isalpha():
                word_full = elem_map[word]
                elem.append(word_full)
    return elem

def read_orbitals(geom, file_in):
    '''
    Read input atomic orbital basis element-wise
    file is in GAMESS US format
    '''

    # GAMESS US file delimiters
    comment = '!'
    start = '$DATA'
    end = '$END'
    
    # Find order of elements in molecule
    elements = what_elements(geom)
   
    # Initialize empty orbital basis per element
    out_basis = {element: {} for element in elements}

    # First pass, extract basis
    with open(file_in, 'r') as file:

        # Read file line by line
        lines = [data.rstrip() for data in file]
        
        # Init line
        nline = 0
        orb_index = 0
        cur_elem = 'NONE'

        # Loop over lines
        while (lines[nline] != end):

            # Current line
            line = lines[nline]
            
            # Ignore line if it is a comment or delimiter
            if line.startswith((comment, start)):
                nline += 1
                continue

            # or if line is empty
            if (len(line) == 0): 
                nline += 1
                continue
            
            # Break line word by word
            raw_line = line.split()

            # If name of element
            if (len(raw_line) == 1):
                cur_elem = raw_line[0]
                nline += 1
                continue

            # If contraction begins in this line
            if raw_line[0].isalpha():

                orb_type, nctr = raw_line
                nctr = int(nctr)
                # Define orbital identifier
                orb = str(orb_index) + " " + line 
                
                # Loop over contraction
                for k in range(nctr):

                    orb_comp = lines[nline + k+1]

                    if (not orb in out_basis[cur_elem].keys()): 
                        out_basis[cur_elem][orb] = []

                    out_basis[cur_elem][orb].append(orb_comp)

                # Go to next orbital
                nline += nctr
                orb_index += 1

            # Go to next line
            nline += 1

    return out_basis

def extract_orbitals(geom, basis, idx, file_out):
    '''
    Extract selected orbitals of read_orbitals output object
    and store to file_out in GAMESS US format

    idx      orbital indices respecting element order in basis
             i.e. orbitals are enumerated for the molecule
    '''

    # Delimiters
    start = '$DATA\n\n'
    end = '$END\n'

    # Find order of elements in molecule
    elements = what_elements(geom)
    outdata = start

    # init orbital counter
    orb_count = 0
    prev_elements = []
    for element in elements:

        # Recover orbitals of element
        orbitals = basis[element]
        norb = len(orbitals)

        # Check if element exists twice in molecule
        if element not in prev_elements:
            outdata += element + '\n'
        else: 
            outdata = outdata[:-1]
    
        for cur_orb in orbitals.keys():

            # Do not double count orbitals across atoms
            if (element in prev_elements) and \
                    (orb_count - norb) in idx:
                continue

            # If the orbital is selected
            if orb_count in idx:

                # Store orbital params
                orbital = basis[element][cur_orb]

                nctr = int(cur_orb[-1])
                orb_param = cur_orb[2:] + '\n'
                outdata += orb_param
                    
                # Write every component to file
                for k in range(nctr): 
                    outdata += orbital[k] + '\n'
    
            # Go to next orbital
            orb_count += 1

        # Store in elements already counted
        prev_elements.append(element)
        outdata += '\n'

    # End file
    outdata += end

    # Finally write to file
    with open(file_out, "w") as file:
        file.write(outdata)

    return 1

def parse_options():

    parser = argparse.ArgumentParser()

    # Monomer geometry
    parser.add_argument('--monomer', '-g', 
            help="filename of monomer geometry in XYZ format", type=str)
    # Cluster geometry
    parser.add_argument('--cluster', '-c', 
            help="filename of cluster geometry in XYZ format", type=str, default=None)
    # Unit of coordinates of input geometries
    parser.add_argument('--unit', '-u', 
            help="unit of coordinates", type=str, default='Angstrom')
    # AO basis
    parser.add_argument('--aobas', '-a', 
            help="atomic orbital basis set", type=str,
            default="")
    # AUX basis
    parser.add_argument('--auxbas', '-x', 
            help="auxiliary basis set", type=str,
            default="")
    # Object to test
    parser.add_argument('--tcase', '-t', 
            help="test case for timer (cases 0-3)", type=int, default=-1)
    # How many times to run the test
    parser.add_argument('--ntest', '-n', 
            help="number of tests to run", type=int,
            default=0)
    # Metric used in the fit
    parser.add_argument('--metric', '-m', 
            help="fit metric J (Coulomb) or S (overlap)", 
            type=str, default='J')
    # Energy component computed with GEM
    parser.add_argument('--ecomp', '-e', 
            help="energy component EE (electron-electron) or EN \
            (electron-nucle) or XR (exchange-repulsion)", 
            type=str, default='')
    # Sum truncation tolerance
    parser.add_argument('--sttol', '-t1', 
    help="sum truncation tolerance value x yields \
            1E-x", type=int, default=-1)
    # Linear dependence tolerance
    parser.add_argument('--ldtol', '-t2', 
    help="linear dependency tolerance \
            value x yields 1E-x (default 13)", type=int, default=13)
    # Filename for output
    parser.add_argument('--out', '-o', 
    help="filename for storing results", type=str, default='')   
    # Verbose for storing results (save interaction arrays or not)
    parser.add_argument('--out_verb', '-v', 
    help="verbose for storing results", type=int, default=0)

    return parser

def init_visu():
    '''
    Set Matplotlib params for plots
    '''

    mpl.rcParams['legend.fontsize'] = 'small'
    mpl.rcParams['lines.linewidth'] = 0.85
    mpl.rcParams["pdf.use14corefonts"] = True
    mpl.rcParams['lines.markersize'] = 4.5
    mpl.rcParams['lines.markerfacecolor'] = 'none'
    mpl.rcParams["legend.edgecolor"] = 'k'
    mpl.rcParams["legend.labelspacing"] = 0.01
    mpl.rcParams["legend.fancybox"] = False
    mpl.rcParams["axes.labelpad"] = 1.4
    mpl.rcParams['axes.linewidth'] = 0.5
    mpl.rcParams['axes.titlesize'] = 'medium'
    # For latex font
    mpl.rcParams['font.size'] = 11.5
    mpl.rcParams['font.family'] = 'STIXGeneral'
    mpl.rcParams['text.usetex'] = True
    mpl.rcParams["text.latex.preamble"].join([
        r"\usepackage{amsmath}",
        ])
    # Number font
    mpl.rcParams['mathtext.fontset'] = 'stixsans'

def parse_test_case(options):
    '''
    Return function and arguments for time benchmark
    This is a wrapper for available tests
    '''
    
    objtag = options.tcase
    system = options.monomer
    unit = options.unit
    aobas = options.aobas
    auxbas = options.auxbas
    metric = options.metric
    rank = basis.format_rank(auxbas)

    if objtag < 0 or objtag > 3: 
        raise ValueError("invalid test")
 
    obj_fun = None
    args = None
    xc = 'b3lyp'

    if (objtag == 0): 
   
        if rank:
            obj_fun = basis.generate_auxbas
            args = (system, aobas, metric, rank, unit)
        else: 
            raise ValueError("test unavailable")

    elif objtag == 1:
        raise NotImplementedError("pending")

    elif objtag == 2:
        
        mol = gto.M(atom=system, basis=aobas, cart=False,
                symmetry=False, unit=unit, verbose=0)
        if len(aobas):
            # the test is for QM basis
            mf = scf.RKS(mol)
        else:
            # The test is for auxiliary basis
            if rank: 
                auxbasis = basis.generate_auxbas(system, aobas, metric, rank)
            else:
                auxbasis = basis.format_basis(auxbas)
            mf = scf.RKS(mol).density_fit(auxbasis=auxbasis)
        
        mf.xc = xc
        obj_fun = lambda mf: mf.kernel()
        args = (mf,)

    elif objtag == 3:
        raise NotImplementedError("pending")

    return (obj_fun, args)

def store_to_file(outfile, key, res):
    '''
    Append existing dictionary in outfile 
    with key and res
    ''' 
    try:
        # Add to entry if file exists 
        with open(outfile, 'rb') as file:
            dic = cpickle.load(file)
            
            if key in dic:
                # append to entry
                cur_res = dic[key]
                cur_res.update(res)
                dic[key] = cur_res
            else:
                # create new entry
                dic[key] = res
    except:
        # Create new dictionary
        dic = {key : res}

    with open(outfile, 'wb') as file:
        cpickle.dump(dic, file)

