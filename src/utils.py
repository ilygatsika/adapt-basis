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

def read_matrix_from_file(filename):
    '''
    Read a matrix stored as a column in a text file
    '''

    with open(filename, 'r') as file:
    
        # Read data line by line as string
        data = file.read().split()

        istart = int(data[0]) + 1
        ncol, nrow = int(data[1]), int(data[2])
        mat = np.array([float(line) for line in data[istart:]], dtype=float)
        mat = mat.reshape(ncol, nrow)

    return mat

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

