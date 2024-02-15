import matplotlib as mpl
import src.basis as basis
from pyscf import gto, scf
import _pickle as cpickle
from os import listdir
from os.path import isfile, join
import numpy as np
import re

'''
    Various utility functions
    (parsers, option setters, orbital formats)
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

                # Split orbital name respecting spaces
                list_cur = re.split(r'(\s+)', cur_orb)
                nctr = int(list_cur[-1])
                outdata += list_cur[2] + list_cur[3] + list_cur[4] + '\n'
                 
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

def count_orbitals(molecule, path):
    '''
    Count number of orbitals on a given molecule
    for all bases in a folder (GAMESS US format)

    molecule is string, eg 'HHO' for water
    '''

    filenames = [f for f in listdir(path) if isfile(join(path, f))]
    # loop bases in directory
    for filename in filenames:
        
        if not (filename[-3:] == 'bas'):
            continue

        with open(path + filename, 'r') as file:

            # Read file line by line
            lines = [data.rstrip() for data in file]
    
        # loop elements
        count = 0
        for element in molecule:
        
            elem_id = elem_map[element]
            iline = 0
            line = lines[iline]
            # loop lines
            while True:
               
                if (line == '$END'): break
                elif (line == elem_id):
                    
                    nctr = 1
                    line = lines[iline + nctr]
                    while (len(line)):

                        # Split orbital name respecting spaces
                        list_cur = re.split(r'(\s+)', line)
                        nctr += int(list_cur[-1]) + 1
                        count += 1
                        line = lines[iline + nctr]

                    break
                else: 
                    iline += 1

                # go to next line
                line = lines[iline]

        print(filename, "\t", count)

    return 1     

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

