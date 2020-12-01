import sys
from contextlib import suppress
import numpy as np
import pandas as pd
import copy
import random
import numba
import itertools

import pathlib



from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.SeqIO import parse

from .utils import isint, choose_dict

import warnings


def gen_rand_seq(length, letters=['A','T','C','G']):
    """
    Generate a random DNA sequence of defined length
    
    Parameters
    ----------
    length : int
        Length of generated sequence
    letters : list
        List of letters to choose from when generating the sequence.
    Returns
    ---------
    seq : string
         Random sequence
    """
    # Confirm argument types
    if type(letters) not in [list, np.ndarray, pd.core.series.Series]:
        raise TypeError("`letters` has to be a list-type.")
        
    if not isint(length):
        raise TypeError("`length` has to be an integer.")
    elif length < 0:
        raise ValueError("`length` has to be non-negative.")
    else:
        length = int(length)
            
    # Generate random sequence
    seq = ''.join(random.choice(letters) for i in range(length))
    return seq


def scramble(sequence, attempts, preserve_content=True):
    """
    Scramble letters of a given sequence. Returns most distant sequence.
    
    Parameters
    ----------
    sequence : string
        
    attempts : int
        Number of scrambles which are created. Most dissimilar one is chosen.
    preserve_content : bool
        If True, shuffles the existing sequence. If False, a completely arbitrary sequence is created.
    Returns
    -------
    scrambles[np.argmax(distances)] : string
        Most distant scrambled sequence
    """
    # Check argument types
    if not isint(attempts):
        raise ValueError("`attempts` has to be an integer.")
    else:
        # If type float, change to int
        attempts = int(attempts)
        
    
    scrambles = np.empty(attempts, dtype='<U{}'.format(len(sequence)))
    
    # Create scrambles
    for i in range(attempts):
        if preserve_content:
            snip = "".join(random.sample(sequence, len(sequence)))
        else:
            snip = "".join(random.choices(["A", "C", "T", "G"], len(sequence)))
            
        scrambles[i] = snip
    
    # Compute number of mismatches
    distances = np.zeros(attempts)
    for i in range(attempts):
        distances[i] = np.sum([x != y for (x, y) in zip(sequence, scrambles[i])])

    return scrambles[np.argmax(distances)]
    
    

def create_scrambles(
    sequence, 
    windowsize, 
    overlap,
    attempts,
    number=1,
    preserve_content=True, 
    ignore_imperfect_scrambling=False):
    """
    Scramble letters in a window of given sequence. 
    
    Parameters
    ----------
    sequence : string
        
    windowsize : int
        Size of the window within letters are scrambled
    overlap : int
        Number of letters that each window overlaps with the neighbors
    attempts : int
        Number of scrambles which are created. Most dissimilar one is chosen.
    number : int
        Number of scrambles per window.
    preserve_content : bool, default True
        If True, shuffles the existing sequence. If Flase, a completely arbitrary sequence is created.
    ignore_imperfect_scrambling : bool, default False 
        If False, returns an error when scrambles are not created evenly (last scramble does not end at
        sequence end).
        
    Returns
    -------
    scrambled_sequences : list
        List of scrambled sequences
    """
    
    # Check argument types
    if not isint(windowsize):
        raise TypeError("`windowsize` has to be an integer.")
    else:
        # If type float, change to int
        windowsize = int(windowsize)
        
    if not isint(overlap):
        raise TypeError("`overlap` has to be an integer.")
    else:
        # If type float, change to int
        overlap = int(overlap)
        
    if not isint(attempts):
        raise TypeError("`attempts` has to be an integer.")
    else:
        # If type float, change to int
        attempts = int(attempts)
    
    if overlap >= windowsize:
        raise ValueError("overlap cannot be equal to or bigger than windowsize.")
        
    if not isint(number):
        raise TypeError("`windowsize` has to be an integer.")
    else:
        # If type float, change to int
        number = int(number)
        
    # Compute number of scrambles
    l_sequence = len(sequence)
    n_scrambles = (l_sequence - windowsize) / (windowsize - overlap) + 1
    
    # Test if scrambles can created throughout whole site
    if not n_scrambles.is_integer():
        if not ignore_imperfect_scrambling:
            raise ValueError("Cannot make scrambles with windowsize {} and overlap {}.".format(windowsize, overlap))
        else:
            n_scrambles = np.floor(n_scrambles)
            warnings.resetwarnings()
            warnings.warn("Imperfect scrambles. Last scramble is omitted.")
        
    
    # Create array to store sequences and add wild type
    scrambled_sequences = np.empty(number * int(n_scrambles) + 1, dtype='<U160')
    scrambled_sequences[0] = sequence
    
    # Function to find right indices for scrambles
    indices = lambda i: slice((windowsize-overlap) * i, (windowsize-overlap) * i + windowsize)
    
    # Take subsequence and get scramble
    for i in range(int(n_scrambles)):
        temp_sequence = list(sequence)
        for j in range(number):
            temp_sequence[indices(i)] = scramble(sequence[indices(i)], attempts, preserve_content=True)
            scrambled_sequences[i * number + j + 1] = "".join(temp_sequence)
    
    return scrambled_sequences



def create_scrambles_df(
    sequence, 
    windowsize, 
    overlap, 
    attempts,
    number=1,
    preserve_content=True, 
    ignore_imperfect_scrambling=False):
    """
    Scramble letters in a window of given sequence. Returns pandas df with metadata.
    
    Parameters
    ----------
    sequence : string
        
    windowsize : int
        Size of the window within letters are scrambled
    overlap : int
        Number of letters that each window overlaps with the neighbors
    attempts : int
        Number of scrambles which are created. Most dissimilar one is chosen.
    number : int
        Number of scrambles per window.
    preserve_content : bool, default True
        If True, shuffles the existing sequence. If Flase, a completely arbitrary sequence is created.
    ignore_imperfect_scrambling : bool, default False 
        If False, returns an error when scrambles are not created evenly (last scramble does not end at
        sequence end).
        
    Returns
    -------
    scrambled_sequences : Pandas.DataFrame
        DataFrame of scrambled sequences.
    """
    
    # Check argument types
    if not isint(windowsize):
        raise TypeError("`windowsize` has to be an integer.")
    else:
        # If type float, change to int
        windowsize = int(windowsize)
        
    if not isint(overlap):
        raise TypeError("`overlap` has to be an integer.")
    else:
        # If type float, change to int
        overlap = int(overlap)
        
    if not isint(attempts):
        raise TypeError("`attempts` has to be an integer.")
    else:
        # If type float, change to int
        attempts = int(attempts)
        
    if not isint(number):
        raise TypeError("`attempts` has to be an integer.")
    else:
        # If type float, change to int
        number = int(number)
    
    # Get scrambles
    scrambled_sequences = create_scrambles(
        sequence, 
        windowsize, 
        overlap, 
        attempts,
        number,
        preserve_content,
        ignore_imperfect_scrambling)
    
    # Read wild type sequence
    wild_type = scrambled_sequences[0]
    
    # Get number of scrambles
    n_scrambles = int(np.floor((len(sequence) - windowsize) / (windowsize - overlap))) + 1
    
    # Compute start and end positions of scrambles
    start_pos = np.arange(0, int(n_scrambles), 1) * (windowsize - overlap)
    stop_pos = np.arange(0,int(n_scrambles), 1) * (windowsize - overlap) + windowsize
    
    # Account for multiple scrambles per window
    if number > 1:
        start_pos = np.array([start_pos[i] * np.ones([n_scrambles, number])[i, :] for i in range(n_scrambles)]).flatten()
        stop_pos = np.array([stop_pos[i] * np.ones([n_scrambles, number])[i, :] for i in range(n_scrambles)]).flatten()
    
    # Store sequences in data frame
    scramble_df = pd.DataFrame(
        {'start_pos':start_pos, 
         'stop_pos':stop_pos, 
         'sequence':scrambled_sequences[1:]}
    )
    
    # Compute center of scrambles
    scramble_df['center_pos'] = scramble_df[['start_pos','stop_pos']].mean(axis = 1)
    
    return scramble_df


rest_sites = {
    
}


def find_restriction_sites(enzyme, sequence_list):
    """Searches for restriction sites of a specific enzyme in a list of sequences
    
    Parameters
    ----------
    enzyme : Bio.Restriction or string
        Name of the enzyme.
    sequence_list: array-type
        
    """
    
    
    if type(enzyme) == str:
        try:
            getattr(sys.modules[__name__], enzyme)
        except AttributeError:
            raise ValueError("There is no restriction enzyme {}. Check spelling, naming is case-sensitive.".format(enzyme))
        enzyme = getattr(sys.modules[__name__], enzyme)
        
    # Enzyme classes in Biopython are weird but this works
    elif not (type(enzyme).__bases__[-1] == Restriction.RestrictionType): 
        raise TypeError("enzyme has to be either string of Bio.Restriction.Restriction.RestrictionType")
    
    # Transform sequence inputs
    sequence_list = _check_sequence_list(sequence_list)
   
    return [enzyme.search(sequence) for sequence in sequence_list]
    
    
def scan_enzymes(sequence_list, enzymes=[]):
    """Compute number if restriction sites in a list of sequences for list of enzymes.
    
    Parameters
    ----------
    sequence_list : array-like
    
    enzymes : 
    
    """
    
    # Check inputs
    if type(enzymes) not in [list, np.ndarray, pd.core.series.Series]:
        raise TypeError("enzymes has to be list, numpy array or pandas series.")
    else:
        if any([(not type(enz) == str) and (not type(enz).__bases__[-1] == Restriction.RestrictionType) for enz in enzymes]):
            raise TypeError("entries in `enzymes` have to be of type string or Bio.RestrictionType")
            
    # Check inputs
    sequence_list = _check_sequence_list(sequence_list)
    
    # Return list of enzymes if none given
    ret_list = False
    # Choose all commercially available enzymes if none given
    if len(enzymes) == 0:
        enzymes = CommOnly
        ret_list = True
        
    num_sites = np.zeros(len(enzymes))
    for i, enz in enumerate(enzymes):
        sites = find_restriction_sites(enz, sequence_list)
        num_sites[i] = len([ele for sub in sites for ele in sub])
    
    if ret_list:
        return num_sites, enzymes
    else:
        return num_sites
        
    
    
def mutations_det(
    sequence, 
    num_mutants=None, 
    mut_per_seq=1, 
    site_start=0, 
    site_end=-1, 
    alph_type="DNA"
):
    """Creates single or double mutants.
    
    
    Parameters
    ----------
    - sequence : string
        DNA sequence that is going to be mutated.
    - num_mutants : int, default None
        Number of mutant sequences. If None, all possible mutatants are created.
    - mut_per_seq : int, default 1
        Number of mutations per sequence.
    - site_start : int, default 0
        Beginning of the site that is about to be mutated.
    - site_end : int, default -1
        End of the site that is about to be mutated.
    - alph_type : string, default "DNA"
        Can either be "DNA" for letter sequences, or "Numeric" for integer sequences.
        
    Returns
    -------
    """
    
    if not (isint(num_mutants) or num_mutants==None):
        raise TypeError("`num_mutants` is of type {} but has to be integer valued.".format(type(num_mutants)))
        
    if not isint(mut_per_seq):
        raise TypeError("`mut_per_seq` is of type {} but has to be integer valued.".format(type(mut_per_seq)))
        
    if not isint(site_start):
        raise TypeError("`site_start` is of type {} but has to be integer valued.".format(type(site_start)))
        
    if not isint(site_end):
        raise TypeError("`site_end` is of type {} but has to be integer valued.".format(type(site_end)))
        
    if alph_type not in ["DNA", "Numeric"]:
        raise ValueError("`alph_type` has to be either \"DNA\" or \"Numeric\" ")
        
    
    # Get site to mutate
    mutation_window = sequence[site_start:site_end]

    # Compute number of possible mutations
    poss_mutants = np.prod([3 * len(mutation_window) - i for i in range(mut_per_seq)])

    if mut_per_seq > 3:
        raise ("For more than triple mutants, you should select random .")

    if num_mutants == None:
        num_mutants = poss_mutants
    elif num_mutants > poss_mutants:
        warnings.resetwarnings()
        warnings.warn("Trying to create more unique mutants than possible. Choosing maximal number of {} mutants.".format(poss_mutants))
        num_mutants = poss_mutants
    elif num_mutants/poss_mutants < 0.1:
        warnings.resetwarnings()
        warnings.warn("Sampling less than 10 percent of possible mutations. Consider generating mutations at random to save computational recources.")

    # Create list and choose wt as first
    mutants = np.empty(num_mutants+1, dtype='U{}'.format(len(sequence)))
    mutants[0] = sequence
    mutant_indeces = create_mutant_index(mutation_window, num_mutants, mut_per_seq)
    if alph_type == "DNA":
        letters = np.array(["A", "C", "G", "T"])
    elif alph_type == "Numeric":
        letters = np.arange(4)
    else:
        raise ValueError("Alphabet type has to be either \"DNA\" or \"Numeric\"")
    
    for i, x in enumerate(mutant_indeces):
        mutants[i+1] = sequence[0:site_start] + mutate_from_index(mutation_window, x, letters) + sequence[site_end:]
        
    return mutants



def create_mutant_index(sequence, num_mutants, mut_per_seq):
    """Create index for mutations."""
    mutants = []
    # Generate single mutants
    mutants = [(j, i) for i in range(3) for j in range(len(sequence))]
    if mut_per_seq > 1:
        somelists = mut_per_seq * [mutants]
        elements = np.array(list(itertools.product(*somelists)))
        mask = np.array([~(np.equal(*element)).all() for element in elements])
        mutants = elements[mask]
    else:
        mutants = [np.array([x]) for x in mutants]
    
        
    if num_mutants < len(mutants * mut_per_seq):
        mutants_inds = np.random.choice(np.arange(len(mutants), dtype=int), num_mutants, replace=False)
        temp_mutants = []
        for index in mutants_inds:
            temp_mutants.append(mutants[index])
            
        mutants = temp_mutants
    return mutants


def mutate_from_index(sequence, index, alph):
    seq_list = list(sequence)
    for loci, mutation in index:
        seq_list[loci] = filter_mutation(alph, seq_list[loci])[mutation]
    return "".join(seq_list)


def _check_sequence_list(sequence_list):
    if type(sequence_list) not in [list, np.ndarray, pd.core.series.Series]:
        raise TypeError("sequence_list has to be list, numpy array or pandas series.")
    else:
        if any([type(seq) not in [str, Seq] for seq in sequence_list]):
            raise TypeError("entries in `sequence_list` have to be of type string or Bio.Seq.Seq.")
            
    for i,seq in enumerate(sequence_list):
        if type(seq) == str:
            sequence_list[i] = Seq(seq)
            
    return sequence_list
    

def add_primers(sequence_list, primer_index, autocomplete=False, len_to_complete=200):
    """Add a orthogonal primers to sequences and possibily add random positions up to certain length.
    
    Parameters
    ----------
    - sequence_list : list of sequences
        List of sequences.
    - primer_index : int
        Index of the primer pair that is added to the sequences.
    - autocomplete : boolean
        If True, adds random base pairs to the end of the sequence until length matches `len_to_complete`.
    - len_to_complete : int
        Length of total sequence after primers are added.
        
    Returns
    - new_seq_list : list of sequences
        List of modified sequences with primer pairs.
    
    """
    
    # Check input
    if type(sequence_list) not in [list, np.ndarray, pd.core.series.Series]:
        raise TypeError("sequence_list has to be list, numpy array or pandas series.")
    else:
        if any([type(seq) not in [str, Seq] for seq in sequence_list]):
            raise TypeError("entries in `sequence_list` have to be of type string or Bio.Seq.Seq.")
            
    if not isint(primer_index):
        raise TypeError("`primer_index` is of type {} but has to be integer valued.".format(type(primer_index)))
        
    if not isint(len_to_complete):
        raise TypeError("`len_to_complete` is of type {} but has to be integer valued.".format(type(len_to_complete)))
    elif len_to_complete < 0:
        raise ValueError("`len_to_complete` has to be non-negative.")
        
    if type(autocomplete) != bool:
        raise TypeError("`len_to_complete` is of type {} but has to be boolean.".format(type(autocomplete)))
        
    # Import primers       
    local_path = pathlib.Path(__file__).parent.absolute()
    kosprimefwd = list(parse(str(local_path) + '/forward_finalprimers.fasta','fasta'))
    kosprimerev = list(parse(str(local_path) + '/reverse_finalprimers.fasta','fasta'))
    
    # Extract primers that is added
    forward = str(kosprimefwd[primer_index].seq)
    reverse = str(kosprimerev[primer_index].seq)
    
    # Copy sequence_list to not override input
    new_seq_list = copy.deepcopy(sequence_list)
    
    # Add primers to sequences and autocomplete if chosen to
    for i, seq in enumerate(sequence_list):
        new_seq_list[i] = forward+seq+reverse
        if autocomplete and len(sequence_list[i]) < len_to_complete:
            new_seq_list[i] = new_seq_list[i] + gen_rand_seq(len_to_complete - len(new_seq_list[i]))
        elif len(new_seq_list[i]) > len_to_complete:
            warnings.resetwarnings()
            warnings.warn("Sequence {} is longer than {}bp.".format(i, len_to_complete))
      
    return new_seq_list
    

@numba.njit
def filter_letter(x, letter):
    return x != letter


@numba.njit
def filter_mutation(alph, letter):
    j = 0
    for i in range(alph.size):
        if filter_letter(alph[i], letter):
            j += 1
    result = np.empty(j, dtype=alph.dtype)
    j = 0
    for i in range(alph.size):
        if filter_letter(alph[i], letter):
            result[j] = alph[i]
            j += 1
    return result