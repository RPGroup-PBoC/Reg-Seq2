import numpy as np
import pandas as pd
import copy
import random
from .gen_utils import isint
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
        raise ValueError("`letters` has to be a list.")
        
    if not isint(length):
        raise ValueError("`length` has to be an integer.")
    else:
        length = int(length)
            
    # Generate random sequence
    seq = ''.join(random.choice(letters) for i in range(length))
    return(seq)


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
        raise ValueError("`windowsize` has to be an integer.")
    else:
        # If type float, change to int
        windowsize = int(windowsize)
        
    if not isint(overlap):
        raise ValueError("`overlap` has to be an integer.")
    else:
        # If type float, change to int
        overlap = int(overlap)
        
    if not isint(attempts):
        raise ValueError("`attempts` has to be an integer.")
    else:
        # If type float, change to int
        attempts = int(attempts)
    
    if overlap >= windowsize:
        raise ValueError("overlap cannot be equal to or bigger than windowsize.")
        
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
    scrambled_sequences = np.empty(int(n_scrambles) + 1, dtype='<U160')
    scrambled_sequences[0] = sequence
    
    # Function to find right indices for scrambles
    indices = lambda i: slice((windowsize-overlap) * i, (windowsize-overlap) * i + windowsize)
    
    # Take subsequence and get scramble
    for i in range(int(n_scrambles)):
        temp_sequence = list(sequence)
        temp_sequence[indices(i)] = scramble(sequence[indices(i)], attempts, preserve_content=True)
        scrambled_sequences[i+1] = "".join(temp_sequence)
    
    return scrambled_sequences



def create_scrambles_df(sequence, windowsize, overlap, attempts, preserve_content=True):
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
    preserve_content : bool
        If True, shuffles the existing sequence. If Flase, a completely arbitrary sequence is created.
        
    Returns
    -------
    scrambled_sequences : Pandas.DataFrame
        DataFrame of scrambled sequences.
    """
    
    # Check argument types
    if not isint(windowsize):
        raise ValueError("`windowsize` has to be an integer.")
    else:
        # If type float, change to int
        windowsize = int(windowsize)
        
    if not isint(overlap):
        raise ValueError("`overlap` has to be an integer.")
    else:
        # If type float, change to int
        overlap = int(overlap)
        
    if not isint(attempts):
        raise ValueError("`attempts` has to be an integer.")
    else:
        # If type float, change to int
        attempts = int(attempts)
    
    # Get scrambles
    scrambled_sequences = create_scrambles(sequence, windowsize, overlap, attempts, preserve_content=True)
    
    # Read wild type sequence
    wild_type = scrambled_sequences[0]
    
    # Get number of scrambles
    n_scrambles = len(scrambled_sequences) - 1
    
    # Compute start and end positions of scrambles
    start_pos = np.arange(0, int(n_scrambles), 1) * (windowsize - overlap)
    stop_pos = np.arange(0,int(n_scrambles), 1) * (windowsize - overlap) + windowsize
    
    # Store sequences in data frame
    scramble_df = pd.DataFrame({'start_pos':start_pos, 'stop_pos':stop_pos, 'sequence':scrambled_sequences[1:]})
    
    # Compute center of scrambles
    scramble_df['center_pos'] = scramble_df[['start_pos','stop_pos']].mean(axis = 1)
    
    return scramble_df




"""
def site_single_mutations(sequence, site_start=0, site_end=-1, alph_type="DNA"):
    
    if alph_type == "DNA":
        letters = np.array(["A", "C", "G", "T"])
        scrambled_sequences = [sequence]
    
        for i in site_start:site_end
            mutations = filter(x-> x[1] != sequence[i], letters)
            for x in mutations
                temp_sequence = collect(sequence)
                temp_sequence[i] = x[1]
                push!(scrambled_sequences, string(temp_sequence...))
            end
        end
    elif alph_type == "Numeric"
        letters = collect(1:4)
        scrambled_sequences = [sequence]
        for i in site_start:site_end
            mutations = filter(x-> x[1] != sequence[i], letters)
            for x in mutations
                temp_sequence = deepcopy(sequence)
                temp_sequence[i] = x[1]
                push!(scrambled_sequences, temp_sequence)
            end
        end
    else:
        raise ValueError("Alphabet type has to be either \"DNA\" or \"Numbers\"")
    
    return scrambled_sequences



function site_single_mutations(sequence; alph_type::String="DNA")
    
    if alph_type == "DNA"
        letters = ["A", "C", "G", "T"]
        scrambled_sequences = [sequence]
    
        for i in 1:length(sequence)
            mutations = filter(x-> x[1] != sequence[i], letters)
            for x in mutations
                temp_sequence = collect(sequence)
                temp_sequence[i] = x[1]
                push!(scrambled_sequences, string(temp_sequence...))
            end
        end
    elseif alph_type == "Numeric"
        letters = collect(1:4)
        scrambled_sequences = [sequence]
        for i in 1:length(sequence)
            mutations = filter(x-> x[1] != sequence[i], letters)
            for x in mutations
                temp_sequence = deepcopy(sequence)
                temp_sequence[i] = x[1]
                push!(scrambled_sequences, temp_sequence)
            end
        end
    else
        throw(ArgumentError("Alphabet type has to be either \"DNA\" or \"Numbers\""))
    end
    
    return scrambled_sequences
end
"""