import numpy as np
import copy
import random


def scramble(sequence, attempts, preserve_content=True):
    """
    Scramble letters of a given sequence. Returns most distant sequence.
    
    Parameters
    ----------
    sequence : string
        
    attempts : int
        Number of scrambles which are created. Most dissimilar one is chosen.
    preserve_content : bool
        If True, shuffles the existing sequence. If Flase, a completely arbitrary sequence is created.
    Returns
    -------
    scrambles[np.argmax(distances)] : string
        Most distant scrambled sequence
    """
    
    if type(attempts) != int:
        raise ValueError("attempts has to be an integer.")
    
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
        distances[i] = np.sum([x == y for (x, y) in zip(sequence, scrambles[i])])

    return scrambles[np.argmax(distances)]
    
    

def create_scrambles(sequence, windowsize, overlap, attempts, preserve_content=True):
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
    preserve_content : bool
        If True, shuffles the existing sequence. If Flase, a completely arbitrary sequence is created.
        
    Returns
    -------
    scrambled_sequences : list
        List of scrambled sequences
    """
    
    if overlap >= windowsize:
        raise ValueError("Overlap cannot be equal to or bigger than windowsize.")
        
    l_sequence = len(sequence)
    n_scrambles = (l_sequence - windowsize) / (windowsize - overlap)
    
    # Test if scrambles can created throughout whole site
    if not n_scrambles.is_integer():
        raise ValueError("Cannot make scrambles with windowsize {} and overlap {}.".format(windowsize, overlap))
    
    scrambled_sequences = np.empty(int(n_scrambles) + 1, dtype='<U160')
    scrambled_sequences[0] = sequence
    
    indices = lambda i: slice((windowsize-overlap) * i, (windowsize-overlap) * i + windowsize)
    
    for i in range(int(n_scrambles)+1):
        temp_sequence = list(sequence)
        temp_sequence[indices(i)] = scramble(sequence, attempts, preserve_content=True)
        scrambled_sequences[i] = "".join(temp_sequence)
    
    return scrambled_sequences