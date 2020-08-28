import numpy as np
import pandas as pd
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
        If True, shuffles the existing sequence. If False, a completely arbitrary sequence is created.
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
    n_scrambles = (l_sequence - windowsize) / (windowsize - overlap) + 1
    
    # Test if scrambles can created throughout whole site
    if not n_scrambles.is_integer():
        raise ValueError("Cannot make scrambles with windowsize {} and overlap {}.".format(windowsize, overlap))
    
    scrambled_sequences = np.empty(int(n_scrambles) + 1, dtype='<U160')
    scrambled_sequences[0] = sequence
    
    indices = lambda i: slice((windowsize-overlap) * i, (windowsize-overlap) * i + windowsize)
    
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
    scrambled_sequences : list
        List of scrambled sequences
    """
    
    if overlap >= windowsize:
        raise ValueError("Overlap cannot be equal to or bigger than windowsize.")
        
    l_sequence = len(sequence)
    n_scrambles = (l_sequence - windowsize) / (windowsize - overlap) + 1
    
    # Test if scrambles can created throughout whole site
    if not n_scrambles.is_integer():
        raise ValueError("Cannot make scrambles with windowsize {} and overlap {}.".format(windowsize, overlap))
    
    scrambled_sequences = np.empty(int(n_scrambles) + 1, dtype='<U160')
    scrambled_sequences[0] = sequence
    
    indices = lambda i: slice((windowsize-overlap) * i, (windowsize-overlap) * i + windowsize)
    
    for i in range(int(n_scrambles)):
        #print((windowsize-overlap)*i)
        #print((windowsize-overlap)*i + windowsize)
        temp_sequence = list(sequence)
        temp_sequence[indices(i)] = scramble(sequence[indices(i)], attempts, preserve_content=True)
        scrambled_sequences[i+1] = "".join(temp_sequence)
    
    start_pos = np.arange(0,int(n_scrambles),1)*(windowsize-overlap)
    stop_pos = np.arange(0,int(n_scrambles),1)*(windowsize-overlap)+windowsize
    
    #print(start_pos)
    #print(stop_pos)
    #print(scrambled_sequences)
    
    scramble_df = pd.DataFrame({'start_pos':start_pos, 'stop_pos':stop_pos, 'sequence':scrambled_sequences[1:]})
    
    #return scrambled_sequences
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