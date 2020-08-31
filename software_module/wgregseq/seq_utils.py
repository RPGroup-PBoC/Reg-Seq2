import numpy as np
import pandas as pd
import copy
import random


def gen_rand_seq(length):
    """
    Generate a random DNA sequence of defined length
    
    Parameters
    ----------
    length : int
    
    Returns
    ---------
    seq : string
         Random DNA sequence
    """

    nt = ['A','T','C','G'] #nucleotides
    seq = ''.join(random.choice(nt) for i in range(length))
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
    scrambled_sequences : Pandas.DataFrame
        DataFrame of scrambled sequences
    """
    
    l_sequence = len(sequence)
    n_scrambles = (l_sequence - windowsize) / (windowsize - overlap) + 1
    
    # Test if scrambles can created throughout whole site
    if not n_scrambles.is_integer():
        raise ValueError("Cannot make scrambles with windowsize {} and overlap {}.".format(windowsize, overlap))
        
    scrambled_sequences = create_scrambles(sequence, windowsize, overlap, attempts, preserve_content=True)
    start_pos = np.arange(0,int(n_scrambles),1)*(windowsize-overlap)
    stop_pos = np.arange(0,int(n_scrambles),1)*(windowsize-overlap)+windowsize
    
    scramble_df = pd.DataFrame({'start_pos':start_pos, 'stop_pos':stop_pos, 'sequence':scrambled_sequences[1:]})
    scramble_df['center_pos'] = scramble_df[['start_pos','stop_pos']].mean(axis = 1)
    
    return scramble_df


def gen_emat_rand(site_size):
    """
    Generate a random energy matrix for a defined sequence length. Arbitrary values for each possible base.
    
    Parameters
    ----------
    site_size : int
        Length of the sequence to generate the energy matrix for, in bp.
    
    Returns
    ----------
    energy matrix : np.array
    """
    return(np.random.normal(1,1,(site_size,4)))

'''
# old version
def gen_emat_single_site(seq, site_start, site_size):
    """
    Generate energy matrix for sequence with one site. Outside of site matrix is zero.
    
    Parameters
    ----------
    seq : string
    
    site_start : int
    
    site_size : int
    
    Returns
    ---------
    seq_emat : np.array
    """
    
    seq_emat = np.zeros((len(seq),4))
    
    seq_emat[site_start:(site_start + site_size),:] = gen_emat_rand(site_size)
    
    return(seq_emat)
'''

def gen_emat_single_site(seq, site_start, site_size, site_mean = 1, site_sd = 1, background_mean = 0, background_sd = 0):
    """
    Generate energy matrix for sequence with one site. 
    Mean and sd values can be set for site and non site positions.
    WT sequence is set to zero.
    
    Parameters
    ----------
    seq : string
    
    site_start : int
    
    site_size : int
    
    site_mean: float
        mean energy for site mutations, for np.random.normal
        
    site_sd: float
        standard deviation of energy for site mutations, for np.random.normal
    
    background_mean: float
        mean energy for non site mutations, for np.random.normal
    
    background_sd: float
        standard deviation of energy for non site mutations, for np.random.normal
    
    Returns
    ---------
    seq_emat : pd.DataFrame
        generated energy matrix
    """
    
    # Set background values
    seq_emat = np.random.normal(background_mean,background_sd,(len(seq),4))
    
    # Set site values
    seq_emat[site_start:(site_start + site_size),:] = np.random.normal(site_mean, site_sd,(site_size,4))
    
    # Convert np.array to pd.DataFrame
    seq_emat = pd.DataFrame(data = seq_emat, columns = ('A','T','C','G'))
    
    # Set WT values = 0
    for ind,char in enumerate(seq, start = 0):
        seq_emat.iloc[ind][char]=0
    
    return(seq_emat)

def sum_emat(seq, emat):
    """
    Retrieve and sum the energy matrix values for a given sequence variant and matrix.
    
    Parameters
    ----------
    seq : string
    
    emat : pd.DataFrame with columns A T C G
    
    Returns
    ---------
    sum : float
    """
    
    mat_vals = []
    
    for ind,char in enumerate(seq, start = 0):
        #print(ind)
        #print(char)
        mat_vals.append(emat.iloc[ind][char])
        
    return(sum(mat_vals))


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