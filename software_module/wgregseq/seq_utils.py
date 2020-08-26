
def scramble(sequence, windowsize, overlap):
    """
    Scramble letters in a window of given sequence. 
    
    Parameters
    ----------
    sequence : string
        
    windowsize : int
        Size of the window within letters are scrambled
    overlap : int
        Number of letters that each window overlaps with the neighbors
        
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
        raise ValueError("Cannot make scrambles with windowsize $windowsize and overlap $overlap.")
    """
    scrambled_sequences = [sequence]
    indices(i) = (windowsize-overlap)*(i-1)+1 : (windowsize-overlap)*(i-1) + windowsize
    for i in 1:Int(n_scrambles)+1
        temp_sequence = collect(sequence)
        temp_sequence[indices(i)] = shuffle(temp_sequence[indices(i)])
        push!(scrambled_sequences, string(temp_sequence...))
    
    return scrambled_sequences
    """