# Functions for generating random sequences and energy matrices

import random
import string
import pandas as pd
import numpy as np

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