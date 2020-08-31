import numpy as np
import pandas as pd
from .gen_utils import isint


def gen_emat_rand(site_size, mean=1, sd=1):
    """
    Generate a random energy matrix for a defined sequence length. Arbitrary values for each possible base, normally distributed around mean 1 with standard deviation 1.
    
    Parameters
    ----------
    site_size : int
        Length of the sequence to generate the energy matrix for, in bp.
    mean : float
        Mean of entries in energy matrix.
    sd : float
        Standard deviation of entries in energy matrix.
    Returns
    ----------
    energy_matrix : np.array
    """
    # Check argument types
    if not isint(site_size):
        raise ValueError("`site_size` has to be an integer.")
    else:
        # If type float, change to int
        site_size = int(site_size)
        
    energy_matrix = np.random.normal(mean, sd, (site_size, 4))
    return energy_matrix



def gen_emat_single_site(
    seq, 
    site_start, 
    site_size, 
    site_mean=1, 
    site_sd=1, 
    background_mean=0, 
    background_sd=0):
    """
    Generate energy matrix for sequence with one site. 
    Mean and sd values can be set for site and non site positions.
    WT sequence is set to zero.
    
    Parameters
    ----------
    seq : string
        Sequence. Used to set entries for wild type binding site to zero.
    site_start : int
        First base of binding site
    site_size : int
        Length of binding site.
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
    
    # Check argument types
    if not isint(site_start):
        raise ValueError("`site_start` has to be an integer.")
    else:
        # If type float, change to int
        site_start = int(site_start)
        
    if not isint(site_size):
        raise ValueError("`site_size` has to be an integer.")
    else:
        # If type float, change to int
        site_size = int(site_size)
        
    if not isinstance(site_mean, int) or isinstance(site_mean, float):
        raise ValueError("`site_mean` has to be an integer or float.")
    
    if not isinstance(site_sd, int) or isinstance(site_sd, float):
        raise ValueError("`site_sd` has to be an integer or float.")
        
    if not isinstance(background_mean, int) or isinstance(background_mean, float):
        raise ValueError("`background_mean` has to be an integer or float.")
        
    if not isinstance(background_sd, int) or isinstance(background_sd, float):
        raise ValueError("`background_sd` has to be an integer or float.")
        
        
    # Set background values
    seq_emat = np.random.normal(background_mean, background_sd, (len(seq), 4))
    
    # Set site values
    seq_emat[site_start:(site_start + site_size), :] = np.random.normal(site_mean, site_sd,(site_size, 4))
    
    # Convert np.array to pd.DataFrame
    seq_emat = pd.DataFrame(data=seq_emat, columns=('A','T','C','G'))
    
    # Set WT values = 0
    for ind,char in enumerate(seq, start=0):
        seq_emat.iloc[ind][char] = 0
    
    return seq_emat 


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
        
    return sum(mat_vals) 