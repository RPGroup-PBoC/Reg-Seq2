import sys

import numpy as np
import pandas as pd

from Bio.Restriction import *
from Bio.Seq import Seq

from .seq_utils import _check_sequence_list
from .utils import isint

def mutation_coverage(wildtype_seq, mutants_list, site_start=0, site_end=None):
    
    if not isint(site_start):
        raise TypeError("`site_start` is of type {} but has to be integer valued.".format(type(site_start)))
        
    if not (isint(site_end) or site_end == None):
        raise TypeError("`site_end` is of type {} but has to be integer valued.".format(type(site_end)))
    
    if site_end==None:
        wildtype_seq = wildtype_seq[site_start:]
        mutations = np.zeros([len(mutants_list), len(wildtype_seq)])
        wildtype_list = list(wildtype_seq)
        
        for i, sequence in enumerate(mutants_list):
            mutations[i, :] = [x != y for (x, y) in zip(list(sequence[site_start:]), wildtype_list)]
    else:
        wildtype_seq = wildtype_seq[site_start:site_end]
        mutations = np.zeros([len(mutants_list), len(wildtype_seq)])
        wildtype_list = list(wildtype_seq)
        for i, sequence in enumerate(mutants_list):
            mutations[i, :] = [x != y for (x, y) in zip(list(sequence[site_start:site_end]), wildtype_list)]
           
    return np.mean(mutations, axis=0)


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
    if type(sequence_list) == pd.core.series.Series:
        sequence_list = copy.deepcopy(sequence_list.values)
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