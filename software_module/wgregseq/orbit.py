import sys
from contextlib import suppress
import numpy as np
import pandas as pd
import copy
import random
import numba
import itertools

from scipy.special import factorial

import pathlib

from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.SeqIO import parse

from .utils import isint, choose_dict

import warnings

def get_replichore(pos, ori = 3923882.5, ter = 1590250.5 ):
    
    """
    Determine the replichore of a bacterial chromosome for a certain position. Requires origin and terminus positions. Assumes E. coli like organization.
    
    pos : int
        Genomic coordinate of interest.
    ori : float
        Genomic coordinate of the origin of replication. 
    ter : float
        Genomic coordinate of the replication terminus.
    """
    
    if((pos<0)| (pos>4641652)):
        raise TypeError("position must be within genome.")
    
    if((pos > ori) | (pos<ter)):
       rep = 1
    elif((pos<ori) & (pos>ter)):
       rep = 2
    
    return rep
       
    
def get_target_oligo(left_pos, right_pos, genome, homology = 90, attB_dir = '+', attB_fwd_seq = 'ggcttgtcgacgacggcggtctccgtcgtcaggatcat', verbose = False):
    """
    Given a set of parameters, get an ORBIT oligo that targets the lagging strand. 
    Left and right positions are absolute genomic coordinates that specify the final nucleotides to keep unmodified in the genome, 
    everything in between will be replaced by attB. In other words the left position nucleotide is the final nt before attB in the oligo.
    The right position nt is the first nt after attB in the oligo.
    
    This function determines the lagging strand by calling `get_replichore()` on the left_pos.
    Typically attB_dir should be set to the same direction as the gene of interest, such that the integrating plasmid will insert with payload facing downstream.
    attB_fwd_seq can be modified, and the total homology can be modified, but should be an even number since homology arms are symmetric. 
    
    Verbose prints helpful statements for testing functionality.
    
    Parameters
    -----------------
    left_pos : int
        Left genomic coordinate of desired attB insertion. attB is added immediately after this nt.
    right_pos : int
        Right genomic coordinate of desired attB insertion. attB is added immediately before this nt.
    genome : str
        Genome as a string.
    homology : int (even)
        Total homology length desired for oligo. Arm length = homology / 2.
    attB_dir : chr ('+' or '-')
        Desired direction of attB  based on genomic strand. Typically same direction as gene.
    attB_fwd_seq : str
        Sequence of attB to insert between homology arms.
    verbose : bool
        If true, prints details about genomic positions and replichore.
    Returns
    ---------------
    oligo : str
        Targeting oligo against lagging strand, including the attB sequence in the correct orientation.
    """
    
    left_pos = int(left_pos)
    
    right_pos = int(right_pos)
    
    # Arm length is 1/2 total homology. Arms are symmetric
    arm_len = int(homology / 2)
    
    # Arms from genome string. Note 0 indexing of string vs. 1 indexing of genomic coordinates.
    # As written, should be inclusive.
    left_arm = genome[(left_pos - arm_len):left_pos]
    
    right_arm = genome[(right_pos - 1):(right_pos - 1 + arm_len)]

    # Generate attB reverse sequence
    seq_attB = Seq(attB_fwd_seq)
    attB_rev_seq = str(seq_attB.reverse_complement())
    
    # Replichore 1
    if get_replichore(left_pos) == 1:
        
        rep = 1
        
        # Reverse complement replichore 1 sequences.
        left_arm_seq = Seq(left_arm)
        left_arm_prime = str(left_arm_seq.reverse_complement())
        
        right_arm_seq = Seq(right_arm)
        right_arm_prime = str(right_arm_seq.reverse_complement())
        
        # Determine attB direction and paste fwd/rev seq accordingly
        if attB_dir == '+':
            
            oligo = right_arm_prime + attB_rev_seq + left_arm_prime
            
        elif attB_dir == '-':
            
            oligo = right_arm_prime + attB_fwd_seq + left_arm_prime
    
    # Replichore 2
    elif get_replichore(left_pos) == 2:
        
        rep = 2
        
        # '+' arm sequence used. Determine attB direction and paste accordingly.
        if attB_dir == '+':
            
            oligo = left_arm + attB_fwd_seq + right_arm
        
        elif attB_dir == '-':
            
            oligo = left_arm + attB_rev_seq + right_arm    
            
    # Verbose print statements
    if verbose:
        
        print('left_arm_coord = ', left_pos - arm_len,' : ', left_pos)
        print('right_arm_coord = ', right_pos - 1, ' : ', right_pos -1 + arm_len)
        print('Replichore = ', rep)
    
    return oligo

def get_target_oligo_df(df, left_pos_col, right_pos_col, dir_col, genome, homology = 90, attB_fwd_seq = 'ggcttgtcgacgacggcggtctccgtcgtcaggatcat'):
    
    """
    Apply get_target_oligo to a dataframe of genomic coordinates and directions. Iterates through df rows calling get_target_oligo given the parameters specified in each column.
    
    Given a set of parameters, get an ORBIT oligo that targets the lagging strand. 
    Left and right positions are absolute genomic coordinates that specify the final nucleotides to keep unmodified in the genome, 
    everything in between will be replaced by attB. In other words the left position nucleotide is the final nt before attB in the oligo.
    The right position nt is the first nt after attB in the oligo.
    
    This function determines the lagging strand by calling `get_replichore()` on the left_pos.
    Typically attB_dir should be set to the same direction as the gene of interest, such that the integrating plasmid will insert with payload facing downstream.
    attB_fwd_seq can be modified, and the total homology can be modified, but should be an even number since homology arms are symmetric. 
        
    Parameters
    -----------------
    df : pd.DataFrame
        Pandas dataframe containing the required genomic coordinates, and gene directions.
    left_pos_col : str
        Column name of left genomic coordinate of desired attB insertion. attB is added immediately after this nt. 
    right_pos_col : str
        Column name of right genomic coordinate of desired attB insertion. attB is added immediately after this nt. 
    dir_col : str
        Column name of desired direction of attB based on genomic strand. Typically same direction as gene.
    genome : str
        Genome as a string.
    homology : int (even)
        Total homology length desired for oligo. Arm length = homology / 2.   
    attB_fwd_seq : str
        Sequence of attB to insert between homology arms.
    verbose : bool
        If true, prints details about genomic positions and replichore.
    Returns
    ---------------
    df_results : pd.DataFrame
        Adds column 'oligo' to input df. 'oligo' contains a string of the targeting oligo sequence against lagging strand, including the attB sequence in the correct orientation.
    """
    
    df_tmp = pd.DataFrame()
    df_results = pd.DataFrame()
    
    for i,row in df.iterrows():
        
        left_pos = row[left_pos_col]
        right_pos = row[right_pos_col]
        attB_dir = row[dir_col]
        
        oligo = get_target_oligo(left_pos, right_pos, genome, homology, attB_dir, attB_fwd_seq)

        df_tmp = df.iloc[[i],:]
        
        df_tmp['oligo'] = oligo
        
        df_results = pd.concat([df_results,df_tmp])
    
    return df_results