"""
this module contains a couple of functions that facilitate the calculation, including:
    1) LHS():    perform 1-D latin hypercube sampling
    2) expand(): expand a community to the whole spatial grid;
    3) export(): export the output object as a file to the local disk
"""

import numpy as np
import pandas as pd
import pickle
from scipy.stats import distributions


def LHS(n,loc,scale,dist):
    
    """
    Inputs:
        n:     integer;size of desired sampling
        loc:   scalar;para 1 for specifying a desired probability distribution
        scale: scalar;para 2 for specifying a distribution
        dist:  string; either 'uniform' or 'normal'
    """
    
    lower_limits  = np.arange(0,n)/n
    higher_limits = np.arange(1,n+1)/n
    
    points = np.random.uniform(low=lower_limits,high=higher_limits,size=n)
    np.random.shuffle(points)
    
    if dist == 'uniform':
        rv = distributions.uniform(loc=loc,scale=scale)
        
    elif dist == 'normal':
        rv = distributions.norm(loc=loc,scale=scale)
    
    lhs = rv.ppf(points)
    
    return lhs


def expand(df,gridsize):
    """
    Put data on a spatial grid
    """
    df_expanded = pd.concat([df]*gridsize)
    
    return df_expanded
    


def export(output):
    """
    
    """
    with open("Output.pickle", "wb") as f:
        pickle.dump(output, f, pickle.HIGHEST_PROTOCOL)
    
