"""
This module, utility.py, contains three functions facilitating the calculation.

    1) LHS():    perform 1-D latin hypercube sampling
    2) expand(): put a community on the spatial grid
    3) export(): export the output object to the local disk
"""

import numpy as np
import pandas as pd
import pickle
from scipy.stats import distributions


def LHS(n,loc,upc,dist):
    """
    Latin hypercube sampling.

    Parameters:
        n:    integer; size of desired sampling
        loc:  scalar; lower bound of desired distribution
        upc:  scalar; upper bound of desired distribution
        dist: string; either 'uniform' or 'normal'
    Returns:
        lhs: 1D array
    """
    
    lower_limits  = np.arange(0,n)/n
    higher_limits = np.arange(1,n+1)/n
    
    points = np.random.uniform(low=lower_limits,high=higher_limits,size=n)
    np.random.shuffle(points)
    
    scale = upc - loc
    if dist == 'uniform':
        rv = distributions.uniform(loc=loc, scale=scale)  
    elif dist == 'normal':
        rv = distributions.norm(loc=loc, scale=scale)
    
    lhs = rv.ppf(points)
    
    return lhs


def expand(df,gridsize):
    """
    Put data (df/series) on a spatial grid.
    
    Parameters:
         df:       dataframe/series
         gridsize: integer
    Return:
         df_expanded: dataframe/series;
    """
    df_expanded = pd.concat([df]*gridsize)
    
    return df_expanded
    

def export(output,name):
    """
    Save the output object as a .pickle file.
    
    Parameters:
        output: object of the class, Output
        name:   naming output file
    """
    
    with open(str(name) + ".pickle", "wb") as f:
        pickle.dump(output, f, pickle.HIGHEST_PROTOCOL)