"""
Script of extracting data from pickle files
By Bin Wang
"""

#import numpy as np
import pandas as pd
import pickle
import sys
import os
import glob
import output

# define a function of extracting data from files in .pickle
def get_pickled_data(key):
    datalist = []
    #single run
    filelist = glob.glob(key+'_'+'20201'+'.pickle')
    
    ## ensemble runs--20
    filelist_19   = glob.glob(key + '_' + '2020' + '[1-9].pickle')
    filelist_1019 = glob.glob(key + '_' + '2020' + '1[0-9].pickle')
    filelist_2020 = glob.glob(key + '_' + '2020' + '20.pickle')
    filelist = filelist_19 + filelist_1019 + filelist_2020

    #filelist_2029 = glob.glob(key+'2[0-9].pickle')
    #filelist_3039 = glob.glob(key+'3[0-9].pickle')
    #filelist_4040 = glob.glob(key+'40.pickle')
    #filelist = filelist_19 + filelist_1019 + filelist_2029 + filelist_3039 + filelist_4040

    filelist.sort(reverse=False)
    for file in filelist:
        with open(file,"rb") as f:
            data = pickle.load(f)
        datalist.append(data)

    return filelist, datalist

def community_drought(data):
    """
    Calculate community-level drought tolerance
    """

    Relative_mass = data.MicrobesSeries.div(data.MicrobesSeries.sum(axis=0),axis=1)
    # enzyme_trait  = data.Microbial_traits['Enz_Induci_Cost'] * data.Microbial_traits['Enz_Gene']
    #enzyme_trait  = (data.Microbial_traits['Enz_Induci_Cost'] + data.Microbial_traits['Enz_Consti_Cost']) * data.Microbial_traits['Enz_Gene']
    drought_tol   = data.Microbial_traits['Drought_tolerance']
    community_drought = Relative_mass.mul(drought_tol,axis=0).sum(axis=0)

    return community_drought

def community_enzyme(data):
    """
    Calculate community-level drought tolerance
    """

    Relative_mass = data.MicrobesSeries.div(data.MicrobesSeries.sum(axis=0),axis=1)
    # enzyme_trait  = data.Microbial_traits['Enz_Induci_Cost'] * data.Microbial_traits['Enz_Gene']
    enzyme_trait  = (data.Microbial_traits['Enz_Induci_Cost'] + data.Microbial_traits['Enz_Consti_Cost']) * data.Microbial_traits['Enz_Gene']
    community_enzyme = Relative_mass.mul(enzyme_trait,axis=0).sum(axis=0)

    return community_enzyme

site = sys.argv[1]    # base site name
key  = sys.argv[2]    # target site

os.chdir('../output_'+site)

filelist, datalist = get_pickled_data(key)
## sub-specific mass
#sub = pd.concat([data.SubstratesSeries for data in datalist], axis=1, sort=False)
## total mass
sub = pd.concat([data.SubstratesSeries.sum(axis=0) for data in datalist], axis=1, sort=False)
# export to csv
sub.to_csv('data/' + 'Sub_' + site +'_'+ key + '.csv')

# microbes
microbes = pd.concat([data.MicrobesSeries for data in datalist], axis=1, sort=False)
microbes.to_csv('data/' + 'Mic_' + site +'_'+ key + '.csv')

# taxon-specific enzyme
Enzyme = pd.concat([data.Enzyme_TaxonSeries for data in datalist], axis=1, sort=False)
Enzyme.to_csv('data/' + 'Enzyme_' + site +'_'+ key + '.csv')

# taxon-specific osmolyte
Osmolyte = pd.concat([data.Osmolyte_TaxonSeries for data in datalist], axis=1, sort=False)
Osmolyte.to_csv('data/' + 'Osmolyte_' + site +'_'+ key + '.csv')

# taxon-specific growth yield
Yield = pd.concat([data.Growth_yield for data in datalist], axis=1, sort=False)
Yield.to_csv('data/' + 'Yield_' + site +'_'+ key + '.csv')
