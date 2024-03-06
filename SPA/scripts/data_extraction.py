"""
Script of extracting data from pickle files

Bin Wang
11/25/2023

"""

#import numpy as np
import pandas as pd
import pickle
import sys
import os
import glob
# import output module
os.chdir('../../src')
sys.path.append(os.getcwd())
import output

# get command line arguments
site = sys.argv[1]    # base site name
key  = sys.argv[2]    # target site

# back to working directory
os.chdir('../output_'+site)

# define a function of extracting data from files in .pickle
def get_pickled_data(key):
    """
    Extract data from .pickle files.

    """

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

#
def community_traits(data):
    """
    Calculate community-level traits


    Returns:
        enzyme:
        drought:
        thermal:
        stress:
    """

    Relative_mass = data.MicrobesSeries.div(data.MicrobesSeries.sum(axis=0),axis=1)

    enzyme_trait  = (data.Microbial_traits['Enz_Induci_Cost'] + data.Microbial_traits['Enz_Consti_Cost']) * data.Microbial_traits['Enz_Gene']
    community_enzyme = Relative_mass.mul(enzyme_trait,axis=0).sum(axis=0)

    drought_tol   = data.Microbial_traits['Drought_tolerance']
    community_drought = Relative_mass.mul(drought_tol,axis=0).sum(axis=0)

    thermal_tol   = data.Microbial_traits['Thermal_tolerance']
    community_thermal = Relative_mass.mul(thermal_tol,axis=0).sum(axis=0)

    stress_tol   = community_drought + community_thermal

    return community_enzyme, community_drought, community_thermal,  stress_tol

####  collect output files
filelist, datalist = get_pickled_data(key)

#### get community tratis
# call function community_trait()
comm_trait_enz = pd.concat([community_traits(data)[0] for data in datalist], axis=1, sort=False)
comm_trait_drt = pd.concat([community_traits(data)[1] for data in datalist], axis=1, sort=False)
comm_trait_thm = pd.concat([community_traits(data)[2] for data in datalist], axis=1, sort=False)
comm_trait_str = pd.concat([community_traits(data)[3] for data in datalist], axis=1, sort=False)
#export traits to independent csv files
comm_trait_enz.to_csv('data/' + 'Enz_' + site +'_'+ key + '.csv')
comm_trait_drt.to_csv('data/' + 'Drt_' + site +'_'+ key + '.csv')
comm_trait_thm.to_csv('data/' + 'Thm_' + site +'_'+ key + '.csv')
comm_trait_str.to_csv('data/' + 'Str_' + site +'_'+ key + '.csv')

#### substrates
## sub-specific mass
subs = pd.concat([data.SubstratesSeries for data in datalist], axis=1, sort=False)
# export to csv
subs.to_csv('data/' + 'Subs_' + site +'_'+ key + '.csv')
## total sub mass
#sub = pd.concat([data.SubstratesSeries.sum(axis=0) for data in datalist], axis=1, sort=False)
#sub.to_csv('data/' + 'Sub_' + site +'_'+ key + '.csv')

##### microbes
microbes = pd.concat([data.MicrobesSeries for data in datalist], axis=1, sort=False)
microbes.to_csv('data/' + 'Mic_' + site +'_'+ key + '.csv')

#### metabolism
# taxon-specific enzyme
Enzyme = pd.concat([data.Enzyme_TaxonSeries for data in datalist], axis=1, sort=False)
Enzyme.to_csv('data/' + 'Enzyme_' + site +'_'+ key + '.csv')
# taxon-specific osmolyte
Osmolyte = pd.concat([data.Osmolyte_TaxonSeries for data in datalist], axis=1, sort=False)
Osmolyte.to_csv('data/' + 'Osmolyte_' + site +'_'+ key + '.csv')
# taxon-specific hsp
Hsp = pd.concat([data.Hsp_TaxonSeries for data in datalist], axis=1, sort=False)
Hsp.to_csv('data/' + 'Hsp_' + site +'_'+ key + '.csv')
# taxon-specific growth yield
Yield = pd.concat([data.Growth_yield for data in datalist], axis=1, sort=False)
Yield.to_csv('data/' + 'Yield_' + site +'_'+ key + '.csv')
