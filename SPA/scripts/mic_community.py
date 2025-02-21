"""
Script of extracting data of microbial community based on taxon-specific mass

By Bin Wang
September 2020
"""

import numpy as np
import pandas as pd
import pickle
import sys
import os
import glob
import output

# define a function of extracting data from files in .pickle 
def get_pickled_data(key):
    datalist = []
    filelist = glob.glob(key+'_'+'20201'+'.pickle')
    filelist.sort(reverse=False)
    for file in filelist:
        with open(file,"rb") as f:
            data = pickle.load(f)
        datalist.append(data)

    return filelist, datalist



site = sys.argv[1]    # site name
key = sys.argv[2]     # 


os.chdir('../output_'+site)

filelist, datalist = get_pickled_data(key)
microbes = pd.concat([data.MicrobesSeries for data in datalist], axis=1, sort=False)

# rename column names with file names
#filelist_new = [int(item[:-7]) for item in filelist]
#sub.columns = filelist_new
#filelist_sorted = sorted(filelist_new,reverse=False)
#sub = sub[filelist_sorted]

# export to csv
microbes.to_csv('data/' + 'Mic_' + site +'_'+ key + '.csv')
