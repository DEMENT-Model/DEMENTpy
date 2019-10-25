"""
-------------------------------------------------------------------------------
    DEMENTpy:Decomposition Model of Enzymatic Traits in Python,v1.0
                   Bin Wang
    Department of Ecology and Evolutionary Biology
           University of California Irvine
                 2019-09-29
-------------------------------------------------------------------------------
Direct correspondance to: 
    Bin Wang (wbwenwu@gmail.com/bwang7@uci.edu)
"""


import os
import sys
#import time
import pandas as pd
import numpy  as np

from initialization import initialize_data
from grid import Grid
from output import Output
from utility import export

def main():
    
    
    print("""
    ----------------------------------------------
    DEMENTpy (DEcomposition Model of Enzymatic Traits in Python) Version 1.0
    Department of Ecology and Evolutionary Biology
    University of California Irvine
    ---------------------------------------------       
    """)
    
    
    np.random.seed(2)
     
    #...obtain the command line argument of runtime_file.txt
    filename = sys.argv[1]  
    runtime = pd.read_csv(filename,header=None,index_col=0,sep='\t')
    
    #... a few system constants
    pulse = int(runtime.loc['pulse',1])
    cycle = int(runtime.loc['end_time',1])
    interval = int(runtime.loc['interval',1])
    mic_reinit = runtime.loc['Mic_reinit',1]
    
    #...set up the working directory
    os.chdir('../input_output')
    
    #...start timing
    #start_time = time.time()
    
    #...Initialize data by calling the Function: Initialize_Data()
    data_initialization = initialize_data(runtime)
    
    #...Prepare for output by creating an instance of the Output class
    Output_init = Output(runtime,data_initialization)
    
    #...Create an instance of the Grid class
    Ecosystem = Grid(runtime,data_initialization)
    
    #...Run the model
    for p in range(pulse):
        
        for i in range(p*cycle, (p+1)*cycle):
        
            # substrates degradation
            Ecosystem.degradation(p,i)
        
            # monomers uptake
            Ecosystem.uptake(p,i)
        
            # microbial metabolism
            Ecosystem.metabolism(i)
        
            # microbial death
            Ecosystem.mortality(i)
        
            # microbial reproduction and dispersal
            Ecosystem.reproduction(i)
        
            # output data using the output method in the Output class
            if i == 0:
                Output_init.output(Ecosystem,i)  # day 1
            elif i%interval==interval-1:
                Output_init.output(Ecosystem,i)  # interval
            
            # if only 1 pusle, skip all following lines within this loop
            if pulse == 1:
                continue
            
            # output every day's microbial mass 
            Output_init.microbes_df(Ecosystem,i)
            
            # re-initialize microbial community
            if i == (p+1)*cycle-1:
                Ecosystem.repopulation(Output_init,i,mic_reinit)
        
        
    #...use the export() funtion from the utility module
    export(Output_init)
    
    #print out time used
    #print('   ',"Cumulative time:", time.time() - start_time)
    
main()
