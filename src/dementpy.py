"""
-------------------------------------------------------------------------------
      DEMENTpy--Decomposition Model of Enzymatic Traits in Python,v1.0
                              Bin Wang, Ph.D.
           Department of Ecology and Evolutionary Biology
                  University of California Irvine
                  Emails: wbwenwu@gmail.com or bwang7@uci.edu
                      Twitter: @bio_atmosphere
-------------------------------------------------------------------------------
"""

#import time
import os
import sys
import pandas as pd
import numpy  as np

from initialization import initialize_data
from grid import Grid
from output import Output
from utility import export

def main():
    
    
    print("""
    ---------------------------------------------------------------------------
         DEMENTpy (DEcomposition Model of Enzymatic Traits in Python)
                               Version 1.0
               Department of Ecology and Evolutionary Biology
                     University of California Irvine
    ---------------------------------------------------------------------------       
    """)
    
     
    #...Obtain the command line arguments of runtime file name and output object name
    runtime = sys.argv[1]
    outname = sys.argv[2]
    
    #...Set up the working directory
    os.chdir('../input')
    
    #...a few system constants
    runtime = pd.read_csv(runtime,header=None,index_col=0,sep='\t')
    pulse = int(runtime.loc['pulse',1])         # number of pulses
    cycle = int(runtime.loc['end_time',1])      # number of time steps in each pulse
    interval = int(runtime.loc['interval',1])   # interval of time step to record outputs
    mic_reinit = runtime.loc['mic_reinit',1]    # indicate whether or not reinitialize microbial community on the spatial grid in a new pulse
    
    #...grow a seed of random number generator
    np.random.seed(2)
    
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
            #if pulse == 1:
            #    continue
            
            # output every day's microbial mass 
            Output_init.microbes_df(Ecosystem,i)
            
            # re-initialize microbial community
            if i == (p+1)*cycle-1:
                Ecosystem.repopulation(Output_init,i,mic_reinit)
        
        
    #...export the Output_init object using the export() funtion in the utility module 
    os.chdir('../output')
    export(Output_init,outname)
    
    #print out time used
    #print('   ',"Cumulative time:", time.time() - start_time)
    
main()
