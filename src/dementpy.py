"""
-------------------------------------------------------------------------------
      DEMENTpy--Decomposition Model of Enzymatic Traits in Python
                                 Bin Wang
                          Email: wbwenwu@gmail.com
                          Twitter: @bioatmo_sphere
-------------------------------------------------------------------------------
"""
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
                                DEMENTpy 
            (DEcomposition Model of Enzymatic Traits in Python)
    ---------------------------------------------------------------------------       
    """)
    
     
    #...Obtain the command line arguments
    input_folder  = sys.argv[1]   # input folder name
    output_folder = sys.argv[2]   # output folder name
    outname       = sys.argv[3]   # output file name and seed of Pseudo-RNG
    site          = sys.argv[4]   # switch-site name
    
    #...Set up the working directory
    os.chdir('../'+input_folder)

    #...seed the pseudorandom number generator
    np.random.seed(int(outname[:]))

    #...a few system constants
    runtime    = pd.read_csv('runtime.txt',header=None,index_col=0,sep='\t')
    pulse      = int(runtime.loc['pulse',1])         # number of pulses
    cycle      = int(runtime.loc['end_time',1])      # number of time steps in each pulse
    interval   = int(runtime.loc['interval',1])      # interval of time step to record outputs
    switch     = int(runtime.loc['switch',1])        # the pulse from which inputs to be changed onwards  
    mode       = int(runtime.loc['dispersal',1])     # 0:'default' or 1:'dispersal'

    #...Initialize data by calling the Function: Initialize_Data()
    data_start = initialize_data(runtime,'.') # starting site

    #...Prepare for output by creating an instance of the Output class
    output_init = Output(runtime,data_start)

    #...Create an instance of the Grid class
    ecosystem = Grid(runtime,data_start)

    #...Run the model
    for p in range(pulse):
        
        for i in range(cycle):
        
            # substrates degradation
            ecosystem.degradation(i)
        
            # monomers uptake
            ecosystem.uptake(i)
        
            # microbial metabolism
            ecosystem.metabolism(i)
        
            # microbial death
            ecosystem.mortality(i)
        
            # microbial reproduction and dispersal
            ecosystem.reproduction(i)
        
            # output basic data using the "output" method in the Output class
            if i == 0:
                output_init.output(ecosystem, p, i)  # day 1
            elif i%interval==interval-1:
                output_init.output(ecosystem, p, i)  # interval
            # output microibal allocation data using the "microbes_tradeoff" method
            output_init.microbes_tradeoff(ecosystem, p, i)
            # output microbial mass of every iteration using the "microbes_abundance" method
            #output_init.microbes_abundance(ecosystem, p, i)

        # re-initialize microbial community in each new pulse
        # e.g., if switch=3 (i.e., switch after a 3-year run)
        if p < switch - 1 or p > switch - 1:
            # stick to the starting site
            ecosystem.reinitialization(data_start, data_start['Microbes_pp'], output_init, mode, p, switch)
        else:
            # switch to data of another site
            data_switch = initialize_data(runtime, site) # switching site
            ecosystem.reinitialization(data_switch,data_start['Microbes_pp'], output_init, mode, p, switch)
    
    #...export the Output_init object to the output_folder using the export() funtion in the utility module 
    os.chdir('../'+output_folder)
    export(output_init, site, outname)
    
main()