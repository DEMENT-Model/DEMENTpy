"""
-------------------------------------------------------------------------------
      DEMENTpy--Decomposition Model of Enzymatic Traits in Python
            by Bin Wang (wbwenwu@gmail.com|@bio_atmosphere)
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
    input_folder  = sys.argv[1]   # input folder name; e.g., a study site or a project name
    output_folder = sys.argv[2]   # output folder name
    outname       = sys.argv[3]   # output file name and seed of Pseudo-RNG
    
    site = input_folder        # no forcing and/or parameter switching
    if len(sys.argv) > 4:
        site      = sys.argv[4]   # switch-site name

    #...Set up the working directory
    os.chdir('../'+input_folder)

    #...seed the pseudorandom number generator
    np.random.seed(int(outname[:]))

    #...a few system constants
    runtime    = pd.read_csv('runtime.txt',header=None,index_col=0,sep='\t')
    pulse      = int(runtime.loc['pulse',1])         # number of pulses
    cycle      = int(runtime.loc['end_time',1])      # number of time steps in each pulse
    interval   = int(runtime.loc['interval',1])      # interval of time step to record outputs
    mode       = int(runtime.loc['dispersal',1])     # 0:'default' or 1:'dispersal'

    #...Initialize data by calling the Function: Initialize_Data()
    data_initialization = initialize_data(runtime,'.') # starting site

    #...Prepare for output by creating an instance of the Output class
    Output_init = Output(runtime,data_initialization)

    #...Create an instance of the Grid class
    Ecosystem = Grid(runtime,data_initialization)

    #...Run the model
    for p in range(pulse):
        
        for i in range(cycle):
        
            # substrates degradation
            Ecosystem.degradation(i)
        
            # monomers uptake
            Ecosystem.uptake(i)
        
            # microbial metabolism
            Ecosystem.metabolism(i)
        
            # microbial death
            Ecosystem.mortality(i)
        
            # microbial reproduction and dispersal
            Ecosystem.reproduction(i)
        
            # output basic data using the "output" method in the Output class
            if i == 0:
                Output_init.output(Ecosystem, p, i)  # day 1
            elif i%interval==interval-1:
                Output_init.output(Ecosystem, p, i)  # interval
            # output microibal allocation data using the "microbes_tradeoff" method
            Output_init.microbes_tradeoff(Ecosystem, p, i)
            # output microbial mass of every iteration using the "microbes_abundance" method
            Output_init.microbes_abundance(Ecosystem, p, i)

        # re-initialize microbial community in each new pulse
        if site == input_folder:
            Ecosystem.reinitialization(data_initialization, data_initialization['Microbes_pp'], Output_init, mode, p)
        else: 
            switch     = int(runtime.loc['switch',1])        # the pulse from which inputs to be changed onwards  
            data_switch = initialize_data(runtime, site)       # switching site
            if p < switch - 1:
                # stick to the starting site
                Ecosystem.reinitialization(data_initialization, data_initialization['Microbes_pp'], Output_init, mode, p, switch) 
            else:
                # switch to data of another site
                Ecosystem.reinitialization(data_switch, data_initialization['Microbes_pp'], Output_init, mode, p, switch)
    
    #...export the Output_init object to the output_folder using the export() funtion in the utility module 
    os.chdir('../'+output_folder)
    export(Output_init, site, outname)
    
main()