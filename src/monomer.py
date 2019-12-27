"""
This module, monomer.py, handles monomer-related properties in a class, Monomer().
"""


import pandas as pd
import numpy as np


class Monomer():
    
    """
    This class deals with all calculations closely related to monomers. Methods include:
    1.monomer_initialization():
    2.monomer_input_rate():   
    3.monomer_uptake_reqenzyme():calculate the required enzymes for taking up monomers in the current system.
    """
    
    
    def __init__(self,runtime,parameters):
        
        
        self.n_monomers  = int(runtime.loc['n_substrates',1]) + 2                  # Number of monomers within the system
        self.n_uptake    = int(runtime.loc['n_uptake',1])                          # Number of uptake transporters for each taxon
        
        self.Monomer_Substrate_Ratio = parameters.loc['Monomer_Substrate_Ratio',1] # default:0; amount of initial monomer relative to substrate per grid box
        self.Uptake_per_monomer      = int(parameters.loc['Uptake_per_monomer',1]) # default:1; number of transporters per monomer 
        self.Init_NH4  = parameters.loc['Init_NH4',1]                              # Initial NH4
        self.Init_PO4  = parameters.loc['Init_PO4',1]                              # Initial PO4
        self.Input_NH4 = parameters.loc['Input_NH4',1]                             # Input of NH4
        self.Input_PO4 = parameters.loc['Input_PO4',1]                             # Input of PO4
        
        
    def monomer_initialization(self,substrates_init):
        
        """
        #...derive the initial Monomers: dataframe; used for later calculations w/ two purposes:
        #...1) get sums over the grid, which is further used to derive the average monomer distri. over the grid
        #...2) provide the data structure for MonomerRatios that stores monomer stoichiometry info.
        ------------------------------------------
        Parameters:
            substrates_init: dataframe; initial substrates pool; from the substrate module
        Return:
            Monomers_df: dataframe (shape:14*3).

        """
        
        # Monomer pool sizes for all elements
        Monomers_array = np.concatenate((np.stack([[0,self.Init_NH4,0],[0,0,self.Init_PO4]]),substrates_init.values*self.Monomer_Substrate_Ratio),axis=0)    
        index = ["NH4","PO4","DeadMic","DeadEnz"] + ["Mon" + str(i) for i in range(3,self.n_monomers-2 + 1)]  #NOTE: index starts from 3
        Monomers_df = pd.DataFrame(data=Monomers_array, index=index, columns=["C","N","P"],dtype='float32')
        
        return Monomers_df
    
    
    def monomer_ratios(self,Monomers_df):
        
        """
        initialization of 'monomer_ratio
        -----------------------------------------------------
        Parameter:
            Monomers_df:
        Return:
            Monomer_ratios
        '"""
        
        is_NH4 = Monomers_df.index == "NH4"
        is_PO4 = Monomers_df.index == "PO4"
        Monomer_ratios = Monomers_df.copy(deep=True)
        Monomer_ratios[:] = 0.0
        Monomer_ratios.loc[is_NH4,'N'] = 1.0
        Monomer_ratios.loc[is_PO4,'P'] = 1.0
        
        return Monomer_ratios
        
        
    def monomer_input_rate(self,sub_mon_input):
        
        """
        Derive the MonInput: monomer input rates; series
        -----------------------------------------------------
        Parameters:
          sub_mon_input:
        Return:
          MonInput:
        """
        
        # load the inputs without NH4 and PO4
        monomer_input = sub_mon_input['Mon']
        # supplement with NH4 and PO4 from the parameters
        MonInput = pd.concat([pd.Series([self.Input_NH4,self.Input_PO4],index=['Input_NH4','Input_PO4']),monomer_input],sort=False)
        MonInput = MonInput.astype('float32')

        return MonInput
       
         
    def monomer_uptake_reqenzyme(self):
        
        """
        derive the Uptake_ReqEnz: dataframe; Rows-monomers;cols-uptake enzymes
        Same number within a row implies redundancy
        Make sure each monomer is taken up by at least one transporter and every transporter takes up at least one monomer
        """
        
        probability_list = [0] * self.n_uptake
        probability_list[0:self.Uptake_per_monomer] = [1]*self.Uptake_per_monomer
        Uptake_ReqEnz_list = [np.random.choice(probability_list,self.n_uptake,replace=False) for i in range(self.n_monomers)]

        index   = ["Mon" + str(i) for i in range(1,self.n_monomers+1)]
        columns = ['Upt' + str(i) for i in range(1,self.n_uptake+1)]
        Uptake_ReqEnz_df = pd.DataFrame(data=np.array(Uptake_ReqEnz_list).reshape(self.n_monomers,self.n_uptake),index=index,columns=columns,dtype='int8')
        # ensure every monomer has a transporter
        probability_list    = [0]* self.n_monomers
        probability_list[0] = 1
        for i in range(self.n_uptake):
            if sum(Uptake_ReqEnz_df.iloc[:,i]) == 0:
                Uptake_ReqEnz_df.iloc[:,i] = np.random.choice(probability_list,self.n_monomers,replace=False).astype('int8')
                
        return Uptake_ReqEnz_df