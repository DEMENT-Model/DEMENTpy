import pandas as pd
import numpy as np
from utility import LHS

class Enzyme():
    
    """
    This module deals with all traits related to exoenzymes and includes methods:
    1. enzyme_pool_initialization():
    2. enzyme_attributes():
    3. enzyme_Ea():
    4. enzyme_uptake_Ea():
    5. enzyme_Vmax():
    6. enzyme_uptake_Vmax():
    7. enzyme_Km():
    8. enzyme_uptake_Km():
    """
    
    
    def __init__(self,runtime,parameters,substrate_index):
        
        """
        Parameters:
            runtime:  runtime specifications by users
            parameters: inputted parameters
        """
        
        self.n_enzymes    = int(runtime.loc['n_enzymes',1])      
        self.n_substrates = int(runtime.loc['n_substrates',1])
        self.n_monomers   = self.n_substrates + 2
        self.n_uptake     = int(runtime.loc['n_uptake',1])
        self.substrate_index = substrate_index                   # index of substrates

        self.Enz_min    = parameters.loc['Enz_min',1]            # Initial min. enzyme present in terms of carbon
        self.Enz_max    = parameters.loc['Enz_max',1]            # Initial max. enzyme present in terms of carbon
        self.Enz_C_cost = parameters.loc['Enz_C_cost',1]         # Per enzyme C cost as a fraction of uptake:1
        self.Enz_N_cost = parameters.loc['Enz_N_cost',1]         # Per enzyme N cost as a fraction of C cost:0.3
        self.Enz_P_cost = parameters.loc['Enz_P_cost',1]         # Per enzyme P cost as a fraction of C cost:0
        self.Enz_Maint_cost = parameters.loc['Enz_Maint_cost',1] # Maintenence cost of enzyme production
        self.Uptake_Ea_min = parameters.loc['Uptake_Ea_min',1]       # Minimum activation energy for uptake
        self.Uptake_Ea_max = parameters.loc['Uptake_Ea_max',1]       # Maximum activation energy for uptake
        self.Vmax0_min     = parameters.loc['Vmax0_min',1]           # Minimum Vmax for enzyme
        self.Vmax0_max     = parameters.loc['Vmax0_max',1]           # Maximum Vmax for enzyme
        self.Uptake_Vmax0_min = parameters.loc['Uptake_Vmax0_min',1] # Minimum uptake Vmax
        self.Uptake_Vmax0_max = parameters.loc['Uptake_Vmax0_max',1] # Maximum uptake Vmax
        self.Specif_factor    = parameters.loc['Specif_factor',1]    # Efficiency-specificity
        self.Vmax_Km            = parameters.loc['Vmax_Km',1]               # slope for Km-Vmax relationship:1
        self.Vmax_Km_int        = parameters.loc['Vmax_Km_int',1]           # intercept for Km-Vmax relationship:0
        self.Km_min             = parameters.loc['Km_min',1]                # Minimum Km: default = 0.01
        self.Uptake_Vmax_Km     = parameters.loc['Uptake_Vmax_Km',1]        # slope for uptake Km-Vmax relationship: 0.02
        self.Uptake_Vmax_Km_int = parameters.loc['Uptake_Vmax_Km_int',1]    # intercept for Uptake Km-Vmax relationship:0
        self.Uptake_Km_min      = parameters.loc['Uptake_Km_min',1]         # Minimum uptake Km: 0.001
        self.Km_error           = parameters.loc['Km_error',1]              # Error term: default = 0
      
    def enzyme_pool_initialization(self):
        
        """
        Initialize the pool sizes of different enzymes
        -----------
        Parameters:
            Enz_min = 0
            Enz_max = 0
        Return:
            Enzymes_df: df
        """
        
        Enzymes_array = np.random.uniform(self.Enz_min,self.Enz_max,self.n_enzymes)
        index = ['Enz'+str(i) for i in range(1,self.n_enzymes+1)]
        Enzymes_df = pd.Series(data = Enzymes_array, index = index, name='C')

        return Enzymes_df
    
    
    
    def enzyme_attributes(self):
        
        """
        Enzyme stoichiometry and maintenence cost
        ----------
        Parameters:
            Enz_C_cost
            Enz_N_cost
            Enz_P_cost
            Enz_Maint_cost
        Return:
            EnzAttrib_df:
        """
        
        EnzAttrib_array = np.tile([self.Enz_C_cost,self.Enz_N_cost,self.Enz_P_cost,self.Enz_Maint_cost],(self.n_enzymes,1))
        index = ["Enz" + str(i) for i in range(1,self.n_enzymes + 1)]
        EnzAttrib_df = pd.DataFrame(data=EnzAttrib_array,index=index, columns=["C_cost","N_cost","P_cost","Maint_cost"])
        
        return EnzAttrib_df
    
    
    def enzyme_Ea(self,Ea_input):
        
        """ 
        Enzyme specificity matrix of activation energies
        -----------
        Input:
            user-supplied activation energy range for each substrate (for now min.= max.)
        Return:
            Ea_df.T: Rows:enzymes; cols: substrates
        """
        
        Ea_series = Ea_input.apply(lambda df: np.random.uniform(df['Ea_min'],df['Ea_max'],self.n_enzymes),axis=1)
        columns = ['Enz' + str(i) for i in range(1,self.n_enzymes + 1)]
        Ea_df = pd.DataFrame(data=Ea_series.tolist(), index=self.substrate_index,columns = columns) # note: .tolist() !!!!!!!!
        
        return Ea_df.T
    
    
    def enzyme_uptake_Ea(self):
        
        """
        Uptake activation energies constant across monomers
        -----------
        Parameters:
            Uptake_Ea_min: 35
            Uptake_Ea_max: 35
        Return:
            Rows are monomers; cols are uptake enzymes
        """
        
        Uptake_Ea_array = np.random.uniform(self.Uptake_Ea_min, self.Uptake_Ea_max, self.n_uptake*self.n_monomers)
        Uptake_Ea_array = Uptake_Ea_array.reshape(self.n_monomers,self.n_uptake)
        
        index = ['Mon'+str(i) for i in range(1,self.n_monomers+1)]
        columns = ['Upt'+str(i) for i in range(1,self.n_uptake+1)]
        Uptake_Ea_df = pd.DataFrame(data = Uptake_Ea_array,index = index,columns = columns)
        
        return Uptake_Ea_df
        
        
    def enzyme_Vmax(self,ReqEnz):
        
        """
        Pre-exponential constants for enzymes
        ------------
        Inputs:
          ReqEnz:     substrate required enzymes from the substrate module
          Vmax0_min:  5  (mg substrate mg-1 enzyme day-1); Minimum Vmax for enzyme
          Vmax0_max:  50 (mg substrate mg-1 enzyme day-1); Maximum Vmax for enzyme
          Specif_factor: Efficiency-specificity(1)
              
        Return:
          Vmax0_df:   Rows are substrates; cols are enzymes (1st used as input for the function below)
          Vmax0_df.T: enzyme*substrate; feeded to expand()
        """
        
        #Vmax0_array = np.random.uniform(self.Vmax0_min, self.Vmax0_max, self.n_substrates*self.n_enzymes)
        Vmax0_array = LHS(self.n_substrates*self.n_enzymes,self.Vmax0_min, self.Vmax0_max-self.Vmax0_min, 'uniform')
        Vmax0_array = Vmax0_array.reshape(self.n_substrates,self.n_enzymes)
        
        #index = ['Sub'+str(i) for i in range(1,self.n_substrates + 1)]
        columns = ['Enz'+str(i) for i in range(1,self.n_enzymes + 1)]
        Vmax0_df = pd.DataFrame(data=Vmax0_array, index=self.substrate_index, columns=columns)
        
        # Account for efficiency-specificity tradeoff by dividing Vmax_0 by the number of substrates (or monomers)
        # targeted and multiplied by a specificity factor
        RE1 = ReqEnz.loc['set1'].iloc[range(self.n_substrates),]
        RE2 = ReqEnz.loc['set2'].iloc[range(self.n_substrates),]
        RE2 = RE2.fillna(0)
        
        #...Total num. of substrates that each enzyme can target: series; by enzyme
        total_substrates = RE1.sum(axis=0) + RE2.sum(axis=0)
        if self.Specif_factor == 0:
            total_substrates[total_substrates>1] = 1
        else:
            total_substrates[total_substrates>1] = total_substrates[total_substrates>1]*self.Specif_factor
        
        Vmax0_df = Vmax0_df.divide(total_substrates,axis=1)
        Vmax0_df = Vmax0_df.fillna(0)
        Vmax0_df[np.isinf(Vmax0_df)] = 0
        
        
        return Vmax0_df, Vmax0_df.T
    
    
        
    def enzyme_uptake_Vmax(self,Uptake_ReqEnz):
        
        """
        Pre-exponential constants for uptake
        -----------------
        
        Input:
            Uptake_ReqEnz: uptake required enzymes; monomers * enzymes; from the monomer module
            Uptake_Vmax0_min: 1	(mg substrate mg-1 substrate day-1)	Minimum uptake Vmax
            Uptake_Vmax0_max: 10 (mg substrate mg-1 substrate day-1)	Maximum uptake Vmax
            Specif_factor: Efficiency-specificity (1)
        Return:
            Uptake_Vmax0_df: Rows are monomers; cols are uptake enzymes
        """
        
        #Uptake_Vmax0_array = np.random.uniform(self.Uptake_Vmax0_min,self.Uptake_Vmax0_max,self.n_uptake*self.n_monomers)
        Uptake_Vmax0_array = LHS(self.n_uptake*self.n_monomers,self.Uptake_Vmax0_min,self.Uptake_Vmax0_max-self.Uptake_Vmax0_min,'uniform')
        Uptake_Vmax0_array = Uptake_Vmax0_array.reshape(self.n_monomers,self.n_uptake)
        
        index = ['Mon'+str(i) for i in range(1,self.n_monomers+1)]
        columns = ['Upt'+str(i) for i in range(1,self.n_uptake+1)]
        Uptake_Vmax0_df = pd.DataFrame(data = Uptake_Vmax0_array,index = index,columns = columns)
        
        # implement the tradeoff with specificity
        total_monomers = Uptake_ReqEnz.sum(axis=0)
        if self.Specif_factor == 0:
            total_monomers[total_monomers>1] = 1
        else:
            total_monomers[total_monomers>1] = total_monomers[total_monomers>1] * self.Specif_factor
            
        Uptake_Vmax0_df = Uptake_Vmax0_df.divide(total_monomers,axis=1)
        Uptake_Vmax0_df = Uptake_Vmax0_df.fillna(0)
        Uptake_Vmax0_df[np.isinf(Uptake_Vmax0_df)] = 0
        
        return Uptake_Vmax0_df
    
    
            
    def enzyme_Km(self,Vmax0):
        
        """ 
        Implement Vmax-Km tradeoff to derive Km as a linear function of Vmax:
        ------------------------------------------
        Input:
          Vmax0:       from the above method
        Parameters:
          Vmax_Km:     slope for Km-Vmax relationship (1)
          Vmax_Km_int: intercept (0)
          Km_error:    error term (0);error term normally distributed with magnitude Km_error.)
          Km_min:      to which the minimum Km is constrained; 0.01
        Return:
          Km_df: substrate * enzyme
        """
        
        mean = Vmax0.mean().mean()*self.Km_error
        Km_array = abs(Vmax0.apply(lambda df:np.random.normal(df*self.Vmax_Km,mean))+self.Vmax_Km_int)
        
        #Km_array = abs(np.random.normal(Vmax0*self.Vmax_Km, Vmax0.mean().mean()*self.Km_error,Vmax0) + self.Vmax_Km_int)
        Km_array[Km_array < self.Km_min] = self.Km_min
        #index = ['Sub'+str(i) for i in range(1,self.n_substrates + 1)]
        columns = ['Enz'+str(i) for i in range(1,self.n_enzymes + 1)]
        Km_df = pd.DataFrame(data = Km_array, index=self.substrate_index, columns = columns)
        
        return Km_df
    
    
    def enzyme_uptake_Km(self,Uptake_Vmax0):
        
        """
        Implement Vmax-Km tradeoff to derive Uptake Km as a linear function of Uptake Vmax:
        ---------------
        Input:
            Uptake_Vmax0:       from the above method
        Parameters:
            Uptake_Vmax_Km:     slope of Km-Vmax correlation (0.2)
            Uptake_Vmax_Km_int: intercept (0)
            Km_error:           error term normally distributed with magnitude (0)
            Uptake_Km_min:      Minimum Km constrained to Km_min (0.001)
        Return:
            Uptake_Km_df
        """
        
        mean = Uptake_Vmax0.mean().mean()*self.Km_error
        Uptake_Km_array = abs(Uptake_Vmax0.apply(lambda df:np.random.normal(df*self.Uptake_Vmax_Km,mean))+self.Uptake_Vmax_Km_int)
        Uptake_Km_array[Uptake_Km_array < self.Uptake_Km_min] = self.Uptake_Km_min
        index = ['Mon'+str(i) for i in range(1,self.n_monomers + 1)]
        columns = ['Upt'+str(i) for i in range(1,self.n_uptake + 1)]
        Uptake_Km_df = pd.DataFrame(data = Uptake_Km_array,index = index,columns = columns)
        
        return Uptake_Km_df