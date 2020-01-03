# substrate.py module holding a class, Substrate()
# by Bin Wang on Dec. 26th, 2019ls

import pandas as pd
import numpy as np


class Substrate():
    """
    Deals with all properties related to substrates.
    
    Includs methods of:
        1) substrate_input(): input rates of varying substrates to the system
        2) substrate_produced_monomer():  substrate-produced monomers
        3) substrate_degradation_enzyme(): substrate-required enzymes
    """
    
    def __init__(self,runtime,parameters,substrates_init):
        """
        The constructor Substrate class.

        Parameters:
            runtime:    dataframe; user-specified parameters when running the model
            parameters: dataframe; parameters related to microbes, substrates,enzymes, and monomers
        """
        
        self.n_substrates      = int(runtime.loc['n_substrates',1])
        self.n_enzymes         = int(runtime.loc['n_enzymes',1])
        self.gridsize          = int(runtime.loc['gridsize',1])
        self.Enzymes_per_sub   = int(parameters.loc['Enzymes_per_sub',1])   # minimum # of enzymes degrading each substrate
        self.Avg_extra_req_enz = int(parameters.loc['Avg_extra_req_enz',1]) # average # of additional enzymes required for substrate degradation
        self.Substrates_start  = substrates_init.astype('float32')          # initial substrate concentrations
    
    
    def substrate_input(self,sub_mon_input):
        """
        Substrate inputs during simulation.

        Parameter:
            sub_mon_input: dataframe; user-provided substrates and monomers input data
        Return:
            SubInput_df: dataframe; datatype: float32
        """
        
        # Substrate input rates
        SubInputC = sub_mon_input['Sub'] # access only the Substrate column
        SubInputN = SubInputC * self.Substrates_start["N"]/self.Substrates_start["C"]
        SubInputP = SubInputC * self.Substrates_start["P"]/self.Substrates_start["C"]
        # Rename series name
        SubInputC.name = 'C' 
        SubInputN.name = 'N'
        SubInputP.name = 'P'
        SubInput_df = pd.concat([SubInputC,SubInputN,SubInputP],axis=1,sort=False)
        SubInput_df['DeadMic'] = SubInput_df['DeadEnz'] = 0  # Change NAs to 0
        SubInput_df = SubInput_df.astype('float32')

        return SubInput_df
        
  
    def substrate_produced_monomer(self):
        """
        Monomers produced by each substrate.

        Return:
            MonomersProduced_df: dataframe; row:substrate; col:monomers; all rows should sum to 1
        """
        
        MonomersProduced_array = np.concatenate((np.array([0]*self.n_substrates*2).reshape((self.n_substrates,2),order='F'),np.diagflat([1]*self.n_substrates)),axis=1)
        index   = ['Sub'+str(i) for i in range(1,self.n_substrates+1)]
        columns = ['Mon'+str(i) for i in range(-1,self.n_substrates+1)]
        MonomersProduced_df = pd.DataFrame(data=MonomersProduced_array,index=index,columns=columns,dtype='int8')
        MonomersProduced_df.rename(columns = {'Mon-1':"NH4",'Mon0':"PO4",'Mon1':"DeadMic",'Mon2':"DeadEnz"},inplace=True)

        return MonomersProduced_df
        
    
    def substrate_degradation_enzyme(self):
        """
        Derive the required enzymes of each substrate.

        Return:
            ReqEnz_df: 3-D DataFrame: sets (2) * substrates * enzymes
        .............................................................
        Strcuture illustration:
        .............................................................
                   Enz1 Enz2 ... Enzn
        Set 1 Sub1 ------------------
        Set 1 Sub1 ------------------
        .     .    .     .       .
        .     .    .     .       .
        Set 1 Subn ------------------
        Set 2 Sub1 ------------------
        Set 2 Sub1 ------------------
        .     .    .     .       .
        .     .    .     .       .
        Set 2 Subn ------------------
        .............................................................
        """

        # Every enzyme has a substrate
        probability_list_enz = [0] * self.n_enzymes
        probability_list_enz[0:self.Enzymes_per_sub] = [1]*self.Enzymes_per_sub
        ReqEnz1_list  = [np.random.choice(probability_list_enz,self.n_enzymes,replace=False) for i in range(self.n_substrates)]
        ReqEnz1_array = np.vstack(ReqEnz1_list)
        # Ensure every substrate has an enzyme
        probability_list_sub = [0] * self.n_substrates
        probability_list_sub[0] = 1
        for i in range(self.n_enzymes):
            if sum(ReqEnz1_array[:,i])==0:
                ReqEnz1_array[:,i] = np.random.choice(probability_list_sub,self.n_substrates,replace=False)
        
        # Choose some substrates that require multiple enzymes        
        probability_list = [0] * self.n_enzymes
        probability_list[0:self.Avg_extra_req_enz] = [1]*self.Avg_extra_req_enz
        ReqEnz2_array = np.random.choice(probability_list,self.n_substrates * self.n_enzymes,replace=True)        
        ReqEnz2_array = ReqEnz2_array.reshape(self.n_substrates,self.n_enzymes)
        ReqEnz2_array[ReqEnz1_array==1] = 0
        for i in range(self.n_substrates):
            if sum(ReqEnz2_array[i,]) == 0:
                ReqEnz2_array[i,] = 0  # no extra enzymes assigned
        
        # achieved via MultiIndex
        ReqEnz_array = np.concatenate((np.tile(ReqEnz1_array,(self.gridsize,1)),np.tile(ReqEnz2_array,(self.gridsize,1))),axis=0)
        index = [np.array(['set1']*self.n_substrates*self.gridsize + ['set2']*self.n_substrates*self.gridsize),
                 np.array(["Sub" + str(i) for i in range(1,self.n_substrates + 1)]*self.gridsize*2)]
        columns = ['Enz' + str(i) for i in range(1,self.n_enzymes + 1)]
        ReqEnz_df = pd.DataFrame(data=ReqEnz_array,index=index,columns=columns,dtype='int8')

        return ReqEnz_df
