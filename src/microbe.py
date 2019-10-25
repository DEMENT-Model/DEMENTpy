"""
This microbe module has one class and two functions:
    
    Microbe(): class
    microbe_osmo_psi(): function
    microbe_mortality_psi(): function
        
-------------------------------------------------------------------------------
Authored by Bin Wang (wbwenwu@gmail.com/bwang7@uci.edu)
Last modified: Sep. 5th, 2019
"""


import numpy as np
import pandas as pd
from utility import LHS

class Microbe():
    
    """
    This class holds all variables related to microbes and methods involving
    composition,stoichiometry, enzyme and gene production, as well as responses
    to environmental factors. Methods include:
        1) microbial_community_initialization(): initialize microbial community on the spatial grid
        2) minimum_cell_quota():          get minimum ratios
        3) microbe_enzyme_gene():         derive taxon-specific genes for enzyme
        4) microbe_osmolyte_gene():       derive genes encoding for osmolytes
        5) microbe_uptake_gene():         derive transporter genes
        6) microbe_uptake_cost():         metabolic cost of producing transporter
        7) microbe_enzproduction_rate():  cost of producing enzymes
        8) microbe_osmoproduction_rate(): cost of producing osmolytes
    
    """
    
    
    def __init__(self,runtime,parameters):
        
        """
        Parameters:
            runtime:     user-specified runtime parameters
            paramteters: input dataframe of all parameters
        """
        
        # system setup parameters
        self.gridsize    = int(runtime.loc['gridsize',1])
        self.n_taxa      = int(runtime.loc['n_taxa',1])
        self.n_enzymes   = int(runtime.loc['n_enzymes',1])
        self.n_substrates= int(runtime.loc['n_substrates',1])
        self.n_monomers  = int(runtime.loc['n_substrates',1])+2  # why+2? b/c two inorganic monomers
        self.n_uptake    = int(runtime.loc['n_uptake',1])        # Number of uptake transporters for each taxon
        self.taxa_per_box= runtime.loc['taxa_per_box',1]         # Probability of each taxon entering a grid cell
        fb = runtime.loc['fb',1]                                 # Probability of fungal taxa 
        self.fb = np.random.choice([1,0], self.n_taxa, replace=True, p=[fb,(1-fb)]) #Index of fungal taxa
        
        # microbial traits
        self.Cfrac_b    = parameters.loc['Cfrac_b',1]     # Bacterial C fraction: 0.825 mg mg-1
        self.Nfrac_b    = parameters.loc['Nfrac_b',1]     # Bacterial N fraction: 0.16  mg mg-1
        self.Pfrac_b    = parameters.loc['Pfrac_b',1]     # Bacterial P fraction: 0.015 mg mg-1
        self.Cfrac_f    = parameters.loc['Cfrac_f',1]     # Fungal C fraction: 0.9	mg mg-1
        self.Nfrac_f    = parameters.loc['Nfrac_f',1]     # Fungal N fraction: 0.09	mg mg-1
        self.Pfrac_f    = parameters.loc['Pfrac_f',1]     # Fungal P fraction: 0.01	mg mg-1
        self.max_size_b = parameters.loc['max_size_b',1]  # C quota threshold for bacterial cell division: 2	mg cm-3
        self.max_size_f = parameters.loc['max_size_f',1]  # C quota threshold for fungal cell division: 50 mg cm-3
        self.Crange     = parameters.loc['Crange',1]      # Tolerance on C fraction: 0.09 mg mg-1
        self.Nrange     = parameters.loc['Nrange',1]      # Tolerance on N fraction: 0.04 mg mg-1
        self.Prange     = parameters.loc['Prange',1]      # Tolerance on P fraction: 0.005 mg mg-1
        
        # metabolism-related traits
        self.Uptake_C_cost_min = parameters.loc['Uptake_C_cost_min',1]  # Minimum C cost of one transporter gene being encoded 
        self.Uptake_C_cost_max = parameters.loc['Uptake_C_cost_max',1]  # Maximum C cost of one transporter gene being encoded
        self.NormalizeUptake   = parameters.loc['NormalizeUptake',1]    # Normalize uptake investment for the number of uptake genes;default:0
        
        self.Enz_per_taxon_min = int(parameters.loc['Enz_per_taxon_min',1]) # Minimum number of enzyme genes a taxon can hold
        self.Enz_per_taxon_max = int(parameters.loc['Enz_per_taxon_max',1]) # Maximum number of enzyme genes a taxon can hold
        self.Constit_Prod_min  = parameters.loc['Constit_Prod_min',1] # =0.00001;Minimum per enzyme production cost as a fraction of biomass
        self.Constit_Prod_max  = parameters.loc['Constit_Prod_max',1] # =0.0001; Maximum per enzyme production cost as a fraction of biomass
        self.Enz_Prod_min      = parameters.loc['Enz_Prod_min',1]     # =0.00001;Minimum per enzyme production cost as a fraction of C uptake rate
        self.Enz_Prod_max      = parameters.loc['Enz_Prod_max',1]     # =0.0001; Maximum ...
        self.NormalizeProd     = parameters.loc['NormalizeProd',1]    # Normalize enzyme production for the number of enzyme genes;default:0
        
        
        self.n_osmolyte = 10            # system-specified number of osmotic compound
        
        self.Osmo_per_taxon_min = 100  # Minimum number of osmotic gene 
        self.Osmo_per_taxon_max = 101  # Max. of osmotic gene
        
        self.Osmo_Consti_Prod_min = 0.00001 # constitutive cost min
        self.Osmo_Consti_Prod_max = 0.0001  # constitutive cost max
        
        self.Osmo_Induci_Prod_min = 0.0001  # inducible cost min
        self.Osmo_Induci_Prod_max = 0.001   # inducible cost max
        
        self.death_rate_bac = 0.001   # Bacterial basal mortality prob.
        self.death_rate_fun = 0.0002  # Fungal basal mortality probability
        self.beta_bac = 10            # mortality coefficient
        self.beta_fun = 10            # mortality coefficient
        
    
    
    def microbial_community_initialization(self):
    
        """
        Initialize a microbial community on the grid with bacteria and/or fungi:
        -> Firstly create a dataframe with only bacteria
        -> Substitute part of bacteria with fungi to create a community comprised of both bacteria and fungi
        -> Randomly place this communtiy on the sptatial grid
        -> Perform stats to the initialized community and output them for record
        
        Parameters:
            Cfrac_b: Bacterial C fraction:    0.825 mg mg-1
            Nfrac_b: Bacterial N fraction:    0.160 mg mg-1
            Pfrac_b: Bacterial P fraction:    0.015 mg mg-1
            Cfrac_f: Fungal C fraction:       0.90	mg mg-1
            Nfrac_f: Fungal N fraction:       0.09	mg mg-1
            Pfrac_f: Fungal P fraction:       0.01	mg mg-1
            max_size_b: C quota threshold for bacterial cell division
            max_size_f: C quota threshold for fungal cell division
        
        Returns:
            -1. microbes_pp: microbial community preceding placement on the spatial grid
            -2. microbes_df: microbial community on the spatial grid
            -3. fb_grid:     index of fungi taxa in the microbial community
            -4. Bac_density: density of bacterial cells
            -5. Fun_density: density of fungal cells
            
        Notes:
        
        References:
            
        """
        
        BacC = 0.5  * self.max_size_b
        BacN = BacC * self.Nfrac_b/self.Cfrac_b
        BacP = BacC * self.Pfrac_b/self.Cfrac_b
        microbes_array = np.tile([BacC,BacN,BacP],(self.n_taxa*self.gridsize,1))
        index = ["Tax" + str(i) for i in range(1,self.n_taxa + 1)] * self.gridsize
        microbes_df = pd.DataFrame(data = microbes_array,
                                   index = index,
                                   columns = ["C","N","P"])
        
        # Obtain the Fungi index by expanding the fb to the spatial grid
        fb_grid = np.tile(self.fb,self.gridsize)
        
        # Fungal pool sizes for all elements
        FunC = 0.5 * self.max_size_f
        FunN = FunC* self.Nfrac_f/self.Cfrac_f
        FunP = FunC* self.Pfrac_f/self.Cfrac_f
        
        #...Substitute with fungal pools of elements: fb_grid == 1
        microbes_df.loc[fb_grid == 1, "C"] = FunC
        microbes_df.loc[fb_grid == 1, "N"] = FunN
        microbes_df.loc[fb_grid == 1, "P"] = FunP
        
        # Export the microbial mass before placement on the grid
        microbes_pp = microbes_df.copy()
        
        #...Derive the number of fungi and bacteria taxa
        Bac_index = microbes_df['C'] == BacC
        Fun_index = microbes_df['C'] == FunC
        
        Bac_taxa = int(sum(Bac_index)/self.gridsize)
        Fun_taxa = int(sum(Fun_index)/self.gridsize)
        print('Before placement--','Bac_taxa=',Bac_taxa,'Fun_taxa=',Fun_taxa)
        
        
        # Randomly place the microbial community created above on the spatial grid to
        #...initialize an spatially explicit microbial community by
        #...creating an array choose_taxa [1:yes;0:no]
        pb = self.taxa_per_box
        choose_taxa = np.random.choice([1,0], self.n_taxa*self.gridsize,replace=True, p=[pb,(1-pb)])
        pf = pb * self.max_size_b/self.max_size_f
        choose_taxa[fb_grid==1] = np.random.choice([1,0], sum(fb_grid),replace=True, p=[pf,(1-pf)])
        microbes_df.loc[choose_taxa==0] = 0
        
        
        #...Derive the number of fungi and bacteria taxa
        Bac_index = microbes_df['C'] == BacC
        Fun_index = microbes_df['C'] == FunC
        
        Bac_taxa = microbes_df[Bac_index].groupby(level=0).sum().shape[0]
        Fun_taxa = microbes_df[Fun_index].groupby(level=0).sum().shape[0]
        
        print('After placement--','Bac_taxa=',Bac_taxa,'Fun_taxa=',Fun_taxa)
        
        Bac_density = sum(Bac_index)*BacC/self.gridsize
        Fun_density = sum(Fun_index)*FunC/self.gridsize
        
        print('After placement--','Bac_density=',Bac_density,'Fun_density=',Fun_density)
        
        return microbes_pp,microbes_df,fb_grid,Bac_density,Fun_density
    
    
    
    def minimum_cell_quota(self):
        
        """
        This will be used in the mortality calculation
        -----------------------
        
        Parameters:
            Cfrac_b: Bacterial C fraction:    0.825 mg mg-1
            Nfrac_b: Bacterial N fraction:    0.160 mg mg-1
            Pfrac_b: Bacterial P fraction:    0.015 mg mg-1
            Cfrac_f: Fungal C fraction:       0.90  mg mg-1
            Nfrac_f: Fungal N fraction:       0.09  mg mg-1
            Pfrac_f: Fungal P fraction:       0.01  mg mg-1
            Crange:  Tolerance on C fraction: 0.090 mg mg-1
            Nrange:  Tolerance on N fraction: 0.040 mg mg-1
            Prange:  Tolerance on P fraction: 0.005 mg mg-1
        
        Output:
            MinRatios: dataframe
            
        """
        
        # First derive the optimal stoichiometry of bacterial taxa
        OptimalRatios_list = [self.Cfrac_b,self.Nfrac_b,self.Pfrac_b] * self.n_taxa
        OptimalRatios_array = np.array(OptimalRatios_list).reshape(self.n_taxa,3)
        index = ["Tax" + str(i) for i in range(1,self.n_taxa + 1)]
        OptimalRatios_df = pd.DataFrame(data = OptimalRatios_array,
                                        index = index,
                                        columns = ["C","N","P"])
        # Then substitute with fungal stoichiometry
        OptimalRatios_df[self.fb==1] = np.tile([self.Cfrac_f,self.Nfrac_f,self.Pfrac_f],(sum(self.fb),1)) 
        
        # Calcualte ratio range
        RangeRatios_list  = [self.Crange,self.Nrange,self.Prange] * self.n_taxa
        RangeRatios_array = np.array(RangeRatios_list).reshape(self.n_taxa,3)
        RangeRatios_df    = pd.DataFrame(data = RangeRatios_array,
                                         index = index,
                                         columns = ["C","N","P"])
  
        # Calculate minimum cell quotas regarding C, N, & P
        MinRatios = OptimalRatios_df - RangeRatios_df
        
        return MinRatios
        
        
   
    def microbe_enzyme_gene(self):
        
        """ 
        Not all taxa need to have an enzyme gene; this method draws the
        number of genes per taxon from a uniform distribution:
        ----------------------------------    
        Parameters:
            Enz_per_taxon_min: minimum # genes a taxon can have (0)
            Enz_per_taxon_max: maximum # genes a taxon can have (40)
            
        Return:
            'EnzGenes_df': Rows are taxa; cols are genes;values: 0/1
        """
        
        genes_per_taxon = np.random.choice(range(self.Enz_per_taxon_min,self.Enz_per_taxon_max+1),self.n_taxa,replace=True)
        
        # Number of genes needed to produce the number of enzymes that system requires 
        n_genes = self.n_enzymes
        
        #EnzGenes_list = [None]*self.n_taxa
        def trial(i):
            probability_list = [0]*n_genes
            probability_list[0:int(genes_per_taxon[i])] = [1]*int(genes_per_taxon[i])
            taxon = np.random.choice(probability_list,n_genes,replace=False)
            return taxon
        
        EnzGenes_list = [trial(i) for i in range(self.n_taxa)]
        
        index = ["Tax" + str(i) for i in range(1,self.n_taxa+1)]
        columns = ['Enz' + str(i) for i in range(1,n_genes+1)]
        EnzGenes_df = pd.DataFrame(data = np.vstack(EnzGenes_list).reshape(self.n_taxa,n_genes),
                                   index = index, 
                                   columns = columns)
        
        return EnzGenes_df
    
    
    def microbe_osmolyte_gene(self):
        
        """
        Derive the taxon-specific number of genes encoding for metabolic production of osmolytes
        this method draws the number of genes per taxon from a uniform distribution.
        ----------------------------------
        Parameters:
            n_osmolyte:
            Osmo_per_taxon_min:
            Osmo_per_taxon_max:
                
        Return:
            OsmoGenes_df: row: taxon; col: genes; values: 0/1
        """
        genes_per_taxon = np.random.choice(range(self.Osmo_per_taxon_min,self.Osmo_per_taxon_max+1),self.n_taxa,replace=True)
        
        # Number of genes needed to produce the system-required number of osmolytes
        n_genes = self.n_osmolyte
        
        # EnzGenes_list = [None]*self.n_taxa
        def trial(i):
            probability_list = [0]*n_genes
            probability_list[0:int(genes_per_taxon[i])] = [1]*int(genes_per_taxon[i])
            taxon = np.random.choice(probability_list,n_genes,replace=False)
            return taxon
        
        OsmoGenes_list = [trial(i) for i in range(self.n_taxa)]
        
        index = ["Tax" + str(i) for i in range(1,self.n_taxa+1)]
        columns = ['Osmo' + str(i) for i in range(1,n_genes+1)]
        OsmoGenes_df = pd.DataFrame(data = np.vstack(OsmoGenes_list).reshape(self.n_taxa,n_genes),
                                    index = index, 
                                    columns = columns)
        
        return OsmoGenes_df
        
        
     
    def microbe_uptake_gene(self,ReqEnz,EnzGenes,MonomersProduced):
        
        """
        Derive the genes encoding uptake enzymes(transporters); all taxa must have at least one uptake gene
        ----------------------
        Inputs:
            ReqEnz:           substrate-required enzymes;3-D dataframe (2 * n_substrates * n_enzymes); from the Substrate class
            EnzGenes:         genes encoding enzymes degrading substrates; 2-D df(taxon*gene); from the 'microbe_enzyme_gene' method
            MonomersProduced: sub * monomers; from the Substrae class
        
        Output:
            UptakeGenes: Rows-taxa; Cols-genes; Values: 0 (no gene) or 1 (yes)
        """
        
        RE2 = ReqEnz.loc['set2'].iloc[0:self.n_substrates,]
        Sub_enz = RE2 + ReqEnz.loc['set1'].iloc[0:self.n_substrates,]
        
        # Matrix multiplication to relate taxa to the monomers they can generate with their enzymes
        UptakeGenes = EnzGenes@np.transpose(Sub_enz)@MonomersProduced
        
        UptakeGenes.iloc[:,0:2] = 1
        UptakeGenes[lambda df: df > 0] = 1   #... = UptakeGenes[UptakeGenes > 0] = 1
        
        # Ensure every taxon is likely to have an uptake enzyme for at least 1 organic monomer
        # Not guaranteed unless uptake_prob = 1 (refer to the R version)
        probability_list = [0]*(self.n_monomers-2)
        probability_list[0] = 1
        for i in range(self.n_taxa):
            if sum(UptakeGenes.iloc[i,2:self.n_monomers]) == 0:
                UptakeGenes.iloc[i,2:self.n_monomers] = np.random.choice(probability_list,self.n_monomers-2,replace=False)

        # Give each taxon a random number of additional uptake genes between the number they have and n_upgenes
        for i in range(self.n_taxa):
            n_zero = len(UptakeGenes.iloc[i][UptakeGenes.iloc[i]==0])
            if n_zero == 0: # has all genes
                continue
            probability_list = [0]*n_zero
            locator = np.random.choice(range(n_zero+1),1) # 0 <= locator <= n_zero
            probability_list[0:int(locator)] = [1]*int(locator)
            UptakeGenes.iloc[i][UptakeGenes.iloc[i] == 0] = np.random.choice(probability_list,n_zero,replace=False)
                
        return UptakeGenes
    
     

    def microbe_uptake_cost(self,UptakeGenesForEnz):
                                 
        """
        Derive the taxon-specific cost of every single gene of uptake transporter (i.e.,fraction of biomass as transporters).
        Note this cost (in terms of Fraction of Biomass C) is same among different genes from the same taxon
        
        Input:
            UptakeGenesForEnz: taxon-specific transporter genes; from the above method
            
        Parameters:
            NormalizeUptake: Normalize uptake investment for the number of uptake genes (0: No, 1: Yes)
            Uptake_C_cost_min: 0.01	transporter mg-1 biomass C
            Uptake_C_cost_max: 0.1	transporter mg-1 biomass C
              
        Outputs:
            UptakeProd_series:
            UptakeGenes_Cost:
        """                         
        
        # If true (== 1), then the uptake potential is normalized to the number of uptake genes (n_uptake)
        if self.NormalizeUptake == 1:
            Normalize = UptakeGenesForEnz.sum(axis=1)/self.n_uptake
        else:
            Normalize = 1
        UptakeGenes = UptakeGenesForEnz.divide(Normalize,axis=0)
        
        # Choose total amount of uptake allocation at random
        # Old method:
        # UptakeProd = np.random.uniform(self.Uptake_C_cost_min,self.Uptake_C_cost_max,self.n_taxa)
        # Use LHS sampling instead:
        index = ["Tax" + str(i) for i in range(1,self.n_taxa+1)]
        UptakeProd_array = LHS(self.n_taxa,self.Uptake_C_cost_min,self.Uptake_C_cost_max-self.Uptake_C_cost_min,'uniform')
        UptakeProd_series= pd.Series(data=UptakeProd_array,index=index)
        
        UptakeGenes_Cost = UptakeGenes.mul(UptakeProd_series,axis=0)  
        
        return UptakeProd_series, UptakeGenes_Cost
    
    

    def microbe_enzproduction_rate(self,EnzGenes,EnzAttrib):
        
        """
        Derive the taxon-specific fraction of 'available C' as enzymes: note that
        this fraction only varies with taxon; every gene within a taxon is same.
        
        Inputs:
            EnzGenes: taxon-specific available genes for enzyme;dataframe (taxon * genes); from the above microbe_enzyme_gene() method
            EnzAttrib: basic enzyme stoichiometry;dataframe [enzyme * 4 (C,N,P,maintenence cost)]
            
        Parameters:
            Constit_Prod_min:
            Constit_Prod_max:
            Enz_Prod_min:
            Enz_Prod_max:
            NormalizeProd: Normalize enzyme production for the number of enzyme genes (0 for no, 1 for yes)
            
        Outputs:
            Tax_Induce_Enzyme_C:  taxon-specific fraction of available C as enzyme for each enzyme
            Tax_Constit_Enzyme_C: taxon-specific fraction of available C as enzyme for each enzyme
        """
        
        
        #...method 1....
        #...Tax_EnzProdInduce: series;storing the production rate (inducible or constituitive) of each taxon (axis=1)
        #index = ["Tax" + str(i) for i in range(1,self.n_taxa+1)]
        #Tax_EnzProdConstit_array = np.random.uniform(self.Constit_Prod_min,self.Constit_Prod_max,self.n_taxa)
        #Tax_EnzProdInduce_array  = np.random.uniform(self.Enz_Prod_min,self.Enz_Prod_max,self.n_taxa)
        #Tax_EnzProdConstit = pd.Series(data=Tax_EnzProdConstit_array,index=index)
        #Tax_EnzProdInduce  = pd.Series(data=Tax_EnzProdInduce_array,index=index)
        
        #...LHS sampling...
        index = ["Tax" + str(i) for i in range(1,self.n_taxa+1)]
        Tax_EnzProdConstit_array = LHS(self.n_taxa,self.Constit_Prod_min,self.Constit_Prod_max-self.Constit_Prod_min,'uniform')
        Tax_EnzProdInduce_array  = LHS(self.n_taxa,self.Enz_Prod_min,self.Enz_Prod_max-self.Enz_Prod_min,'uniform')
        Tax_EnzProdConstit = pd.Series(data=Tax_EnzProdConstit_array,index=index)
        Tax_EnzProdInduce  = pd.Series(data=Tax_EnzProdInduce_array,index=index)
        
        #...method 2....
        #...Tax_EnzProdInduce: series;storing the production rate (inducible or constituitive) of each taxon (axis=1)  
        #Tax_EnzProdConstit= EnzGenes.apply(lambda df: np.asscalar(np.random.uniform(self.Constit_Prod_min,self.Constit_Prod_max,1)),axis=1)
        #Tax_EnzProdInduce = EnzGenes.apply(lambda df: np.asscalar(np.random.uniform(self.Enz_Prod_min,self.Enz_Prod_max,1)),axis=1)
        
        #...derive the production rate of every single gene of each taxon
        EnzProdInduce = EnzGenes.mul(Tax_EnzProdInduce,axis=0)
        EnzProdConstit= EnzGenes.mul(Tax_EnzProdConstit,axis=0)
        #...Normalize enzyme production of Inductive and constituitive enzyme production
        if self.NormalizeProd == 1:
            Normalize = EnzGenes.sum(axis=1)/self.Enz_per_taxon_max
        else:
            Normalize = 1
        EnzProdInduce  = EnzProdInduce.divide(Normalize,axis=0)
        EnzProdConstit = EnzProdConstit.divide(Normalize,axis=0)
        #EnzProdInduce.fillna(0)  
        #EnzProdConstit.fillna(0)
                                        
        #...Relative enzyme carbon cost for each enzyme gene of each taxon
        Tax_Induce_Enzyme_C  = EnzProdInduce.mul(EnzAttrib["C_cost"], axis=1)  #C_cost = 1; so it doesn't matter
        Tax_Constit_Enzyme_C = EnzProdConstit.mul(EnzAttrib["C_cost"],axis=1)
    
        return Tax_EnzProdConstit, Tax_EnzProdInduce, Tax_Constit_Enzyme_C, Tax_Induce_Enzyme_C 
    
    
    def microbe_osmoproduction_rate(self,OsmoGenes_df):
        
        """
        Distribution of osmolyte production rate (i.e.,proportion of available C as osmolytes) 
        
        Parameters:
            OsmoGenes_df: taxon-specific available genes for osmolyte;
                          dataframe (taxon * genes); Generated by the above microbe_osmolyte_gene() method
            
        Parameters:
            Osmo_Constit_Prod_min:
            Osmo_Constit_Prod_max:
            Osmo_Induci_Prod_min:
            Osmo_Induci_Prod_max:
            NormalizeProd:  0
            
        Output:
            Tax_OsmoProd_Consti_series:
            Tax_OsmoProd_Induci_series:
            Tax_Consti_Osmo_C: taxon-specific fraction of available C as osmolytes for each gene
            Tax_Induce_Osmo_C: taxon-specific fraction of available C as osmolytes for each gene
        """
        
        #n_osmo = int(1)
        #Tax_OsmoProd_Constit = np.random.uniform(self.Osmo_Constit_Prod_min,self.Osmo_Constit_Prod_max,self.n_taxa)
        #Tax_OsmoProd_Induci  = np.random.uniform(self.Osmo_Induci_Prod_min,self.Osmo_Induci_Prod_max,self.n_taxa)
        
        Tax_OsmoProd_Consti = LHS(self.n_taxa,self.Osmo_Consti_Prod_min,self.Osmo_Consti_Prod_max-self.Osmo_Consti_Prod_min,'uniform')
        Tax_OsmoProd_Induci = LHS(self.n_taxa,self.Osmo_Induci_Prod_min,self.Osmo_Induci_Prod_max-self.Osmo_Induci_Prod_min,'uniform')
        
        
        index = ["Tax" + str(i) for i in range(1,self.n_taxa+1)]
        #columns = ['Osmo'+ str(i) for i in range(1,n_osmo+1)]
        Tax_OsmoProd_Consti_series = pd.Series(data=Tax_OsmoProd_Consti,index = index)
        Tax_OsmoProd_Induci_series = pd.Series(data=Tax_OsmoProd_Induci,index = index)
                                              
        #...derive the production rate of every single gene of each taxon
        Tax_Consti_Osmo_C = OsmoGenes_df.mul(Tax_OsmoProd_Consti_series,axis=0)
        Tax_Induci_Osmo_C = OsmoGenes_df.mul(Tax_OsmoProd_Induci_series,axis=0)
        
        
        return Tax_OsmoProd_Consti_series, Tax_OsmoProd_Induci_series, Tax_Consti_Osmo_C, Tax_Induci_Osmo_C
    
    
    
    def microbe_drought_tol(self,Tax_Induci_Osmo_C):
        
        """
        Drought tolerance is postively correlated with taxon-specific inducible
        osmotic allocation efficiency.
        --------------------------
        input:
            Tax_Induci_Osmo_C:
        
        """
        Tax_Osmo_Alloc = Tax_Induci_Osmo_C.sum(axis=1) #+ Tax_Consti_Osmo_C.sum(axis=1)
        #Tax_tolerance = Tax_Osmo_Alloc.rank(axis=0,method='min')/self.n_taxa
        Tax_tolerance = (Tax_Osmo_Alloc - Tax_Osmo_Alloc.min())/(Tax_Osmo_Alloc.max()-Tax_Osmo_Alloc.min())
        
        Tax_tolerance = Tax_tolerance.fillna(0)
        Tax_tolerance[np.isinf(Tax_tolerance)] = 0
        
        return Tax_tolerance
    
    
    
    def microbe_mortality(self,fb_grid):
        
        """
        Derive taxon-specific microbial mortality parameters
        
        Parameters:
            fb_grid: index of bacteria (0) vs fungi (1) over the grid
            
        Returns:
            death_rate: spatal taxon-specific basal death probability;array
            beta:       spatial beta;array
            
        """
        death_rate = (1-fb_grid) * self.death_rate_bac + fb_grid * self.death_rate_fun
        beta = (1-fb_grid) * self.beta_bac + fb_grid * self.beta_fun
        
        
        return beta, death_rate


def microbe_osmo_psi(wp,alfa,wp_fc,wp_th):
    
    """
    Inducible osmolytes' production triggered when PSI declines to a **threshold** value,wp_fc,
    below which the production increaseas and reaches maxima at water potential of wp_th
    --------------------------
    Parameters:
        wp: water potential at a daiy step 
        alfa: shape factor quantifying curve concavity; could be distinguished btw bacteria and fungi
        wp_fc: water potential at field capacity
        wp_th: water potential threshold
        
    Returns:
        f_osmo: scalar
        
    References:
        Modified from Manzoni et al. 2012 Ecology 
    """
    
    if wp >= wp_fc:
        f_osmo = 0.0
    elif wp <= wp_th:
        f_osmo = 1.0
    else:
        x = np.log(wp/wp_fc)
        y = np.log(wp_th/wp_fc)
        
        f_osmo = (x/y)**alfa

    return f_osmo


def microbe_mortality_prob(wp,wp_fc,death_rate,beta,tolerance):
    """
    microbial mortality probability as a function of water potential and drought tolerance
    
    Paramters:
        wp:    water potential;scalar
        wp_fc: field capacity;scalar
        death_rate: basal mortality prob. distinguished btw fungi and bacteria;array
        beta:       distinguished btw fungi and bacteria;array
        tolerance: taxon-specific drought tolerance;dataframe;values: 0 - 1
        
    Returns:
        mortality_rate: taxon-specific mortality probability
    
    References:
        Allison and Goulden,2017,Soil Biology and Biochemistry
    """
    
    if wp >= wp_fc:
        mortality_rate = death_rate
    else:
        tolerance = tolerance.to_numpy()
        # option 1
        mortality_rate = death_rate*(1 - beta*(wp-wp_fc)*(1-tolerance))
        # option 2
        mortality_rate = death_rate*(1 - beta*(wp-wp_fc)) * (1/np.exp(tolerance))

    
    return mortality_rate
    



    
    
        
