# output.py module dealing with outputs of DEMENTpy.
# Bin Wang in January, 2020

import numpy as np
import pandas as pd

class Output():
    """
    This class deals with outputs.
    
    Accepts data derived from the initialization.py and grid.py modules, and have two methods:
        output():             stores all time series
        microbes_abundance(): a special method
    """
    
    def __init__(self,runtime,data_init):
        """
        The constructor of Output class.
        
        Parameters:
            runtime:   user-specified parameters when running the model
            data_init: data dictionary;all inititialized data from the 'initialization.py' module   
        Returns:
            Initialization:   all data initialized preceding the execution of grid.py: dictionary
            Microbial_traits: microbial traits only pulled out from Initialization: dataframe
            SubstratesSeries: substrate-specific total mass over the grid: Substrate * C
            Substrates_Sum:   total substrates over the grid: day * C
            MonomersSeries:   monomoer-specific total mass over the grid: Monomer * C
            Monomers_Sum:     total monomers over the grid: day * C
            MicrobesSeries:   taxon-specific total biomass over the grid: Taxon * (C,N,P)
            Microbes_Sum:     total biomass over the grid: day * (C,N,P)
            TransporterSeries: taxon-specific total transporter production over the grid
            EnzymesSeries:    enzyme-specific total mass over the grid: Enzyme * C
            Enzymes_Sum:      total enzyme summed up over the grid: list
            OsmolyteSeries:   taxon-specific total osmolyte production over the grid
            RespSeries:       total respiration over the grid
            CUE:
        Notes:
            Variables that are not tracked but instead by DEMENT(R version):
                NH4Series:
                PO4Series:   
        """

        # A couple of vars used in processing outputs 
        n_taxa    = int(runtime.loc['n_taxa',1])                  # number of taxa
        Mic_index = ["Tax" + str(i) for i in range(1,n_taxa + 1)] # microbial taxa index
        

        # Pass all runtime parameters to Runtime
        self.Runtime = runtime
        # Pass all initialized data (a dictionary) to 'Initialization'
        self.Initialization = data_init
        
        # Pull out microbial traits data and put them in a dataframe:Microbial_traits
        data = np.concatenate((data_init['fb'][0:n_taxa],
                               data_init['Microbes_pp']['C'][0:n_taxa],
                               data_init['Microbes_pp']['N'][0:n_taxa],
                               data_init['Microbes_pp']['P'][0:n_taxa],
                               data_init['UptakeGenes'].sum(axis=1)[0:n_taxa],
                               data_init['EnzGenes'].sum(axis=1)[0:n_taxa],
                               data_init['OsmoGenes'].sum(axis=1)[0:n_taxa],
                               data_init['UptakeGenes_trait'][0:n_taxa],
                               data_init['EnzProdConsti_trait'][0:n_taxa],
                               data_init['EnzProdInduci_trait'][0:n_taxa],
                               data_init['OsmoProdConsti_trait'][0:n_taxa],
                               data_init['OsmoProdInduci_trait'][0:n_taxa],
                               data_init['TaxDroughtTol'][0:n_taxa])).reshape(13,n_taxa).T
        columns = ['F/B','C','N','P','Uptake_Gene','Enz_Gene','Osmo_Gene','Uptake_Cost',
                   'Enz_Consti_Cost','Enz_Induci_Cost','Osmo_Consti_Cost','Osmo_Induci_Cost',
                   'Drought_tolerance']
        self.Microbial_traits = pd.DataFrame(data=data, index=Mic_index, columns=columns)
        
        
        # Account for inputs in mass balance
        # self.Cum_Substrate      = (data_init['Substrates'].groupby(level=0,sort=False).sum()).sum(axis=0)
        # self.Cum_Monomer        = (data_init['Monomers'].groupby(level=0,sort=False).sum()).sum(axis=0)
        # self.Cum_Monomer_ratios = ecosystem.Cum_Monomer_ratios
        
        # Degradation rates
        #DecayRates_list       = [0] * n_substrates * gridsize
        #DecayRates_array      = np.array(DecayRates_list).reshape(n_substrates * gridsize,1)
        #DecayRates_df         = pd.DataFrame(data = DecayRates_array, index= data_init['Substrates'].index, columns= ['C'])
        #self.DecayRatesSeries = DecayRates_df.groupby(level=0,sort=False).sum()
        
        # Leaching
        #self.Cum_Leaching_N   = pd.Series([0],index=[0])
        #self.Cum_Leaching_P   = pd.Series([0],index=[0])

        # Substrates
        Substrates_grid           = data_init['Substrates'].groupby(level=0,sort=False).sum()
        Substrates_grid['C'].name = 0  # set the series name to 0
        self.SubstratesSeries     = Substrates_grid['C']
        #self.Substrates_Sum      = pd.Series([Substrates_grid['C'].sum()],index=[0])
        
        # Monomers
        Monomers_grid           = data_init['Monomers'].groupby(level=0,sort=False).sum()
        Monomers_grid['C'].name = 0
        self.MonomersSeries     = Monomers_grid["C"]
        #self.Monomers_Sum      = pd.Series([sum(Monomers_grid["C"])],index=[0])
        #self.NH4Series         = pd.Series([Monomers_grid.loc["NH4","N"]],index=[0])
        #self.PO4Series         = pd.Series([Monomers_grid.loc["PO4","P"]],index=[0])
        
        # Microbes
        Microbes_grid           = data_init['Microbes'].groupby(level=0,sort=False).sum() 
        Microbes_grid['C'].name = 0                                                     # Rename the series name to 0
        self.MicrobesSeries     = Microbes_grid['C']                                    # Taxon-specific total biomass summed over the grid
        #self.Microbes_Sum      = pd.Series([Microbes_grid['C'].sum(axis=0)],index=[0]) # Total biomass summed over the grid
        #self.Microbes_Interim  = Microbes_grid['C']                                    # Microbes right before executing metabolism calcualtion
        #self.MicrobesSeries_repop = Microbes_grid['C']                                 # Originally created for reinitialization

        # Number of individuals of Taxon count
        Taxon_index            = data_init['Microbes']['C'] > 0
        Taxon_index.name       = 0
        self.Taxon_count       = Taxon_index.groupby(level=0,sort=False).sum().astype('uint32')
        self.Taxon_count_repop = Taxon_index.groupby(level=0,sort=False).sum().astype('uint32')
        
        # Transporters: taxon-specific production summed over the grid by taxon
        #self.TransporterSeries = pd.Series(data=[0]*n_taxa,index=Mic_index)
        
        # Enyzmes: taxon-specific production summed over the grid by taxon
        #self.EnzymeConSeries    = pd.Series(data=[0]*n_taxa,index=Mic_index)
        #self.EnzymeIndSeries    = pd.Series(data=[0]*n_taxa,index=Mic_index)
        #self.Enzyme_TaxonSeries = pd.Series(data=[0]*n_taxa,index=Mic_index)
        
        # Osmolytes: taxon-specific production summed over the grid by taxon
        #self.OsmolyteConSeries    = pd.Series(data=[0]*n_taxa,index=Mic_index)
        #self.OsmolyteIndSeries    = pd.Series(data=[0]*n_taxa,index=Mic_index)
        #self.Osmolyte_TaxonSeries = pd.Series(data=[0]*n_taxa,index=Mic_index)
        
        # Growth yield
        #self.Growth_Yield = pd.Series(data=[0]*n_taxa,index=Mic_index)

        # Taxon-specific CUE
        #self.CUE_TaxonSeries = pd.Series(data=[0]*n_taxa,index=Mic_index)
        
        # Enzymes
        Enzymes_grid       = data_init['Enzymes'].groupby(level=0,sort=False).sum() # total of each enzyme summed over the spatial grid
        Enzymes_grid.name  = 0
        self.EnzymesSeries = Enzymes_grid 
        #self.Enzymes_Sum  = pd.Series([Enzymes_grid.sum()],index=[0])
        
        # Emergent properties over the grid
        self.RespSeries = pd.Series([0],index=[0], dtype='float32')  # respiration
        self.CUE_system = pd.Series([0],index=[0], dtype='float32')  # emergent CUE   
        self.Kill       = pd.Series([0],index=[0], dtype='uint32')   # stochastic death toll
        
       
    def output(self,ecosystem,day):
        """
        Records outputs in various variables of each iteration.

        Parameters:
            ecosystem: object from the grid.py module
            day:       the day to record to outputs
        Returns:
        """
        
        # Account for inputs in mass balance
        # self.Cum_Substrate = self.Cum_Substrate + ecosystem.SubInput.sum(axis = 0)
        # self.Cum_Monomer = self.Cum_Monomer + (ecosystem.MonInput.mul(ecosystem.Cum_Monomer_ratios,axis=0)).sum(axis=0)
        
        # DecayRates
        # DecayRates_grid = ecosystem.DecayRates.groupby(level=0,sort=False).sum()
        # self.DecayRatesSeries = pd.concat([self.DecayRatesSeries,DecayRates_grid],axis=1,sort=False)
        
        # Substrates
        Substrates_grid           = ecosystem.Substrates.groupby(level=0,sort=False).sum()
        Substrates_grid['C'].name = day + 1 # index the output by day
        self.SubstratesSeries     = pd.concat([self.SubstratesSeries,Substrates_grid['C']], axis=1, sort=False)
        #self.Substrates_Sum   = pd.concat([self.Substrates_Sum,pd.Series([Substrates_grid['C'].sum()],index=[day+1])],axis=0,sort=False)
        
        # Monomers
        Monomers_grid           = ecosystem.Monomers.groupby(level=0,sort=False).sum()
        Monomers_grid['C'].name = day + 1
        self.MonomersSeries     = pd.concat([self.MonomersSeries,Monomers_grid['C']], axis=1, sort=False)
        #self.Monomers_Sum   = pd.concat([self.Monomers_Sum, pd.Series([sum(Monomers_grid["C"])],  index=[day+1])], axis=0, sort=False)
        #self.NH4Series      = pd.concat([self.NH4Series, pd.Series([Monomers_grid.loc["NH4","N"]],index=[day+1])], axis=0, sort=False)
        #self.PO4Series      = pd.concat([self.PO4Series, pd.Series([Monomers_grid.loc["PO4","P"]],index=[day+1])], axis=0, sort=False)
        
        # Interim Microbes
        #Microbes_interim_grid = ecosystem.Microbes_interim.groupby(level=0,sort=False).sum()
        #self.Microbes_Interim = pd.concat([self.MicrobesSeries,Microbes_interim_grid['C']],axis=1,sort=False)
        # Count taxon for averaging taxon CUE
        #Taxon_index = (ecosystem.Microbes_interim)['C'] > 0
        #Taxon_index.name = day + 1
        #taxon_count = Taxon_index.groupby(level=0,sort=False).sum()
        #self.Taxon_count = pd.concat([self.Taxon_count,taxon_count],axis=1,sort=False)
        
        # Microbe
        ## Taxon abundance
        Taxon_index      = ecosystem.Microbes['C'] > 0
        Taxon_index.name = day + 1
        self.Taxon_count = pd.concat([self.Taxon_count, Taxon_index.groupby(level=0,sort=False).sum().astype('uint32')], axis=1, sort=False)
        ## Taxon biomass
        Microbes_grid           = ecosystem.Microbes.groupby(level=0,sort=False).sum()
        Microbes_grid['C'].name = day + 1
        self.MicrobesSeries     = pd.concat([self.MicrobesSeries, Microbes_grid['C']], axis=1, sort=False)
        #self.Microbes_Sum   = pd.concat([self.Microbes_Sum,pd.Series([Microbes_grid['C'].sum()],index=[day+1])],axis=0,sort=False)
        
        # Microbe-Transporters
        #Transporter_grid = ecosystem.Transporters.groupby(level=0,sort=False).sum() #taxon-specific transpotors summed over the grid by taxon
        #Transporter_grid.name = day + 1
        #self.TransporterSeries = pd.concat([self.TransporterSeries,Transporter_grid],axis=1,sort=False)
        
        # Microbe-Enzymes
        # Constitutive
        #Enzyme_Con_grid = ecosystem.Enzyme_Con.groupby(level=0,sort=False).sum() # Taxon-specific osmolytes summed over the grid by taxon
        #Enzyme_Con_grid.name = day + 1
        #self.EnzymeConSeries = pd.concat([self.EnzymeConSeries,Enzyme_Con_grid],axis=1,sort=False)
        # Inducible
        #Enzyme_Ind_grid = ecosystem.Enzyme_Ind.groupby(level=0,sort=False).sum()
        #Enzyme_Ind_grid.name = day + 1
        #self.EnzymeIndSeries = pd.concat([self.EnzymeIndSeries,Enzyme_Ind_grid],axis=1,sort=False)
        # Total
        #self.Enzyme_TaxonSeries = pd.concat([self.Enzyme_TaxonSeries,Enzyme_Con_grid+Enzyme_Ind_grid],axis=1,sort=False)
        
        # Microbe-Osmolytes
        # Constitutive
        #Osmolyte_Con_grid = ecosystem.Osmolyte_Con.groupby(level=0,sort=False).sum()
        #Osmolyte_Con_grid.name = day + 1
        #self.OsmolyteConSeries = pd.concat([self.OsmolyteConSeries,Osmolyte_Con_grid],axis=1,sort=False)
        # Inducible
        #Osmolyte_Ind_grid = ecosystem.Osmolyte_Ind.groupby(level=0,sort=False).sum()
        #Osmolyte_Ind_grid.name = day + 1
        #self.OsmolyteIndSeries = pd.concat([self.OsmolyteIndSeries,Osmolyte_Ind_grid],axis=1,sort=False)
        # Total
        #self.Osmolyte_TaxonSeries = pd.concat([self.Osmolyte_TaxonSeries,Osmolyte_Con_grid+Osmolyte_Ind_grid],axis=1,sort=False)
        
        # Growth yield by Taxon
        #GY_grid = ecosystem.Growth_Yield.groupby(level=0,sort=False).sum()
        #GY_grid.name = day + 1
        #self.Growth_Yield = pd.concat([self.Growth_Yield,GY_grid],axis=1,sort=False)
        
        # Microbe-Taxon-specific CUE
        #CUE_Taxon_grid = ecosystem.CUE_Taxon.groupby(level=0,sort=False).sum()
        #CUE_Taxon_grid = CUE_Taxon_grid/taxon_count
        #CUE_Taxon_grid.name = day + 1
        #self.CUE_TaxonSeries = pd.concat([self.CUE_TaxonSeries,CUE_Taxon_grid],axis=1,sort=False)
        
        # Enzymes over the grid
        # Derive the enzyme-specific enzyme production summed over the grid by enzyme
        Enzymes_grid       = ecosystem.Enzymes.groupby(level=0,sort=False).sum()
        Enzymes_grid.name  = day + 1
        self.EnzymesSeries = pd.concat([self.EnzymesSeries,Enzymes_grid], axis=1, sort=False)
        #self.Enzymes_Sum  = pd.concat([self.Enzymes_Sum,pd.Series([Enzymes_grid.sum()],index=[day+1])],axis=0,sort=False)
        
        # Respiration
        self.RespSeries = pd.concat([self.RespSeries, pd.Series([ecosystem.Respiration],index=[day+1],dtype='float32')], axis=0, sort=False)
        # Carbon use efficiency
        self.CUE_system = pd.concat([self.CUE_system, pd.Series([ecosystem.CUE_system], index=[day+1],dtype='float32')], axis=0, sort=False)
        # Death toll of stochasticity origin 
        self.Kill       = pd.concat([self.Kill, pd.Series([ecosystem.Kill],index=[day+1],dtype='uint32')], axis=0, sort=False)
    
    
    def microbes_abundance(self,ecosystem,day):
        """
        Seperately output Microbes from each time step and put them in a dataframe.

        Aims to eventually deal with reinitializing microbial community on the
        grid in a new pulse via calculating the cumulative frequency of different
        taxa over each cycle/pulse. The reason to have this seperate method instead of using the method
        above is b/c of the need to track outputs in every iteration, whereas the
        method above tracks outputs with a certain time interval (unless the interval set to 1)
        
        Parameters:
            ecosystem: an instance of the Grid object, in which only the Microbes is used.
            day:       the iteration index
        Returns:
            MicrobesSeries_repop: dataframe; tracking total biomass of differing taxa
            Taxon_count_repop:    dataframe; tracking abundances of differing taxa
        """
        # Track biomass of every taxon
        #Microbes_grid             = ecosystem.Microbes.groupby(level=0,sort=False).sum()
        #Microbes_grid['C'].name   = day + 1
        #self.MicrobesSeries_repop = pd.concat([self.MicrobesSeries_repop,Microbes_grid['C']],axis=1,sort=False)
        
        # Track abundance of every taxon
        Taxon_index      = ecosystem.Microbes['C'] > 0
        Taxon_index.name = day + 1
        taxon_count      = Taxon_index.groupby(level=0, sort=False).sum().astype('uint32')
        self.Taxon_count_repop = pd.concat([self.Taxon_count_repop,taxon_count], axis=1, sort=False)