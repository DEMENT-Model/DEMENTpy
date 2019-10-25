import numpy as np
import pandas as pd

class Output():

    """
    This module sets up the output properties by accepting data derived from
    the initialization.py and grid.py modules and uses two methods to record
    data generated from each iteration:
        output():
        microbes_df():
    ---------------------------------------------------------------------------
    Last modified by Bin Wang on 10/22/2019
    """
    
    def __init__(self,runtime,data_init):
        
        """
        Define and initialize variables in the Output class
        -----------------------------------------------------------------------
        Parameters:
            runtime:   user-specified parameters when running the model
            data_init: data dictionary;all inititialed data from the 'initialization.py' module
            
        Returns:
            Initialization:   all data initialized preceding the execution of grid.py: dictionary
            Microbial_traits: microbial traits only pulled out from Initialization: dataframe
            SubstratesSeries: substrate-specific total mass over the grid: Substrate * C
            Substrates_Sum:   total substrates over the grid: day * C
            MonomersSeries: monomoer-specific total mass over the grid: Monomer * C
            Monomers_Sum:   total monomers over the grid: day * C
            NH4Series:
            PO4Series:  
            MicrobesSeries: taxon-specific total biomass over the grid: Taxon * (C,N,P)
            Microbes_Sum:   total biomass over the grid: day * (C,N,P)
            TransporterSeries: taxon-specific total transporter production over the grid
            EnzymesSeries: enzyme-specific total mass over the grid: Enzyme * C
            Enzymes_Sum:   total enzyme summed up over the grid: list
            OsmolyteSeries: taxon-specific total osmolyte production over the grid
            RespSeries: total respiration over the grid
            CUE:
        """
        
        # A few vars used in processing outputs 
        # n_substrates = int(runtime.loc['n_substrates',1])
        # self.gridsize = int(runtime.loc['gridsize',1])
        n_taxa = int(runtime.loc['n_taxa',1])
        Mic_index = ["Tax" + str(i) for i in range(1,n_taxa + 1)] # Microbial taxa index
        
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
        self.Microbial_traits = pd.DataFrame(data=data,index=Mic_index,columns=columns)
        
        
        # Account for inputs in mass balance
        # self.Cum_Substrate = (data_init['Substrates'].groupby(level=0,sort=False).sum()).sum(axis=0)
        # self.Cum_Monomer = (data_init['Monomers'].groupby(level=0,sort=False).sum()).sum(axis=0)
        # self.Cum_Monomer_ratios = ecosystem.Cum_Monomer_ratios
        
        
        # Degradation rates
        #DecayRates_list = [0] * n_substrates * gridsize
        #DecayRates_array= np.array(DecayRates_list).reshape(n_substrates * gridsize,1)
        #DecayRates_df   = pd.DataFrame(data = DecayRates_array,
        #                               index= data_init['Substrates'].index,
        #                               columns= ['C'])
        #self.DecayRatesSeries = DecayRates_df.groupby(level=0,sort=False).sum()
        
        
        # Substrates
        Substrates_grid = data_init['Substrates'].groupby(level=0,sort=False).sum()
        Substrates_grid['C'].name = 0 # reset the index to 0
        self.SubstratesSeries = Substrates_grid['C']
        self.Substrates_Sum   = pd.Series([Substrates_grid['C'].sum()],index=[0])
        
        # Monomers
        Monomers_grid = data_init['Monomers'].groupby(level=0,sort=False).sum()
        self.MonomersSeries   = Monomers_grid["C"]
        self.Monomers_Sum     = pd.Series([sum(Monomers_grid["C"])],index=[0])
        self.NH4Series        = pd.Series([Monomers_grid.loc["NH4","N"]],index=[0])
        self.PO4Series        = pd.Series([Monomers_grid.loc["PO4","P"]],index=[0])
        
        # Leaching
        #self.Cum_Leaching_N   = pd.Series([0],index=[0])
        #self.Cum_Leaching_P   = pd.Series([0],index=[0])
        
        
        # Total microbial biomass summed over the grid by taxon
        Microbes_grid = data_init['Microbes'].groupby(level=0,sort=False).sum()
        Microbes_grid['C'].name = 0
        
        self.MicrobesSeries_repop = Microbes_grid['C'] #
        self.MicrobesSeries = Microbes_grid['C']       #
        self.Microbes_Interim = Microbes_grid['C']     # Microbes right before executing metabolism calcualtion
        self.Microbes_Sum = pd.Series([Microbes_grid['C'].sum(axis=0)],index=[0])
        
        # Number of individuals of Taxon count
        Taxon_index = data_init['Microbes']['C'] > 0
        Taxon_index.name = 0
        self.Taxon_count = Taxon_index.groupby(level=0,sort=False).sum()
        self.Taxon_count_repop = Taxon_index.groupby(level=0,sort=False).sum()
        
        # Transporters: taxon-specific production summed over the grid by taxon
        self.TransporterSeries = pd.Series(data=[0]*n_taxa,index=Mic_index)
        
        # Enyzmes: taxon-specific production summed over the grid by taxon
        self.EnzymeConSeries = pd.Series(data=[0]*n_taxa,index=Mic_index)
        self.EnzymeIndSeries = pd.Series(data=[0]*n_taxa,index=Mic_index)
        self.Enzyme_TaxonSeries = pd.Series(data=[0]*n_taxa,index=Mic_index)
        
        # Osmolytes: taxon-specific production summed over the grid by taxon
        self.OsmolyteConSeries = pd.Series(data=[0]*n_taxa,index=Mic_index)
        self.OsmolyteIndSeries = pd.Series(data=[0]*n_taxa,index=Mic_index)
        self.Osmolyte_TaxonSeries = pd.Series(data=[0]*n_taxa,index=Mic_index)
        
        # Growth yield
        self.Growth_Yield =pd.Series(data=[0]*n_taxa,index=Mic_index)

        # Taxon-specific CUE
        self.CUE_TaxonSeries = pd.Series(data=[0]*n_taxa,index=Mic_index)
        
        
        # Enzymes
        # Total enzyme summed over the spatial grid by enzyme 
        Enzymes_grid = data_init['Enzymes'].groupby(level=0,sort=False).sum()
        self.EnzymesSeries = Enzymes_grid 
        self.Enzymes_Sum   = pd.Series([Enzymes_grid.sum()],index=[0])
        
        
        # Emergent properties over the grid
        self.RespSeries         = pd.Series([0],index=[0])
        self.CUE_System_Series  = pd.Series([0],index=[0])
        self.Kill               = pd.Series([0],index=[0]) # number of taxon killed
        
        
        
    def output(self,ecosystem,day):
        
        """
        This method records outputs in various variables at a daily time step.
        
        Input:
            ecosystem: object from the grid.py module
            day:       the day to record to outputs
        """
        
        
        # Account for inputs in mass balance
        # self.Cum_Substrate = self.Cum_Substrate + ecosystem.SubInput.sum(axis = 0)
        # self.Cum_Monomer = self.Cum_Monomer + (ecosystem.MonInput.mul(ecosystem.Cum_Monomer_ratios,axis=0)).sum(axis=0)
        
        # DecayRates
        # DecayRates_grid = ecosystem.DecayRates.groupby(level=0,sort=False).sum()
        # self.DecayRatesSeries = pd.concat([self.DecayRatesSeries,DecayRates_grid],axis=1,sort=False)
        
        
        # Substrates
        Substrates_grid = ecosystem.Substrates.groupby(level=0,sort=False).sum()
        Substrates_grid['C'].name = day + 1 # index the output by day
        self.SubstratesSeries = pd.concat([self.SubstratesSeries,Substrates_grid['C']],axis=1,sort=False)
        self.Substrates_Sum = pd.concat([self.Substrates_Sum,pd.Series([Substrates_grid['C'].sum()],index=[day+1])],axis=0,sort=False)
        
        
        # Monomers: note leaching is missing
        Monomers_grid = ecosystem.Monomers.groupby(level=0,sort=False).sum()
        self.MonomersSeries = pd.concat([self.MonomersSeries,Monomers_grid["C"]],axis=1,sort=False)
        self.NH4Series = pd.concat([self.NH4Series,pd.Series([Monomers_grid.loc["NH4","N"]],index=[day+1])],axis=0,sort=False)
        self.PO4Series = pd.concat([self.PO4Series,pd.Series([Monomers_grid.loc["PO4","P"]],index=[day+1])],axis=0,sort=False)
        self.Monomers_Sum = pd.concat([self.Monomers_Sum,pd.Series([sum(Monomers_grid["C"])],index=[day+1])],axis=0,sort=False)
        
        
        # Interim Microbes
        Microbes_interim_grid = ecosystem.Microbes_interim.groupby(level=0,sort=False).sum()
        self.Microbes_Interim = pd.concat([self.MicrobesSeries,Microbes_interim_grid['C']],axis=1,sort=False)
        # Count taxon for averaging taxon CUE
        Taxon_index = (ecosystem.Microbes_interim)['C'] > 0
        Taxon_index.name = day + 1
        taxon_count = Taxon_index.groupby(level=0,sort=False).sum()
        self.Taxon_count = pd.concat([self.Taxon_count,taxon_count],axis=1,sort=False)
        
        # Microbe-Transporters
        # Derive the taxon-specific transpotors summed over the grid by taxon
        Transporter_grid = ecosystem.Transporters.groupby(level=0,sort=False).sum()
        #Transporter_grid = Transporter_grid/Microbes_w_grid['C']
        Transporter_grid.name = day + 1
        self.TransporterSeries = pd.concat([self.TransporterSeries,Transporter_grid],axis=1,sort=False)
        
        # Microbe-Enzymes
        # Derive the taxon-specific osmolytes summed over the grid by taxon
        # Constitutive
        Enzyme_Con_grid = ecosystem.Enzyme_Con.groupby(level=0,sort=False).sum()
        #Enzyme_Con_grid = Enzyme_Con_grid/Microbes_w_grid['C']
        Enzyme_Con_grid.name = day + 1
        self.EnzymeConSeries = pd.concat([self.EnzymeConSeries,Enzyme_Con_grid],axis=1,sort=False)
        # Inducible
        Enzyme_Ind_grid = ecosystem.Enzyme_Ind.groupby(level=0,sort=False).sum()
        #Enzyme_Ind_grid = Enzyme_Ind_grid/Microbes_w_grid['C']
        Enzyme_Ind_grid.name = day + 1
        self.EnzymeIndSeries = pd.concat([self.EnzymeIndSeries,Enzyme_Ind_grid],axis=1,sort=False)
        # Total
        self.Enzyme_TaxonSeries = pd.concat([self.Enzyme_TaxonSeries,Enzyme_Con_grid+Enzyme_Ind_grid],axis=1,sort=False)
        
        # Microbe-Osmolytes
        # Derive the taxon-specific osmolytes summed over the grid by taxon
        # Constitutive
        Osmolyte_Con_grid = ecosystem.Osmolyte_Con.groupby(level=0,sort=False).sum()
        #Osmolyte_Con_grid = Osmolyte_Con_grid/Microbes_w_grid['C']
        Osmolyte_Con_grid.name = day + 1
        self.OsmolyteConSeries = pd.concat([self.OsmolyteConSeries,Osmolyte_Con_grid],axis=1,sort=False)
        # Inducible
        Osmolyte_Ind_grid = ecosystem.Osmolyte_Ind.groupby(level=0,sort=False).sum()
        #Osmolyte_Ind_grid = Osmolyte_Ind_grid/Microbes_w_grid['C']
        Osmolyte_Ind_grid.name = day + 1
        self.OsmolyteIndSeries = pd.concat([self.OsmolyteIndSeries,Osmolyte_Ind_grid],axis=1,sort=False)
        # Total
        self.Osmolyte_TaxonSeries = pd.concat([self.Osmolyte_TaxonSeries,Osmolyte_Con_grid+Osmolyte_Ind_grid],axis=1,sort=False)
        
        # Growth yield by Taxon
        GY_grid = ecosystem.Growth_Yield.groupby(level=0,sort=False).sum()
        GY_grid.name = day + 1
        self.Growth_Yield = pd.concat([self.Growth_Yield,GY_grid],axis=1,sort=False)
        
        
        # Microbe-Taxon-specific CUE
        CUE_Taxon_grid = ecosystem.CUE_Taxon.groupby(level=0,sort=False).sum()
        CUE_Taxon_grid = CUE_Taxon_grid/taxon_count
        CUE_Taxon_grid.name = day + 1
        self.CUE_TaxonSeries = pd.concat([self.CUE_TaxonSeries,CUE_Taxon_grid],axis=1,sort=False)
        
        
        # Microbe
        Microbes_grid = ecosystem.Microbes.groupby(level=0,sort=False).sum()
        
        # for outputs
        Microbes_grid['C'].name = day + 1
        self.MicrobesSeries = pd.concat([self.MicrobesSeries,Microbes_grid['C']],axis=1,sort=False)
        self.Microbes_Sum = pd.concat([self.Microbes_Sum,pd.Series([Microbes_grid['C'].sum()],index=[day+1])],axis=0,sort=False)
        
        # Enzymes over the grid
        # Derive the enzyme-specific enzyme production summed over the grid by enzyme
        Enzymes_grid = ecosystem.Enzymes.groupby(level=0,sort=False).sum()
        Enzymes_grid.name = day + 1
        self.EnzymesSeries = pd.concat([self.EnzymesSeries,Enzymes_grid],axis=1,sort=False)
        self.Enzymes_Sum = pd.concat([self.Enzymes_Sum,pd.Series([Enzymes_grid.sum()],index=[day+1])],axis=0,sort=False)
        
        # Respiration
        #self.RespSeries.append(pd.Series([ecosystem.Respiration],index=[day]))
        self.RespSeries = pd.concat([self.RespSeries,pd.Series([ecosystem.Respiration],index=[day+1])],axis=0,sort=False)
        # Carbon use efficiency
        self.CUE_System_Series = pd.concat([self.CUE_System_Series,pd.Series([ecosystem.CUE_System],index=[day+1])],axis=0,sort=False)
        
        self.Kill = pd.concat([self.Kill,pd.Series([ecosystem.Kill],index=[day+1])],axis=0,sort=False)
    
    
    def microbes_df(self,ecosystem,day):
        
        """
        seperately output Microbes from each time step and put them in a dataframe:
            MicrobesSeries_repop
        
        aiming to eventually deal with relocating microbial community on the grid after a pulse.
        via computing the average frequency of different taxa over each cycle,
        based on which initializing a new community on the grid at the very start of a new cycle
        
        The reason of having this seperate method instead of using the method above is b/c of
        the need to track very timestep's output, whereas the method above only track outputs wiith a certain interval.
        -----------------------------------------------
        input:
            ecosystem: an instance of the Grid object, in which the Microbes property is intended to be used here.
        """
        
        Microbes_grid = ecosystem.Microbes.groupby(level=0,sort=False).sum()
        Microbes_grid['C'].name = day + 1
        self.MicrobesSeries_repop = pd.concat([self.MicrobesSeries_repop,Microbes_grid['C']],axis=1,sort=False)
        
        # Count each taxon
        Taxon_index = (ecosystem.Microbes)['C'] > 0
        Taxon_index.name = day + 1
        taxon_count = Taxon_index.groupby(level=0,sort=False).sum()
        self.Taxon_count_repop = pd.concat([self.Taxon_count_repop,taxon_count],axis=1,sort=False)
        
        
        
        