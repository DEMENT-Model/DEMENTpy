# ------------------
# this module, grid.py, deal with calculations of all microbe-related activites on a spatial grid with a class, Grid().
# by Bin Wang on December 27th, 2019 
# ------------------

import sys
import numpy  as np
import pandas as pd

from microbe import microbe_osmo_psi
from microbe import microbe_mortality_prob as MMP
from enzyme  import Arrhenius
from utility import expand

class Grid():
    """
    This class holds all variables related to microbe, substrate, monomer, and enzyme over the spatial grid.
    
    Accepts returns from the module 'initialization.py' and includes methods as follows:
        1) degrdation():   explicit substrate degradation
        2) uptake():       explicit monomers uptake
        3) metabolism():   cellular processes and emergent CUE and respiration
        4) mortality():    determine mortality of microbial cells based on mass thresholds
        5) reproduction(): compute cell division and dispersal
        6) repopulation(): resample taxa from the microbial pool and place them on the grid
    Coding philosophy:
        Each method starts with passing some global variables to local ones and creating
        some indices facilitating dataframe index/column processing and ends up with updating
        state variables and passing them back to the global ones. All computation stays in between.   
    Reminder:
        Keep a CLOSE EYE on the indexing throughout the matrix/dataframe operations
    """
    
    def __init__(self,runtime,data_init): 
        """
        The constructor of Grid class.

        Parameters:
            runtime:   user-specified parameters
            data_init: dictionary;initialized data from the module 'initialization.py'
        """

        self.cycle          = int(runtime.loc['end_time',1])
        self.gridsize       = int(runtime.loc['gridsize',1])
        self.n_taxa         = int(runtime.loc["n_taxa",1])
        self.n_substrates   = int(runtime.loc["n_substrates",1])
        self.n_enzymes      = int(runtime.loc["n_enzymes",1])
        self.n_monomers     = self.n_substrates + 2
        
        #Degradation
        self.Substrates_init = data_init['Substrates']                    # Substrates initialized
        self.Substrates      = data_init['Substrates'].copy(deep=True)    # Substrates;df; w/ .copy() avoiding mutation
        self.SubInput        = data_init['SubInput']                      # Substrate inputs
        self.Enzymes_init    = data_init['Enzymes']                       # Initial pool of Enzymes
        self.Enzymes         = data_init['Enzymes'].copy(deep=True)       # Enzymes
        self.ReqEnz          = data_init['ReqEnz']                        # Enzymes required by each substrate
        self.Ea              = data_init['Ea']                            # Enzyme activatin energy
        self.Vmax0           = data_init['Vmax0']                         # Max. reaction speed
        self.Km0             = data_init['Km0']                           # Half-saturation constant
        self.SubstrateRatios = np.float32('nan')                          # Substrate stoichiometry
        self.DecayRates      = np.float32('nan')                          # Substrate decay rate

        #Uptake
        self.Microbes_init   = data_init['Microbes_pp']                   # microbial community before placement
        self.Microbes        = data_init['Microbes'].copy(deep=True)      # microbial community after placement
        self.Monomers_init   = data_init['Monomers']                      # Monomers initialized
        self.Monomers        = data_init['Monomers'].copy(deep=True)      # Monomers
        self.MonInput        = data_init['MonInput']                      # Inputs of monomers
        self.Uptake_Ea       = data_init['Uptake_Ea']                     # transporter enzyme Ea
        self.Uptake_Vmax0    = data_init['Uptake_Vmax0']                  # transporter Vmax
        self.Uptake_Km0      = data_init['Uptake_Km0']                    # transporter Km
        self.Monomer_ratios  = data_init['Monomer_ratio'].copy(deep=True) # monomer stoichiometry
        self.Uptake_ReqEnz   = data_init['Uptake_ReqEnz']                 # Enzymes required by monomers 
        self.Uptake_Enz_Cost = data_init['UptakeGenesCost']               # Cost of encoding each uptake gene
        self.Taxon_Uptake_C  = np.float32('nan')                          # taxon uptake of C 
        self.Taxon_Uptake_N  = np.float32('nan')                          # taxon uptake of N 
        self.Taxon_Uptake_P  = np.float32('nan')                          # taxon uptake of P
        
        #Metabolism
        self.Consti_Enzyme_C   = data_init["EnzProdConstit"]    # C cost of encoding constitutive enzyme
        self.Induci_Enzyme_C   = data_init["EnzProdInduce"]     # C Cost of encoding inducible enzyme 
        self.Consti_Osmo_C     = data_init['OsmoProdConsti']    # C Cost of encoding constitutive osmolyte
        self.Induci_Osmo_C     = data_init['OsmoProdInduci']    # C Cost of encoding inducible osmolyte 
        self.Uptake_Maint_Cost = data_init['Uptake_Maint_cost'] # Respiration cost of uptake transporters: 0.01	mg C transporter-1 day-1     
        self.Enz_Attrib        = data_init['EnzAttrib']         # Enzyme attributes; dataframe
        self.AE_ref            = data_init['AE_ref']            # Reference AE:constant of 0.5;scalar
        self.AE_temp           = data_init['AE_temp']           # AE sensitivity to temperature;scalar
        self.Respiration       = np.float32('nan')              # Respiration
        self.CUE_system        = np.float32('nan')              # emergent CUE
        #self.Transporters = float('nan')
        #self.Osmolyte_Con = float('nan')
        #self.Osmolyte_Ind = float('nan')
        #self.Enzyme_Con   = float('nan')
        #self.Enzyme_Ind   = float('nan')
        #self.CUE_Taxon    = float('nan')
        #self.Growth_Yield = float('nan')

        #Mortality
        self.MinRatios        = data_init['MinRatios']           # Minimal cell quotas
        self.C_min            = data_init['C_min']               # C threshold value of living cell
        self.N_min            = data_init['N_min']               # N threshold value of living cell
        self.P_min            = data_init['P_min']               # P threshold value of living cell
        self.basal_death_prob = data_init['basal_death_prob']    # Basal death probability of microbes
        self.death_rate       = data_init['death_rate']          # Change rate of mortality with water potential
        self.tolerance        = data_init['TaxDroughtTol']       # taxon drought tolerance
        self.wp_fc            = data_init['wp_fc']               # -1.0
        self.wp_th            = data_init['wp_th']               # -6.0
        self.alpha            = data_init['alpha']               # integer;1
        self.Kill             = np.float32('nan')                # total number of cells stochastically killed
        
        # Reproduction
        self.fb         =  data_init['fb']                       # index of fungal taxa (=1)
        self.max_size_b =  data_init['max_size_b']               # threshold of cell division
        self.max_size_f =  data_init['max_size_f']               # threshold of cell division
        self.x          =  int(runtime.loc['x',1])               # x dimension of grid
        self.y          =  int(runtime.loc['y',1])               # y dimension of grid
        self.dist       =  int(runtime.loc['dist',1])            # maximum dispersal distance: 1 cell
        self.direct     =  int(runtime.loc['direct',1])          # dispersal direction: 0.95
        
        # Climate data
        self.temp = data_init['Temp']       # series; temperature
        self.psi  = data_init['Psi']        # series; water potential
        
        # Global constants
        self.Km_Ea = np.float32(20)         # kj mol-1;activation energy for both enzyme and transporter
        self.Tref  = np.float32(293)        # reference temperature of 20 celcius
    

    def degradation(self,pulse,day):
        """
        Explicit degradation of different substrates.

        Calculation procedure:
        1. Determine substates pool: incl. inputs
        2. Compute Vmax & Km and make them follow the index of Substrates
        3. Follow the Michaelis-Menten equation to compute full degradation rate
        4. Impose the substrate-required enzymes upon the full degradation rate
        5. Adjust cellulose rate with LCI(lignocellulose index)
        """
        
        # constant of lignocellulose index--LCI
        LCI_slope = np.float32(-0.8)  
        # Substrates index by subtrate names
        Sub_index = self.Substrates.index
        
        
        # total mass of each substrate: C+N+P
        rss = self.Substrates.sum(axis=1) 
        # substrate stoichiometry; NOTE:ensure NA(b/c of 0/0 in df) = 0
        SubstrateRatios = self.Substrates.divide(rss,axis=0)
        SubstrateRatios = SubstrateRatios.fillna(0)
        
        # moisture effects on enzymatic kinetics
        if self.psi[day] >= self.wp_fc:
            f_psi = np.float32(1.0)
        else:
            f_psi = np.exp(0.25*(self.psi[day] - self.wp_fc)).astype('float32')
        
        # Arrhenius equation for Vmax and Km multiplied by exponential decay for Psi sensitivity
        #Vmax = self.Vmax0 * np.exp((-self.Ea/0.008314)*(1/(self.temp[day]+273) - 1/self.Tref)) * f_psi  # Vmax: (enz*gridsize) * sub
        #Km   = self.Km0 * np.exp((-self.Km_Ea/0.008314)*(1/(self.temp[day]+273) - 1/self.Tref))         # Km: (sub*gridsize) * enz
        Vmax = self.Vmax0 * Arrhenius(self.Ea,   self.temp[day]) * f_psi  # Vmax: (enz*gridsize) * sub
        Km   = self.Km0   * Arrhenius(self.Km_Ea,self.temp[day])          # Km:   (sub*gridsize) * enz

        # Multiply Vmax by enzyme concentration
        tev_transition       = Vmax.mul(self.Enzymes,axis=0)                                          # (enz*gridsize) * sub
        tev_transition.index = [np.arange(self.gridsize).repeat(self.n_enzymes),tev_transition.index] # create a MultiIndex
        tev                  = tev_transition.stack().unstack(1).reset_index(level=0,drop=True)       # (sub*gridsize) * enz
        tev                  = tev[Km.columns]                                                        # ensure to re-order the columns b/c of python's default alphabetical ordering

        # Michaelis-Menten equation
        Decay = tev.mul(rss,axis=0)/Km.add(rss,axis=0)
        
        # Pull out each batch of required enzymes and sum across redundant enzymes
        batch1 = (self.ReqEnz.loc['set1'].values * Decay).sum(axis=1)
        #batch2 = (self.ReqEnz.loc['set2'].values * Decay).sum(axis=1)
        
        # Assess the rate-limiting enzyme and set decay to that rate
        #DecaySums = pd.concat([batch1, batch2],axis=1)
        #DecayRates0 = DecaySums.min(axis=1, skipna=True)
        
        # Compare to substrate available and take the min, allowing for a tolerance of 1e-9
        DecayRates = pd.concat([batch1,rss],axis=1,sort=False).min(axis=1,skipna=True)
        
        # Adjust cellulose rate by linking cellulose degradation to lignin concentration (LCI) 
        ss7 = self.Substrates.loc[Sub_index=='Lignin'].sum(axis=1).values
        DecayRates.loc[Sub_index=='Cellulose'] *= np.float32(1) + (ss7/(ss7 + self.Substrates.loc[Sub_index=='Cellulose','C'])) * LCI_slope
        
        # Update Substrates Pool by removing decayed C, N, & P. Depending on specific needs, adding inputs of substrates can be done here
        self.Substrates -= SubstrateRatios.mul(DecayRates,axis=0) #+ self.SubInput 
        
        # Pass these two back to the global variables to be used in the next method
        self.SubstrateRatios = SubstrateRatios
        self.DecayRates      = DecayRates


    def uptake(self,pulse,day):
        """
        Explicit uptake of different monomers by transporters following the Michaelis-Menten equation.

        Calculaton procedure:
        -> Average monomers across the grid:
        -> Determine pool of monomers: add degradation and input, update stoichimoetry
        -> Maximum uptake:
        -> Uptake by Monomer:
        -> Uptake by Taxon:
        """
        
        # Every monomer averaged over the grid in each time step
        self.Monomers = expand(self.Monomers.groupby(level=0,sort=False).sum()/self.gridsize,self.gridsize)
        

        # Indices
        is_org = (self.Monomers.index != "NH4") & (self.Monomers.index != "PO4") # organic monomers
        #is_mineral = (Monomers.index == "NH4") | (Monomers.index == "PO4")

        # Update monomer ratios in each time step with organic monomers following the substrates
        self.Monomer_ratios[is_org] = self.SubstrateRatios.values
        
        # Determine monomer pool from decay and input
        # Organic monomers derived from substrate-decomposition
        Decay_Org = self.Monomer_ratios[is_org].mul(self.DecayRates.values,axis=0)
        # inputs of organic and mineral monomers
        #Input_Org = MR_transition[is_org].mul(self.MonInput[is_org].tolist(),axis=0)
        #Input_Mineral = MR_transition[is_mineral].mul((self.MonInput[is_mineral]).tolist(),axis=0)
        # Monomer pool determined
        self.Monomers.loc[is_org] += Decay_Org #+ Input_Org
        #self.Monomers.loc[is_mineral] += Input_Mineral
        
        # Get the total mass of each monomer: C+N+P
        rsm = self.Monomers.sum(axis=1)
        # Recalculate monomer ratios after updating monomer pool and before uptake calculation
        self.Monomer_ratios.loc[is_org] = self.Monomers.loc[is_org].divide(rsm[is_org],axis=0)
        self.Monomer_ratios             = self.Monomer_ratios.fillna(0)


        # Start calculating monomer uptake
        # Moisture impacts on uptake, mimicking the diffusivity implications
        if self.psi[day] >= self.wp_fc:
            f_psi = np.float32(1.0)
        else:
            f_psi = np.exp(0.5*(self.psi[day] - self.wp_fc)).astype('float32')
        
        # Caculate uptake enzyme kinetic parameters; monomer * Upt
        #Uptake_Vmax = self.Uptake_Vmax0 * np.exp((-self.Uptake_Ea/0.008314)*(1/(self.temp[day]+273) - 1/self.Tref)) * f_psi
        #Uptake_Km   = self.Uptake_Km0 * np.exp((-self.Km_Ea/0.008314)*(1/(self.temp[day]+273) - 1/self.Tref))
        Uptake_Vmax = self.Uptake_Vmax0 * Arrhenius(self.Uptake_Ea,self.temp[day]) * f_psi
        Uptake_Km   = self.Uptake_Km0   * Arrhenius(self.Km_Ea,self.temp[day])

        # Equation for hypothetical potential uptake (per unit of compatible uptake protein)
        Potential_Uptake = (self.Uptake_ReqEnz * Uptake_Vmax).mul(rsm.values,axis=0)/Uptake_Km.add(rsm.values,axis=0)
        
        # Derive the mass of each transporter of each taxon NOTE: transpose the df to Upt*(Taxa*grid)
        MicCXGenes = (self.Uptake_Enz_Cost.mul(self.Microbes.sum(axis=1),axis=0)).T
        
        # Define Max_Uptake: (Monomer*gridsize) * Taxon
        Max_Uptake_array = np.array([0]*self.gridsize*self.n_monomers*self.n_taxa).reshape(self.gridsize*self.n_monomers,self.n_taxa).astype('float32')
        Max_Uptake       = pd.DataFrame(data=Max_Uptake_array, index=self.Monomers.index, columns=self.Microbes.index[0:self.n_taxa])
        # Matrix multiplication to get max possible uptake by monomer(extract each grid point separately for operation)
        for i in range(self.gridsize):
            i_monomer = np.arange(i * self.n_monomers, (i+1) * self.n_monomers)
            i_taxa    = np.arange(i * self.n_taxa,     (i+1) * self.n_taxa)
            Max_Uptake.iloc[i_monomer,:] = Potential_Uptake.iloc[i_monomer,:].values @ MicCXGenes.iloc[:,i_taxa].values
        
        # Take the min of the monomer available and the max potential uptake, and scale the uptake to what's available
        csmu                = Max_Uptake.sum(axis=1)  # total potential uptake of each monomer
        Uptake              = Max_Uptake.mul(pd.concat([csmu,rsm],axis=1).min(axis=1,skipna=True)/csmu,axis=0) #(Monomer*gridsize) * Taxon
        Uptake.loc[csmu==0] = np.float32(0)
        # End computing monomer uptake
        

        # Update Monomers
        # By monomer: total uptake (monomer*gridsize) * 3(C-N-P)
        self.Monomers -= self.Monomer_ratios.mul(Uptake.sum(axis=1),axis=0)
        
        # Derive Taxon-specific total uptake of C, N, & P
        # By taxon: total uptake; (monomer*gridsize) * taxon
        C_uptake_df = Uptake.mul(self.Monomer_ratios["C"],axis=0)
        N_uptake_df = Uptake.mul(self.Monomer_ratios["N"],axis=0)
        P_uptake_df = Uptake.mul(self.Monomer_ratios["P"],axis=0)
        # generic multi-index
        C_uptake_df.index = N_uptake_df.index = P_uptake_df.index = [np.arange(self.gridsize).repeat(self.n_monomers),C_uptake_df.index]
        TUC_df = C_uptake_df.groupby(level=[0]).sum()
        TUN_df = N_uptake_df.groupby(level=[0]).sum()
        TUP_df = P_uptake_df.groupby(level=[0]).sum()
        # Update these 3 global variables
        self.Taxon_Uptake_C = TUC_df.stack().values     # spatial C uptake: array
        self.Taxon_Uptake_N = TUN_df.stack().values     # spatial N uptake: array
        self.Taxon_Uptake_P = TUP_df.stack().values     # spatial P uptake: array

        
    def metabolism(self,day):
        """
        Explicitly calculate intra-cellular production of metabolites.
        
        Handles both constitutive (standing biomass) and inducible (immediate monomers uptake) pathways following:
          1. constitutive enzyme and osmolyte production
          2. inducible enzyme and osmolyte production
          3. emergent CUE & Respiration
          4. update both Enzymes (with production & loss) and Substrates (with dead enzymes)
        """
        
        # index of dead enzyme in Substrates
        is_deadEnz = self.Substrates.index == "DeadEnz"
        # Constants
        Osmo_N_cost      = np.float32(0.3)   # N cost per unit of osmo-C production 
        Osmo_Maint_cost  = np.float32(5.0)   # C loss per unit of osmo-C production
        Enzyme_Loss_Rate = np.float32(0.04)  # enzyme turnover rate(=0.04; Allison 2006)

        # Scalar of water potential impact: call the function microbe_osmo_psi()
        f_psi = microbe_osmo_psi(self.psi[day],self.alpha,self.wp_fc,self.wp_th)

        #---------------------------------------------------------------------#
        #......................constitutive processes.........................#
        #---------------------------------------------------------------------#
        # Transporters' maintenence
        # taxon-specific uptake cost determined by total biomass C: 0.1 - 0.01
        # Taxon_Transporter_Cost  = (self.Uptake_Enz_Cost.mul(self.Microbes['C'],axis=0)).sum(axis=1) #NOTE Microbes['C'] vs Microbes.sum(axis=1)
        # taxon-specific respiration cost of producing transporters: self.uptake_maint_cost = 0.01
        Taxon_Transporter_Maint = (self.Uptake_Enz_Cost.mul(self.Microbes['C'],axis=0)).sum(axis=1) * self.Uptake_Maint_Cost
        
        #...............................................
        # Variable Acronyms:
        # OECCN : Osmo_Enzyme_Consti_Cost_N
        # ARROEC: Avail_Req_ratio_osmo_enzyme_consti
        # MNAOEC: Min_N_Avail_Osmo_Enzyme_Consti
        #...............................................
        # Osmolyte before adjustment
        Taxon_Osmo_Consti          = self.Consti_Osmo_C.mul(self.Microbes['C'],axis=0) * f_psi
        Taxon_Osmo_Consti_Cost_N   = (Taxon_Osmo_Consti * Osmo_N_cost).sum(axis=1)
        # Enzyme before adjustment
        Taxon_Enzyme_Consti        = self.Consti_Enzyme_C.mul(self.Microbes['C'],axis=0)
        Taxon_Enzyme_Consti_Cost_N = (Taxon_Enzyme_Consti.mul(self.Enz_Attrib['N_cost'],axis=1)).sum(axis=1)
        # Adjust osmolyte & enzyme production based on available N in microbial biomass
        OECCN  = Taxon_Osmo_Consti_Cost_N + Taxon_Enzyme_Consti_Cost_N                                    # Total N cost
        MNAOEC = (pd.concat([OECCN[OECCN>0],self.Microbes['N'][OECCN>0]],axis=1)).min(axis=1,skipna=True) # get the minimum value
        ARROEC = (MNAOEC/OECCN[OECCN>0]).fillna(0)                                                        # Derive ratio of availabe N to required N
        # Osmolyte adjusted
        Taxon_Osmo_Consti[OECCN>0] = Taxon_Osmo_Consti[OECCN>0].mul(ARROEC,axis=0)            # adjusted osmolyte
        Taxon_Osmo_Consti_Maint    = (Taxon_Osmo_Consti * Osmo_Maint_cost).sum(axis=1)        # maintenece
        Taxon_Osmo_Consti_Cost_N   = (Taxon_Osmo_Consti * Osmo_N_cost).sum(axis=1)            # N cost (no P)
        Taxon_Osmo_Consti_Cost_C   = Taxon_Osmo_Consti.sum(axis=1) + Taxon_Osmo_Consti_Maint  # total C consumption
        # Enzyme adjusted
        Taxon_Enzyme_Consti.loc[OECCN>0] = Taxon_Enzyme_Consti.loc[OECCN>0].mul(ARROEC,axis=0)                         # adjusted enzyme
        Taxon_Enzyme_Consti_Maint        = (Taxon_Enzyme_Consti.mul(self.Enz_Attrib['Maint_cost'],axis=1)).sum(axis=1) # maintinence
        Taxon_Enzyme_Consti_Cost_N       = (Taxon_Enzyme_Consti.mul(self.Enz_Attrib['N_cost'],    axis=1)).sum(axis=1) # N cost
        Taxon_Enzyme_Consti_Cost_P       = (Taxon_Enzyme_Consti.mul(self.Enz_Attrib['P_cost'],    axis=1)).sum(axis=1) # P cost
        Taxon_Enzyme_Consti_Cost_C       = Taxon_Enzyme_Consti.sum(axis=1) + Taxon_Enzyme_Consti_Maint                 # C cost (total)

        #---------------------------------------------------------------------#
        #.....Inducible processes.............................................#
        #---------------------------------------------------------------------#
        # Assimilation efficiency constrained by temperature
        Taxon_AE  = self.AE_ref + (self.temp[day] - (self.Tref - np.float(273))) * self.AE_temp  #scalar
        # Taxon growth respiration
        Taxon_Growth_Respiration = self.Taxon_Uptake_C * (np.float32(1) - Taxon_AE)
        
        #.................................................
        # Variable Acronyms:
        # OEICN : Osmo_Enzyme_Induci_Cost_N
        # OEIAN : Osmo_Enzyme_Induci_Avail_N
        # ARROEI: Avail_Req_ratio_osmo_enzyme_induci
        # MNAOEI: Min_N_Avail_Osmo_Enzyme_Induci
        #..................................................
        # Inducible Osmolyte production only when psi reaches below wp_fc
        Taxon_Osmo_Induci          = self.Induci_Osmo_C.mul(self.Taxon_Uptake_C * Taxon_AE,axis=0) * f_psi
        Taxon_Osmo_Induci_Cost_N   = (Taxon_Osmo_Induci * Osmo_N_cost).sum(axis=1)            # Total osmotic N cost of each taxon (.sum(axis=1))
        # Inducible enzyme production
        Taxon_Enzyme_Induci        = self.Induci_Enzyme_C.mul(self.Taxon_Uptake_C * Taxon_AE,axis=0)
        Taxon_Enzyme_Induci_Cost_N = (Taxon_Enzyme_Induci.mul(self.Enz_Attrib['N_cost'],axis=1)).sum(axis=1) # Total enzyme N cost of each taxon (.sum(axis=1))
        # Adjust production based on N availabe
        OEICN  = Taxon_Osmo_Induci_Cost_N + Taxon_Enzyme_Induci_Cost_N                         # Total N cost of osmolyte and enzymes
        OEIAN  = pd.Series(data=self.Taxon_Uptake_N, index=self.Microbes.index)                # N available
        MNAOEI = (pd.concat([OEICN[OEICN>0],OEIAN[OEICN>0]],axis=1)).min(axis=1,skipna=True)   # Get the minimum value by comparing N cost to N available
        ARROEI = (MNAOEI/OEICN[OEICN>0]).fillna(0)                                             # Ratio of Available to Required
        # Osmolyte adjusted: accompanying maintenence and N cost
        Taxon_Osmo_Induci[OEICN>0] = Taxon_Osmo_Induci.loc[OEICN>0].mul(ARROEI,axis=0)
        Taxon_Osmo_Induci_Maint    = (Taxon_Osmo_Induci * Osmo_Maint_cost).sum(axis=1) 
        Taxon_Osmo_Induci_Cost_N   = (Taxon_Osmo_Induci * Osmo_N_cost).sum(axis=1)
        Taxon_Osmo_Induci_Cost_C   = Taxon_Osmo_Induci.sum(axis=1) + Taxon_Osmo_Induci_Maint
        # Enzyme adjusted: Total enzyme carbon cost (+ CO2 loss), N cost, and P cost for each taxon
        Taxon_Enzyme_Induci[OEICN>0] = Taxon_Enzyme_Induci.loc[OEICN>0].mul(ARROEI,axis=0)
        Taxon_Enzyme_Induci_Maint    = (Taxon_Enzyme_Induci.mul(self.Enz_Attrib["Maint_cost"],axis=1)).sum(axis=1)
        Taxon_Enzyme_Induci_Cost_N   = (Taxon_Enzyme_Induci.mul(self.Enz_Attrib["N_cost"],    axis=1)).sum(axis=1)
        Taxon_Enzyme_Induci_Cost_P   = (Taxon_Enzyme_Induci.mul(self.Enz_Attrib["P_cost"],    axis=1)).sum(axis=1)
        Taxon_Enzyme_Induci_Cost_C   = Taxon_Enzyme_Induci.sum(axis=1) + Taxon_Enzyme_Induci_Maint
        # Derive C, N, & P deposited as biomass from Uptake; ensure no negative values
        Microbe_C_Gain = self.Taxon_Uptake_C - Taxon_Growth_Respiration   - Taxon_Enzyme_Induci_Cost_C - Taxon_Osmo_Induci_Cost_C
        Microbe_N_Gain = self.Taxon_Uptake_N - Taxon_Enzyme_Induci_Cost_N - Taxon_Osmo_Induci_Cost_N
        Microbe_P_Gain = self.Taxon_Uptake_P - Taxon_Enzyme_Induci_Cost_P

        #---------------------------------------------------------------------#
        #...............................Integration...........................#
        #---------------------------------------------------------------------#
        # Update Microbial pools with GAINS (from uptake) and LOSSES (from constitutive production)
        self.Microbes.loc[:,'C'] += Microbe_C_Gain - Taxon_Enzyme_Consti_Cost_C - Taxon_Osmo_Consti_Cost_C - Taxon_Transporter_Maint
        self.Microbes.loc[:,'N'] += Microbe_N_Gain - Taxon_Enzyme_Consti_Cost_N - Taxon_Osmo_Consti_Cost_N 
        self.Microbes.loc[:,'P'] += Microbe_P_Gain - Taxon_Enzyme_Consti_Cost_P
        self.Microbes[self.Microbes<0] = np.float32(0)  # avoid negative values
        # Growth_yield = Microbe_C_Gain - Taxon_Enzyme_Consti_Cost_C - Taxon_Osmo_Consti_Cost_C - Taxon_Transporter_Maint

        # Taxon-specific emergent CUE
        #CUE_taxon = Microbes['C'].copy() # create a dataframe and set all vals to 0
        #CUE_taxon[:] = 0
        #pos_uptake_index = self.Taxon_Uptake_C > 0
        #CUE_taxon[pos_uptake_index] = Microbe_C_Gain[pos_uptake_index]/self.Taxon_Uptake_C[pos_uptake_index]
        
        # System-level emergent CUE
        Taxon_Uptake_C_grid = self.Taxon_Uptake_C.sum()  # Total C Uptake
        if Taxon_Uptake_C_grid == 0:
            self.CUE_system = np.float32(0)
        else:
            self.CUE_system = Microbe_C_Gain.sum()/Taxon_Uptake_C_grid
        
        # Respiration from Constitutive + Inducible(NOTE: missing sum(MicLoss[,"C"]) in the Mortality below)
        self.Respiration = (Taxon_Transporter_Maint + Taxon_Growth_Respiration + Taxon_Osmo_Consti_Maint + Taxon_Osmo_Induci_Maint + Taxon_Enzyme_Consti_Maint + Taxon_Enzyme_Induci_Maint).sum(axis=0)
        
        # Derive Enzyme production
        Taxon_Enzyme_Production       = Taxon_Enzyme_Consti + Taxon_Enzyme_Induci  # gene-specific prod of enzyme of each taxon: (taxon*gridsize) * enzyme
        Taxon_Enzyme_Production.index = [np.arange(self.gridsize).repeat(self.n_taxa),Taxon_Enzyme_Production.index] # create a multi-index
        EP_df = Taxon_Enzyme_Production.groupby(level=0).sum() # enzyme-specific production in each grid cell
        Enzyme_Production = EP_df.stack().values # 1-D array
        # Derive Enzyme turnover
        Enzyme_Loss = self.Enzymes * Enzyme_Loss_Rate 

        # Update Enzyme pools by adding enzymes produced and substracting the 'dead' enzymes
        self.Enzymes += Enzyme_Production - Enzyme_Loss

        # Update Substrates pools with dead enzymes
        DeadEnz_df       = pd.concat([Enzyme_Loss,Enzyme_Loss.mul(self.Enz_Attrib['N_cost'].tolist()*self.gridsize,axis=0),Enzyme_Loss.mul(self.Enz_Attrib['P_cost'].tolist()*self.gridsize,axis=0)],axis=1)
        DeadEnz_df.index = [np.arange(self.gridsize).repeat(self.n_enzymes), DeadEnz_df.index] # create a multi-index
        DeadEnz_gridcell = DeadEnz_df.groupby(level=0).sum()  # total dead mass across taxa in each grid cell
        self.Substrates.loc[is_deadEnz] += DeadEnz_gridcell.values

    
    def mortality(self,day):       
        """
        Calculate microbial mortality and update stoichiometry of the alive and microbial pools.
        
        Kill microbes that are starving deterministically and microbes that are drought intolerant stochastically
        Also update Substrates with input from dead microbes, monomers(with leaching loss), and respiration
        """

        # Constants
        Leaching        = np.float32(0.1)  # Abiotic monomer loss rate
        Psi_slope_leach = np.float32(0.5)  # Mositure sensivity of abiotic monomer loss rate
        # Indices
        Mic_index  = self.Microbes.index
        is_DeadMic = self.Substrates.index == "DeadMic"
        is_NH4     = self.Monomers.index   == "NH4"
        is_PO4     = self.Monomers.index   == "PO4"
        
        # Reset the index to arabic numerals from taxa series 
        self.Microbes  = self.Microbes.reset_index(drop=True)
        MinRatios      = self.MinRatios.reset_index(drop=True)
        
        # Create a blank dataframe, Death, having the same structure as Microbes
        Death    = self.Microbes.copy(deep=True)
        Death[:] = np.float32(0)
        # Create a series, kill, holding boolean value of False
        kill = pd.Series([False]*self.n_taxa*self.gridsize)
        
        # Start of calcualtion of mortality first with THRESHOLD
        # Kill microbes deterministically based on threshold values: C_min: 0.086; N_min:0.012; P_min: 0.002
        starve_index = (self.Microbes["C"]>0) & ((self.Microbes["C"]<self.C_min)|(self.Microbes["N"]<self.N_min)|(self.Microbes["P"]<self.P_min))
        # Index the dead, put them in Death, and set them to 0 in Microbes 
        Death.loc[starve_index]         = self.Microbes[starve_index]
        self.Microbes.loc[starve_index] = np.float32(0)
        # Index the locations where microbial cells remain alive
        mic_index = self.Microbes["C"] > 0
        
        # Kill microbes stochastically based on mortality prob as a function of water potential and drought tolerance
        # call the function, MMP:microbe_mortality_psi() 
        r_death             = MMP(self.psi[day],self.wp_fc,self.basal_death_prob,self.death_rate,self.tolerance)
        kill.loc[mic_index] = r_death[mic_index] > np.random.uniform(0,1,sum(mic_index)).astype('float32')
        # Index the dead, put them in Death, and set them to 0 in Microbes 
        Death.loc[kill]         = self.Microbes[kill]
        self.Microbes.loc[kill] = np.float32(0)
        # Index locations where microbes remain alive
        mic_index = self.Microbes['C']>0

        # Calculate the total dead mass (threshold & drought) across taxa in each grid cell
        Death_gridcell = Death.groupby(Death.index//self.n_taxa).sum(axis=0)
        
        # Distinguish between conditions of complete death VS partial death
        # All cells die
        if sum(mic_index) == 0:
            
            #...Update Substrates pool by adding dead microbial biomass
            self.Substrates.loc[is_DeadMic] += Death_gridcell.values
        
        # Partly die and adjust stoichiometry of those remaining alive
        else:
            
            # Index only those taxa in Microbes that have below-minimum quotas: Mic_subset
            MicrobeRatios = self.Microbes[mic_index].divide(self.Microbes[mic_index].sum(axis=1),axis=0)
            mic_index_sub = (MicrobeRatios["C"]<MinRatios[mic_index]["C"])|(MicrobeRatios["N"]<MinRatios[mic_index]["N"])|(MicrobeRatios["P"]<MinRatios[mic_index]["P"])
            rat_index     = self.Microbes.index.map(mic_index_sub).fillna(False)
            # Derive the Microbes wanted
            Mic_subset    = self.Microbes[rat_index]
            StartMicrobes = Mic_subset.copy(deep=True)

            # Derive new ratios and Calculate difference between actual and min ratios  
            MicrobeRatios = Mic_subset.divide(Mic_subset.sum(axis=1),axis=0)
            MinRat        = MinRatios[rat_index]  
            Ratio_dif     = MicrobeRatios - MinRat
            # Create a df recording the ratio differences < 0
            Ratio_dif_0 = Ratio_dif.copy(deep=True)
            Ratio_dif_0[Ratio_dif>0] = np.float32(0)  
            # Create a df recording the ratio differences > 0
            Excess = Ratio_dif.copy(deep=True)
            Excess[Ratio_dif<0] = np.float32(0)

            # Determine the limiting nutrient that will be conserved
            Limiting = (-Ratio_dif/MinRat).idxmax(axis=1) # Series of index of the first occurrence of maximum in each row
            # Set all deficient ratios to their minima
            MicrobeRatios[Ratio_dif<0] = MinRat[Ratio_dif<0]
            # Reduce the mass fractions for non-deficient elements in proportion to the distance from the minimum
            # ....Partition the total deficit to the excess element(s) in proportion to their distances from their minima
            MicrobeRatios[Ratio_dif>0] += Excess.mul((Ratio_dif_0.sum(axis=1)/Excess.sum(axis=1)),axis=0)[Ratio_dif>0]
            
            # Construct hypothetical nutrient quotas for each possible minimum nutrient
            MC  = Mic_subset["C"]
            MN  = Mic_subset["N"]
            MP  = Mic_subset["P"]
            MRC = MicrobeRatios["C"]
            MRN = MicrobeRatios["N"]
            MRP = MicrobeRatios["P"]

            new_C = pd.concat([MC, MN*MRC/MRN, MP*MRC/MRP],axis=1)
            new_C = new_C.fillna(0)
            new_C[np.isinf(new_C)] = np.float32(0)
            new_C.columns = ['C','N','P']

            new_N = pd.concat([MC*MRN/MRC, MN, MP*MRN/MRP],axis=1)
            new_N = new_N.fillna(0)
            new_N[np.isinf(new_N)] = np.float32(0)
            new_N.columns = ['C','N','P']
            
            new_P = pd.concat([MC*MRP/MRC, MN*MRP/MRN, MP],axis=1)
            new_P = new_P.fillna(0)
            new_P[np.isinf(new_P)] = np.float32(0)
            new_P.columns = ['C','N','P']
            
            # Insert the appropriate set of nutrient quotas scaled to the minimum nutrient
            C = [new_C.loc[i,Limiting[i]] for i in Limiting.index] #list
            N = [new_N.loc[i,Limiting[i]] for i in Limiting.index] #list
            P = [new_P.loc[i,Limiting[i]] for i in Limiting.index] #list
            
            # Update Microbes
            self.Microbes.loc[rat_index] = np.vstack((C,N,P)).transpose()

            # Sum up the element losses from biomass across whole grid and calculate average loss
            MicLoss = StartMicrobes - self.Microbes[rat_index]
            # Update total respiration by adding ...
            self.Respiration += sum(MicLoss['C'])
            # Update monomer pools 
            self.Monomers.loc[is_NH4,"N"] += sum(MicLoss["N"])/self.gridsize
            self.Monomers.loc[is_PO4,"P"] += sum(MicLoss["P"])/self.gridsize
            
            # Update Substrates pool by adding dead microbial biomass            
            self.Substrates.loc[is_DeadMic] += Death_gridcell.values
        # End of if else clause
        
        # Leaching of monomers
        Leaching_N = self.Monomers.loc[is_NH4,"N"] * Leaching * np.exp(Psi_slope_leach * (self.psi[day]-self.wp_fc))
        Leaching_P = self.Monomers.loc[is_PO4,"P"] * Leaching * np.exp(Psi_slope_leach * (self.psi[day]-self.wp_fc))
        # Update Monomers
        self.Monomers.loc[is_NH4,"N"] -= Leaching_N
        self.Monomers.loc[is_PO4,"P"] -= Leaching_P
        
        # Restore the index to taxa series
        self.Microbes.index = Mic_index

        # Update the death toll of cells
        self.Kill = kill.sum().astype('uint32')

    
    def reproduction(self,day):           
        """
        Calculate reproduction and dispersal, and update microbial composition/distrituion on the spatial grid.
 
        Parameters:
            fb         : index of fungal taxa
            max_size_b : threshold of cell division
            max_size_f : threshold of cell division
            x,y        : x,y dimension of grid
            dist       : maximum dispersal distance: 1 cell
            direct     : dispersal direction: 0.95    
        """
        
        # Microbes' index
        Mic_index = self.Microbes.index
        # Set up the Colonization dataframe: taxon * 3(C,N,&P)
        Colonization    = self.Microbes.copy(deep=True)
        Colonization    = Colonization.reset_index(drop=True)
        Colonization[:] = np.float32(0)
        
        
        #STEP 1: Fungal translocation by calculating average biomass within fungal taxa 
        # Count the fungal taxa before cell division
        # Add one or two fungi to the count series based on size
        Fungi_df = pd.Series(data=[0]*self.n_taxa*self.gridsize, index=Mic_index, name='Count', dtype='int8')
        Fungi_df.loc[(self.fb==1)&(self.Microbes['C']>0)]               = np.int8(1)
        Fungi_df.loc[(self.fb==1)&(self.Microbes['C']>self.max_size_f)] = np.int8(2)
        Fungi_count   = Fungi_df.groupby(level=0,sort=False).sum()
        # Derive average biomass of fungal taxa
        Microbes_grid              = self.Microbes.groupby(level=0,sort=False).sum()
        Mean_fungi                 = Microbes_grid.divide(Fungi_count,axis=0)
        Mean_fungi[Fungi_count==0] = np.float32(0)
        # Expand the fungal average across the grid
        eMF = expand(Mean_fungi,self.gridsize)
        
        
        #STEP 2: Cell division & translocate nutrients
        MicrobesBeforeDivision = self.Microbes.copy(deep=True)
        #bacterial cell division
        bac_index = (self.fb==0) & (self.Microbes['C']>self.max_size_b)
        self.Microbes[bac_index] = self.Microbes[bac_index]/2
        #fungal cell division
        fun_index = (self.fb==1) & (self.Microbes['C']>self.max_size_f)
        self.Microbes[fun_index] = self.Microbes[fun_index]/2
        # Put daughter cells into a seperate dataframe, Reprod
        Reprod = MicrobesBeforeDivision - self.Microbes
        # Translocate nutrients within fungal taxa after reproduction
        self.Microbes[(self.fb==1)&(self.Microbes['C']>0)] = eMF[(self.fb==1)&(self.Microbes['C']>0)]
        # Index the daughter cells of fungi vs bacteria
        daughters_b = (Reprod['C']>0) & (self.fb==0)
        daughters_f = (Reprod['C']>0) & (self.fb==1)
        # set all fungi equal to their grid averages for translocation before colonization
        Reprod[daughters_f] = eMF[daughters_f]
        

        #STEP 3: dispersal calculation
        num_b   = sum(daughters_b)
        num_f   = sum(daughters_f)
        shift_x = pd.Series(data=[0] * self.gridsize*self.n_taxa, index=Mic_index, dtype='int8')
        shift_y = pd.Series(data=[0] * self.gridsize*self.n_taxa, index=Mic_index, dtype='int8')
        # Bacterial dispersal movements in X & Y direction  
        shift_x.loc[daughters_b] = np.random.choice([i for i in range(-self.dist, self.dist+1)],num_b,replace=True).astype('int8')
        shift_y.loc[daughters_b] = np.random.choice([i for i in range(-self.dist, self.dist+1)],num_b,replace=True).astype('int8')
        # Fungi always move positively in x direction, and in y direction constrained to one box away determined by probability "direct"                
        shift_x.loc[daughters_f] = np.int8(1)
        shift_y.loc[daughters_f] = np.random.choice([-1,0,1], num_f, replace=True, p=[0.5*(1-self.direct),self.direct,0.5*(1-self.direct)]).astype('int8')
        # Calculate x,y coordinates of dispersal destinations (% remainder of x/x)
        new_x           = (list(np.repeat(range(1,self.x+1),self.n_taxa)) * self.y + shift_x + self.x) % self.x
        new_y           = (list(np.repeat(range(1,self.y+1),self.n_taxa*self.x))   + shift_y + self.y) % self.y
        new_x[new_x==0] = self.x  # Substitute coordinates when there is no shift 
        new_y[new_y==0] = self.y  # Substitute coordinates when there is no shift
        # Convert x,y coordinates to a Series of destination locations
        index_series = (self.n_taxa * ((new_y-1)*self.x + (new_x-1))) + list(range(1,self.n_taxa+1)) * self.gridsize - 1
        

        #Step 4: colonization of dispersed microbes
        # Transfer reproduced cells to new locations and sum when two or more of the same taxa go to same location
        Colonization.iloc[index_series[daughters_b],] = Reprod[daughters_b].values
        Colonization.iloc[index_series[daughters_f],] = Reprod[daughters_f].values
        # Colonization of dispersing microbes
        self.Microbes += Colonization.values

    
    def repopulation(self,output,day,mic_reinit):
        """
        Reinitialize a microbial community on the spatial grid in a new pulse.
        
        Meanwhile, start with new subsrates, monomers, and enzymes on the grid that are initialized from the very beginning

        Parameters:
            output:     an instance of the Output class, from which the var, MicrobesSeries_repop,
                        referring to taxon-specific total mass over the grid is retrieved
            pulse:      the pulse index
            day:        the day index
            mic_reinit: 0/1; 1 means reinitialization
        Returns:
            update Substrates, Monomers, and Microbes, as well as Enzymes
        """

        # reinitialize substrates, monomers, and enzymes in a new pulse
        self.Substrates = self.Substrates_init.copy(deep=True)
        self.Monomers   = self.Monomers_init.copy(deep=True)
        self.Enzymes    = self.Enzymes_init.copy(deep=True)

        # reinitialize microbial community in a new pulse if True
        if mic_reinit == True:
            
            self.Microbes = self.Microbes_init.copy(deep=True) #NOTE copy()
            #fb = self.fb[0:self.n_taxa]
            #max_size_b = self.max_size_b
            #max_size_f = self.max_size_f
            
            # cumulative abundance; note the column index
            # option 1: mass-based.
            #cum_abundance = output.MicrobesSeries_repop.iloc[:,(day+2-self.cycle):(day+2)].sum(axis=1)
            # option 2: abundance-based
            cum_abundance = output.Taxon_count_repop.iloc[:,(day+2-self.cycle):(day+2)].sum(axis=1)
            
            # account for different cell mass sizes of bacteria and fungi
            #if sum(fb==1) == 0: # no fungal taxa
            #    frequencies = cum_abundance/cum_abundance.sum()
            #else:
            #    cum_abundance.loc[fb==1] = cum_abundance[fb==1]*max_size_b/max_size_f
            #    frequencies = cum_abundance/cum_abundance.sum()
            
            # Switched to taxon abundance-based, so no more adjustments
            frequencies = cum_abundance/cum_abundance.sum()
            frequencies = frequencies.fillna(0)
            probs       = pd.concat([frequencies,1-frequencies],axis=1,sort=False)
            # Randomly assign microbes to each grid box based on prior densities
            choose_taxa = np.array([0]* self.gridsize * self.n_taxa).reshape(self.n_taxa,self.gridsize)
            for i in range(self.n_taxa):
                # Alternative 1
                choose_taxa[i,:] = np.random.choice([1,0],self.gridsize,replace=True,p=probs.iloc[i,:])
                        
            # Note order='F'
            self.Microbes.loc[np.ravel(choose_taxa,order='F')==0] = np.float32(0)