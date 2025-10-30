"""
DEMENTpy Library Interface

This module provides a clean API for using DEMENTpy as a library
rather than as a standalone command-line application.

Key features:
- Initialize Grid from dictionaries (no file I/O required)
- Simplified parameter interface
- State save/restore
- Output extraction for coupling with other models
"""

import numpy as np
import pandas as pd
from typing import Dict, Any, Optional, Tuple
import os
import sys

# Add current directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from grid import Grid
from initialization import initialize_data
from utility import expand
import tempfile
import shutil


class DEMENTLibrary:
    """
    Library interface for DEMENTpy Grid class.

    Provides methods to:
    1. Create Grid instances without file I/O
    2. Run daily timesteps
    3. Extract outputs for coupling
    4. Save/restore state
    """


    @staticmethod
    def create_default_runtime(
        end_time: int = 365,
        gridsize: int = 100,
        n_taxa: int = 5,
        n_substrates: int = 3,
        n_enzymes: int = 3,
        x: int = 10,
        y: int = 10,
        dist: int = 1,
        direct: float = 0.95
    ) -> pd.DataFrame:
        """
        Create a default runtime configuration DataFrame.

        This mimics the format expected by Grid.__init__() but allows
        creation without reading from file.

        Args:
            end_time: Number of daily timesteps (default: 365)
            gridsize: Total number of grid cells (default: 100 = 10x10)
            n_taxa: Number of microbial taxa
            n_substrates: Number of substrate types
            n_enzymes: Number of enzyme types
            x: X dimension of grid (default: 10)
            y: Y dimension of grid (default: 10)
            dist: Maximum dispersal distance in cells (default: 1)
            direct: Dispersal directionality (default: 0.95)

        Returns:
            DataFrame in the format expected by Grid
        """
        # Validate that gridsize matches x*y if provided
        if gridsize != x * y:
            raise ValueError(f"gridsize ({gridsize}) must equal x*y ({x}*{y}={x*y})")

        runtime = pd.DataFrame({
            0: [
                'pulse', 'end_time', 'x', 'y', 'gridsize', 'n_taxa', 'fb',
                'taxa_per_box', 'n_enzymes', 'n_substrates', 'n_uptake',
                'n_osmolytes', 'NormalizeUptake', 'NormalizeProd',
                'Enzymes_per_sub', 'Avg_extra_req_enz', 'Uptake_per_monomer',
                'Monomer_Substrate_Ratio', 'Enz_min', 'Enz_max',
                'Init_NH4', 'Init_PO4', 'Input_NH4', 'Input_PO4',
                'direct', 'dist', 'interval'
            ],
            1: [
                10,         # pulse
                end_time,   # end_time
                x,          # x dimension
                y,          # y dimension
                gridsize,   # gridsize
                n_taxa,     # n_taxa
                0.0,        # fb (fungal fraction)
                0.01,       # taxa_per_box
                n_enzymes,  # n_enzymes
                n_substrates,  # n_substrates
                n_substrates + 2,  # n_uptake (monomers)
                10,         # n_osmolytes
                0,          # NormalizeUptake
                0,          # NormalizeProd
                1,          # Enzymes_per_sub
                0,          # Avg_extra_req_enz
                1,          # Uptake_per_monomer
                0,          # Monomer_Substrate_Ratio
                0,          # Enz_min
                0,          # Enz_max
                0,          # Init_NH4
                0,          # Init_PO4
                0,          # Input_NH4
                0,          # Input_PO4
                direct,     # direct dispersal directionality
                dist,       # dist dispersal distance
                1           # interval
            ]
        })
        runtime = runtime.set_index(0)
        return runtime

    @staticmethod
    def create_minimal_initialization(
        gridsize: int = 10,
        n_taxa: int = 5,
        n_substrates: int = 3,
        substrate_c: float = 500.0,  # g C/m²
        substrate_n: float = 20.0,   # g N/m²
        substrate_p: float = 5.0,    # g P/m²
        microbial_c: float = 50.0,   # g C/m²
        temperature: float = 15.0,   # °C
        moisture_psi: float = -0.5,  # MPa
        fb: float = 0.1,             # Fraction of fungal taxa
        parameters_dict: Dict = None # Optional override parameters
    ) -> Dict[str, Any]:
        """
        Create a complete initialization dictionary for Grid with all 60+ required parameters.

        This provides default values for all required data structures,
        allowing Grid to be instantiated without file I/O.

        Args:
            gridsize: Grid dimension
            n_taxa: Number of microbial taxa
            n_substrates: Number of substrate types
            substrate_c: Initial substrate carbon (g/m²)
            substrate_n: Initial substrate nitrogen (g/m²)
            substrate_p: Initial substrate phosphorus (g/m²)
            microbial_c: Initial microbial biomass C (g/m²)
            temperature: Initial temperature (°C)
            moisture_psi: Water potential (MPa)
            fb: Fraction of fungal taxa (0-1)
            parameters_dict: Optional dictionary to override default parameter values

        Returns:
            Dictionary with all keys expected by Grid.__init__()
        """
        # Load default parameters (from typical parameters.csv)
        if parameters_dict is None:
            parameters_dict = {}

        # Default parameter values from DEMENTpy/input/parameters.csv
        params = {
            'max_size_b': 2.0,
            'max_size_f': 50.0,
            'Cfrac_b': 0.825,
            'Nfrac_b': 0.16,
            'Pfrac_b': 0.015,
            'Cfrac_f': 0.9,
            'Nfrac_f': 0.09,
            'Pfrac_f': 0.01,
            'C_min': 0.086,
            'N_min': 0.012,
            'P_min': 0.002,
            'Uptake_Maint_cost': 0.01,
            'Enz_Prod_min': 0.00001,
            'Enz_Prod_max': 0.0001,
            'Constit_Prod_min': 0.00001,
            'Constit_Prod_max': 0.0001,
            'Osmo_Consti_Prod_min': 0.0000001,
            'Osmo_Consti_Prod_max': 0.000001,
            'Osmo_Induci_Prod_min': 0.01,
            'Osmo_Induci_Prod_max': 0.1,
            'CUE_ref': 0.5,
            'CUE_temp': -0.005,
            'death_rate_bac': 0.001,
            'death_rate_fun': 0.001,
            'wp_fc': -2.0,
            'wp_th': -6.0,
            'alpha': 0.01,
            'Vmax0_min': 5.0,
            'Vmax0_max': 50.0,
            'Uptake_Vmax0_min': 1.0,
            'Uptake_Vmax0_max': 10.0,
            'Uptake_Ea_min': 35.0,
            'Uptake_Ea_max': 35.0,
            'Km_min': 0.01,
            'Uptake_Km_min': 0.001,
        }
        # Override with user-provided parameters
        params.update(parameters_dict)

        n_cells = gridsize * gridsize
        n_enzymes = n_substrates  # One enzyme per substrate type
        n_monomers = n_substrates + 2  # C, N, P

        # ===================================================================
        # 1. SUBSTRATES
        # ===================================================================
        # Substrates must include 'DeadMic' and 'DeadEnz' for recycling dead biomass and enzymes
        special_substrates = ['DeadMic', 'DeadEnz']
        organic_substrate_names = [f'Sub{i+1}' for i in range(n_substrates)]
        substrate_names = special_substrates + organic_substrate_names
        n_substrates_total = len(substrate_names)

        # Initial substrate pool - DataFrame with columns ['C', 'N', 'P'] and substrate names as index
        # DeadMic and DeadEnz start at 0; organic substrates get the input carbon
        c_per_substrate = substrate_c / n_substrates if n_substrates > 0 else 0
        n_per_substrate = substrate_n / n_substrates if substrate_n > 0 else 0
        p_per_substrate = substrate_p / n_substrates if substrate_p > 0 else 0

        substrate_c_values = [0, 0] + [c_per_substrate] * n_substrates
        substrate_n_values = [0, 0] + [n_per_substrate] * n_substrates
        substrate_p_values = [0, 0] + [p_per_substrate] * n_substrates

        substrates_base = pd.DataFrame({
            'C': substrate_c_values,
            'N': substrate_n_values,
            'P': substrate_p_values
        }, index=substrate_names)

        # Expand across grid: (n_substrates_total, 3) → (n_substrates_total*gridsize, 3)
        substrates_df = expand(substrates_base, gridsize)

        # Substrate inputs (daily) - only organic substrates get inputs, not DeadMic/DeadEnz
        sub_input_df = pd.DataFrame(
            0.1,  # Small daily input (g/m²/day)
            index=range(365),
            columns=organic_substrate_names
        )

        # Substrate degradation required enzymes (ReqEnz matrix)
        # Each organic substrate requires one specific enzyme (diagonal for organic subs)
        # DeadMic and DeadEnz don't require enzymes for degradation
        # Shape: (n_substrates_total, n_enzymes) → (n_substrates_total*gridsize, n_enzymes)
        req_enz_matrix = np.zeros((n_substrates_total, n_enzymes), dtype=int)
        # Only organic substrates need enzymes (diagonal pattern)
        for i in range(min(n_substrates, n_enzymes)):
            req_enz_matrix[2 + i, i] = 1  # offset by 2 for DeadMic and DeadEnz

        req_enz_base = pd.DataFrame(
            req_enz_matrix,
            index=['set1'] * n_substrates_total,
            columns=[f'Enz{i+1}' for i in range(n_enzymes)]
        )
        req_enz = expand(req_enz_base, gridsize)

        # Substrates-produced monomers (MonomersProduced)
        # All substrates produce C, N, P monomers when degraded
        monomers_produced = pd.DataFrame({
            'C': [0.5] * n_substrates_total,  # 50% of C goes to monomers
            'N': [0.05] * n_substrates_total, # 5% of mass is N
            'P': [0.005] * n_substrates_total  # 0.5% of mass is P
        }, index=substrate_names)

        # ===================================================================
        # 2. ENZYMES
        # ===================================================================
        enzyme_names = [f'Enz{i+1}' for i in range(n_enzymes)]

        # Enzymes initialization - Series with enzyme names as index, expanded across grid
        # Shape: (n_enzymes,) → (n_enzymes*gridsize,)
        # Index will have enzyme names repeated gridsize times
        enzymes_base = pd.Series(
            0.0,
            index=enzyme_names,
            dtype=np.float32
        )
        enzymes_series = expand(enzymes_base, gridsize)

        # Enzyme Vmax0 - DataFrame with enzyme names as index, ALL substrate names as columns
        # Shape: (n_enzymes, n_substrates_total) → (n_enzymes*gridsize, n_substrates_total)
        # Index: enzyme names repeated gridsize times
        vmax_base = pd.DataFrame(
            np.full((n_enzymes, n_substrates_total),
                   (params['Vmax0_min'] + params['Vmax0_max']) / 2,
                   dtype=np.float32),
            index=enzyme_names,
            columns=substrate_names
        )
        vmax_values = expand(vmax_base, gridsize)

        # Enzyme Km0 - DataFrame with ALL substrate names as index, enzyme names as columns
        # Shape: (n_substrates_total, n_enzymes) → (n_substrates_total*gridsize, n_enzymes)
        # Index: substrate names repeated gridsize times
        km_base = pd.DataFrame(
            np.full((n_substrates_total, n_enzymes),
                   params['Km_min'],
                   dtype=np.float32),
            index=substrate_names,
            columns=enzyme_names
        )
        km_values = expand(km_base, gridsize)

        # Enzyme Ea - Same structure as Vmax0
        # Shape: (n_enzymes, n_substrates_total) → (n_enzymes*gridsize, n_substrates_total)
        ea_base = pd.DataFrame(
            np.full((n_enzymes, n_substrates_total), 37.0, dtype=np.float32),
            index=enzyme_names,
            columns=substrate_names
        )
        ea_values = expand(ea_base, gridsize)

        # Enzyme attributes (constitutive vs inducible)
        enz_attrib_data = []
        for i in range(1, n_enzymes + 1):
            enz_attrib_data.append({
                'Enzyme': f'Enz{i}',
                'Induci': 1,  # Inducible
                'Consti': 0   # Not constitutive (for simplicity)
            })
        enz_attrib_df = pd.DataFrame(enz_attrib_data).set_index('Enzyme')

        # ===================================================================
        # 3. MONOMERS
        # ===================================================================
        # Total monomers = inorganic (NH4, PO4) + organic monomers from substrates
        # Organic monomers are named after substrates (Mon1, Mon2, ...) or use substrate names
        inorganic_monomers = ['NH4', 'PO4']
        organic_monomer_names = [f'Mon{i+1}' for i in range(n_substrates)]
        all_monomer_names = inorganic_monomers + organic_monomer_names
        n_monomers_actual = len(all_monomer_names)

        # Monomers initialization - DataFrame with columns ['C', 'N', 'P'] and monomer names as index
        # Shape: (n_monomers_total, 3) → (n_monomers_total*gridsize, 3)
        monomers_base = pd.DataFrame({
            'C': [0.01] * n_monomers_actual,  # Small initial pools
            'N': [0.01] * n_monomers_actual,
            'P': [0.01] * n_monomers_actual
        }, index=all_monomer_names)
        monomers_df = expand(monomers_base, gridsize)

        # Monomer inputs (daily) - only for inorganic + general monomers
        mon_input_df = pd.DataFrame(
            0.01,  # Small daily input
            index=range(365),
            columns=inorganic_monomers + ['C', 'N', 'P']  # NH4, PO4, C, N, P inputs
        )

        # Monomer stoichiometry ratios - start from Monomers structure and set to zeros
        # This mimics what the real initialization does
        monomer_ratios_df = monomers_df.copy()
        monomer_ratios_df[:] = 0.0
        # Set NH4 to have N=1, PO4 to have P=1
        monomer_ratios_df.loc[monomer_ratios_df.index == 'NH4', 'N'] = 1.0
        monomer_ratios_df.loc[monomer_ratios_df.index == 'PO4', 'P'] = 1.0
        # Organic monomers will get updated during uptake() from SubstrateRatios

        # Monomer uptake required enzymes (identity matrix - each monomer needs one transporter)
        # Shape: (n_monomers, n_monomers) → (n_monomers*gridsize, n_monomers)
        uptake_req_enz_base = pd.DataFrame(
            np.eye(n_monomers_actual, dtype=int),
            index=all_monomer_names,
            columns=all_monomer_names
        )
        uptake_req_enz = expand(uptake_req_enz_base, gridsize)

        # Uptake transporter parameters - DataFrames with monomer names as index, transporter names as columns
        # Shape: (n_monomers, n_monomers) → (n_monomers*gridsize, n_monomers)
        uptake_transporter_names = [f'Upt{i+1}' for i in range(n_monomers_actual)]

        uptake_ea_base = pd.DataFrame(
            np.full((n_monomers_actual, n_monomers_actual), params['Uptake_Ea_min'], dtype=np.float32),
            index=all_monomer_names,
            columns=uptake_transporter_names
        )
        uptake_ea = expand(uptake_ea_base, gridsize)

        uptake_vmax_base = pd.DataFrame(
            np.full((n_monomers_actual, n_monomers_actual),
                   (params['Uptake_Vmax0_min'] + params['Uptake_Vmax0_max']) / 2,
                   dtype=np.float32),
            index=all_monomer_names,
            columns=uptake_transporter_names
        )
        uptake_vmax = expand(uptake_vmax_base, gridsize)

        uptake_km_base = pd.DataFrame(
            np.full((n_monomers_actual, n_monomers_actual), params['Uptake_Km_min'], dtype=np.float32),
            index=all_monomer_names,
            columns=uptake_transporter_names
        )
        uptake_km = expand(uptake_km_base, gridsize)

        # ===================================================================
        # 4. MICROBES
        # ===================================================================
        # Determine number of bacterial vs fungal taxa
        n_fungi = int(n_taxa * fb)
        n_bacteria = n_taxa - n_fungi
        fungi_indices = list(range(n_bacteria + 1, n_taxa + 1)) if n_fungi > 0 else []

        # Fungal index array (1 if taxon is fungal, 0 if bacterial)
        fb_array = np.zeros(n_taxa, dtype=np.float32)
        if n_fungi > 0:
            fb_array[-n_fungi:] = 1.0

        # Create initial microbial community
        microbes_data = []
        cells_per_taxon = max(1, n_cells // n_taxa)

        for taxon in range(1, n_taxa + 1):
            is_fungal = (taxon in fungi_indices)
            c_frac = params['Cfrac_f'] if is_fungal else params['Cfrac_b']
            n_frac = params['Nfrac_f'] if is_fungal else params['Nfrac_b']
            p_frac = params['Pfrac_f'] if is_fungal else params['Pfrac_b']

            for i in range(cells_per_taxon):
                cell_idx = (taxon - 1) * cells_per_taxon + i
                if cell_idx >= n_cells:
                    break

                x = cell_idx % gridsize
                y = cell_idx // gridsize

                c_biomass = microbial_c / n_cells
                total_biomass = c_biomass / c_frac

                microbes_data.append({
                    'x': x,
                    'y': y,
                    'Taxon': taxon,
                    'C': c_biomass,
                    'N': total_biomass * n_frac,
                    'P': total_biomass * p_frac,
                    **{f'Enz{i}': 0.5 for i in range(1, n_enzymes+1)}  # Initial enzyme gene expression
                })

        microbes_df = pd.DataFrame(microbes_data)

        # Microbial pool (before placement on grid) - same as microbes_df
        microbes_pp = microbes_df.copy()

        # Bacterial and fungal densities (cells per grid cell)
        bac_density = len(microbes_df[microbes_df['Taxon'] <= n_bacteria]) / n_cells
        fun_density = len(microbes_df[microbes_df['Taxon'] > n_bacteria]) / n_cells if n_fungi > 0 else 0.0

        # ===================================================================
        # 5. MICROBIAL TRAITS AND GENES
        # ===================================================================
        # All gene DataFrames use taxon names as index
        taxon_names = [f'Tax{i}' for i in range(1, n_taxa + 1)]

        # Enzyme genes (which taxa produce which enzymes)
        # For simplicity: all taxa can produce all enzymes (with some random variation)
        enzyme_genes = pd.DataFrame(
            np.random.randint(0, 2, size=(n_taxa, n_enzymes)),
            index=taxon_names,
            columns=enzyme_names
        ).astype(np.float32)

        # Osmolyte genes (for drought tolerance)
        osmolyte_genes = pd.DataFrame(
            np.random.uniform(0, 1, size=(n_taxa, 2)),  # Constitutive and inducible
            index=taxon_names,
            columns=['OsmoConsti', 'OsmoInduci']
        ).astype(np.float32)

        # Uptake genes (for monomer uptake) - all taxa can uptake all monomers
        uptake_genes = pd.DataFrame(
            np.ones((n_taxa, n_monomers_actual)),  # All taxa can uptake all monomers
            index=taxon_names,
            columns=all_monomer_names
        ).astype(np.float32)

        # Gene costs and production rates - all need taxon names as index and expand across grid
        # Trait arrays (not expanded)
        uptake_genes_trait = np.full(n_monomers_actual, params['Uptake_Maint_cost'], dtype=np.float32)
        osmo_prod_consti_trait = np.full(2, (params['Osmo_Consti_Prod_min'] + params['Osmo_Consti_Prod_max']) / 2, dtype=np.float32)
        osmo_prod_induci_trait = np.full(2, (params['Osmo_Induci_Prod_min'] + params['Osmo_Induci_Prod_max']) / 2, dtype=np.float32)
        enz_prod_consti_trait = np.full(n_enzymes, (params['Constit_Prod_min'] + params['Constit_Prod_max']) / 2, dtype=np.float32)
        enz_prod_induci_trait = np.full(n_enzymes, (params['Enz_Prod_min'] + params['Enz_Prod_max']) / 2, dtype=np.float32)

        # Taxon-specific cost distributions (need to be expanded)
        # UptakeGenesCost: (n_taxa, n_monomers) → (n_taxa*gridsize, n_monomers)
        uptake_genes_cost_base = pd.DataFrame(
            np.full((n_taxa, n_monomers_actual), params['Uptake_Maint_cost'], dtype=np.float32),
            index=taxon_names,
            columns=all_monomer_names
        )
        uptake_genes_cost = expand(uptake_genes_cost_base, gridsize)

        # OsmoProdConsti: (n_taxa, 2) → (n_taxa*gridsize, 2)
        osmo_prod_consti_base = pd.DataFrame(
            np.tile(osmo_prod_consti_trait, (n_taxa, 1)),
            index=taxon_names,
            columns=['OsmoConsti1', 'OsmoConsti2']
        ).astype(np.float32)
        osmo_prod_consti = expand(osmo_prod_consti_base, gridsize)

        # OsmoProdInduci: (n_taxa, 2) → (n_taxa*gridsize, 2)
        osmo_prod_induci_base = pd.DataFrame(
            np.tile(osmo_prod_induci_trait, (n_taxa, 1)),
            index=taxon_names,
            columns=['OsmoInduci1', 'OsmoInduci2']
        ).astype(np.float32)
        osmo_prod_induci = expand(osmo_prod_induci_base, gridsize)

        # EnzProdConstit: (n_taxa, n_enzymes) → (n_taxa*gridsize, n_enzymes)
        enz_prod_constit_base = pd.DataFrame(
            np.tile(enz_prod_consti_trait, (n_taxa, 1)),
            index=taxon_names,
            columns=enzyme_names
        ).astype(np.float32)
        enz_prod_constit = expand(enz_prod_constit_base, gridsize)

        # EnzProdInduce: (n_taxa, n_enzymes) → (n_taxa*gridsize, n_enzymes)
        enz_prod_induce_base = pd.DataFrame(
            np.tile(enz_prod_induci_trait, (n_taxa, 1)),
            index=taxon_names,
            columns=enzyme_names
        ).astype(np.float32)
        enz_prod_induce = expand(enz_prod_induce_base, gridsize)

        # TaxDroughtTol: (n_taxa, 1) → (n_taxa*gridsize, 1)
        tax_drought_tol_base = pd.DataFrame(
            np.random.uniform(params['wp_th'], params['wp_fc'], size=n_taxa),
            index=taxon_names,
            columns=['DroughtTol']
        ).astype(np.float32)
        tax_drought_tol = expand(tax_drought_tol_base, gridsize)

        # ===================================================================
        # 6. MORTALITY PARAMETERS
        # ===================================================================
        # Basal death probability and death rate
        basal_death_prob = np.full(n_taxa, params['death_rate_bac'], dtype=np.float32)
        basal_death_prob[fungi_indices if fungi_indices else []] = params['death_rate_fun']

        death_rate = params['death_rate_bac']  # Scalar

        # ===================================================================
        # 7. CLIMATE AND ENVIRONMENTAL PARAMETERS
        # ===================================================================
        # Temperature and moisture series (365 days)
        temp_series = np.full(365, temperature, dtype=np.float32)
        psi_series = np.full(365, moisture_psi, dtype=np.float32)

        # ===================================================================
        # 8. OTHER PARAMETERS
        # ===================================================================
        # Assimilation efficiency
        ae_ref = params['CUE_ref']
        ae_temp = params['CUE_temp']

        # Minimum cell quotas
        min_ratios_base = pd.DataFrame({
            'C:N': [params['Cfrac_b'] / params['Nfrac_b']],
            'C:P': [params['Cfrac_b'] / params['Pfrac_b']],
            'N:P': [params['Nfrac_b'] / params['Pfrac_b']]
        })
        min_ratios = expand(min_ratios_base, gridsize)

        # Expand all the gene-related parameters that haven't been expanded yet
        uptake_genes_expanded = expand(uptake_genes, gridsize)
        osmolyte_genes_expanded = expand(osmolyte_genes, gridsize)
        enzyme_genes_expanded = expand(enzyme_genes, gridsize)

        uptake_genes_trait_expanded = expand(pd.Series(uptake_genes_trait), gridsize)
        osmo_prod_consti_trait_expanded = expand(pd.Series(osmo_prod_consti_trait), gridsize)
        osmo_prod_induci_trait_expanded = expand(pd.Series(osmo_prod_induci_trait), gridsize)
        enz_prod_consti_trait_expanded = expand(pd.Series(enz_prod_consti_trait), gridsize)
        enz_prod_induci_trait_expanded = expand(pd.Series(enz_prod_induci_trait), gridsize)

        # Note: uptake_genes_cost, osmo_prod_consti, osmo_prod_induci, enz_prod_constit,
        # enz_prod_induce, tax_drought_tol are already expanded above

        # ===================================================================
        # 9. ASSEMBLE THE COMPLETE DATA DICTIONARY
        # ===================================================================
        data_init = {
            # Substrates
            'Substrates': substrates_df,
            'SubInput': sub_input_df,
            'ReqEnz': req_enz,
            'MonomersProduced': monomers_produced,

            # Monomers
            'Monomers': monomers_df,
            'Monomer_ratio': monomer_ratios_df,
            'MonInput': mon_input_df,
            'Uptake_ReqEnz': uptake_req_enz,

            # Enzymes
            'Enzymes': enzymes_series,
            'Km0': km_values,
            'Uptake_Km0': uptake_km,
            'Uptake_Ea': uptake_ea,
            'Uptake_Vmax0': uptake_vmax,
            'Ea': ea_values,
            'Vmax0': vmax_values,
            'EnzAttrib': enz_attrib_df,

            # Microbes
            'Microbes_pp': microbes_pp,
            'fb': fb_array,
            'Microbes': microbes_df,
            'Bac_density': bac_density,
            'Fun_density': fun_density,

            # Microbial traits and genes (all expanded)
            'MinRatios': min_ratios,
            'UptakeGenes': uptake_genes_expanded,
            'OsmoGenes': osmolyte_genes_expanded,
            'EnzGenes': enzyme_genes_expanded,
            'UptakeGenes_trait': uptake_genes_trait_expanded,
            'OsmoProdConsti_trait': osmo_prod_consti_trait_expanded,
            'OsmoProdInduci_trait': osmo_prod_induci_trait_expanded,
            'EnzProdConsti_trait': enz_prod_consti_trait_expanded,
            'EnzProdInduci_trait': enz_prod_induci_trait_expanded,
            'UptakeGenesCost': uptake_genes_cost,  # Already expanded above
            'OsmoProdConsti': osmo_prod_consti,    # Already expanded above
            'OsmoProdInduci': osmo_prod_induci,    # Already expanded above
            'EnzProdConstit': enz_prod_constit,    # Already expanded above
            'EnzProdInduce': enz_prod_induce,      # Already expanded above
            'TaxDroughtTol': tax_drought_tol,      # Already expanded above

            # Mortality
            'basal_death_prob': basal_death_prob,
            'death_rate': death_rate,

            # Metabolism
            'AE_ref': ae_ref,
            'AE_temp': ae_temp,
            'Uptake_Maint_cost': params['Uptake_Maint_cost'],

            # Thresholds
            'C_min': params['C_min'],
            'N_min': params['N_min'],
            'P_min': params['P_min'],
            'max_size_b': params['max_size_b'],
            'max_size_f': params['max_size_f'],

            # Environmental
            'wp_fc': params['wp_fc'],
            'wp_th': params['wp_th'],
            'alpha': params['alpha'],

            # Climate
            'Temp': temp_series,
            'Psi': psi_series
        }

        return data_init

    @staticmethod
    def create_grid(
        end_time: int = 365,
        gridsize: int = 100,
        n_taxa: int = 5,
        n_substrates: int = 3,
        substrate_c: float = 500.0,
        substrate_n: float = 20.0,
        substrate_p: float = 5.0,
        microbial_c: float = 50.0,
        temperature: float = 15.0,
        moisture_psi: float = -0.5,
        fb: float = 0.1,
        x: int = 10,
        y: int = 10,
        parameters_dict: Dict = None,
        use_file_based_init: bool = True,
        **kwargs
    ) -> Grid:
        """
        Create a Grid instance with complete parameters.

        This is the main entry point for creating a DEMENTpy Grid
        with minimal external file dependencies.

        Args:
            end_time: Number of days to simulate
            gridsize: Total number of grid cells (default: 100 = 10x10)
            n_taxa: Number of microbial taxa
            n_substrates: Number of substrate types
            substrate_c: Initial substrate C (g/m²)
            substrate_n: Initial substrate N (g/m²)
            substrate_p: Initial substrate P (g/m²)
            microbial_c: Initial microbial biomass C (g/m²)
            temperature: Initial temperature (°C)
            moisture_psi: Water potential (MPa)
            fb: Fraction of fungal taxa (0-1)
            x: X dimension of grid (default: 10)
            y: Y dimension of grid (default: 10)
            parameters_dict: Optional dictionary to override default parameters
            use_file_based_init: Use DEMENTpy's file-based initialization (recommended)
            **kwargs: Additional parameters passed to runtime configuration

        Returns:
            Initialized Grid instance ready for simulation
        """
        # Ensure gridsize matches x*y
        if gridsize != x * y:
            gridsize = x * y

        if use_file_based_init:
            # Hybrid approach: Write temporary files and use DEMENTpy's proven initialization
            return DEMENTLibrary._create_grid_from_temp_files(
                end_time=end_time,
                gridsize=gridsize,
                x=x,
                y=y,
                n_taxa=n_taxa,
                n_substrates=n_substrates,
                substrate_c=substrate_c,
                substrate_n=substrate_n,
                substrate_p=substrate_p,
                microbial_c=microbial_c,
                temperature=temperature,
                moisture_psi=moisture_psi,
                fb=fb,
                parameters_dict=parameters_dict
            )
        else:
            # Direct initialization (experimental - may have index issues)
            # n_substrates_total includes DeadMic and DeadEnz (2 special substrates)
            n_substrates_total = n_substrates + 2
            runtime = DEMENTLibrary.create_default_runtime(
                end_time=end_time,
                gridsize=gridsize,
                n_taxa=n_taxa,
                n_substrates=n_substrates_total,
                n_enzymes=n_substrates,
                x=x,
                y=y,
                **kwargs
            )

            data_init = DEMENTLibrary.create_minimal_initialization(
                gridsize=gridsize,
                n_taxa=n_taxa,
                n_substrates=n_substrates,
                substrate_c=substrate_c,
                substrate_n=substrate_n,
                substrate_p=substrate_p,
                microbial_c=microbial_c,
                temperature=temperature,
                moisture_psi=moisture_psi,
                fb=fb,
                parameters_dict=parameters_dict
            )

            grid = Grid(runtime, data_init)
            return grid

    @staticmethod
    def _create_grid_from_temp_files(
        end_time: int,
        gridsize: int,
        x: int,
        y: int,
        n_taxa: int,
        n_substrates: int,
        substrate_c: float,
        substrate_n: float,
        substrate_p: float,
        microbial_c: float,
        temperature: float,
        moisture_psi: float,
        fb: float,
        parameters_dict: Dict
    ) -> Grid:
        """
        Create Grid using temporary parameter files.

        This hybrid approach uses our parameter mapping but leverages
        DEMENTpy's proven file-based initialization.
        """
        # Create temporary directory
        temp_dir = tempfile.mkdtemp(prefix='dementpy_')

        try:
            # Write parameter files
            DEMENTLibrary._write_parameter_files(
                temp_dir,
                gridsize=gridsize,
                n_taxa=n_taxa,
                n_substrates=n_substrates,
                substrate_c=substrate_c,
                substrate_n=substrate_n,
                substrate_p=substrate_p,
                microbial_c=microbial_c,
                temperature=temperature,
                moisture_psi=moisture_psi,
                parameters_dict=parameters_dict
            )

            # Create runtime DataFrame
            # n_substrates_total includes DeadMic and DeadEnz (2 special substrates)
            n_substrates_total = n_substrates + 2
            runtime = DEMENTLibrary.create_default_runtime(
                end_time=end_time,
                gridsize=gridsize,
                n_taxa=n_taxa,
                n_substrates=n_substrates_total,  # Total includes special substrates
                n_enzymes=n_substrates,  # Only organic substrates need enzymes
                x=x,
                y=y
            )

            # Use DEMENTpy's initialize_data function
            data_init = initialize_data(runtime, temp_dir)

            # Create Grid
            grid = Grid(runtime, data_init)

            return grid

        finally:
            # Clean up temporary files
            shutil.rmtree(temp_dir, ignore_errors=True)

    @staticmethod
    def _write_parameter_files(
        temp_dir: str,
        gridsize: int,
        n_taxa: int,
        n_substrates: int,
        substrate_c: float,
        substrate_n: float,
        substrate_p: float,
        microbial_c: float,
        temperature: float,
        moisture_psi: float,
        parameters_dict: Dict
    ):
        """
        Write minimal parameter files needed by DEMENTpy's initialize_data.
        """
        # Get default parameters
        if parameters_dict is None:
            parameters_dict = {}

        params = {
            'max_size_b': 2.0,
            'max_size_f': 50.0,
            'Cfrac_b': 0.825,
            'Nfrac_b': 0.16,
            'Pfrac_b': 0.015,
            'Cfrac_f': 0.9,
            'Nfrac_f': 0.09,
            'Pfrac_f': 0.01,
            'Crange': 0.09,
            'Nrange': 0.04,
            'Prange': 0.005,
            'C_min': 0.086,
            'N_min': 0.012,
            'P_min': 0.002,
            'Uptake_C_cost_min': 0.01,
            'Uptake_C_cost_max': 0.1,
            'Uptake_Maint_cost': 0.01,
            'Enz_per_taxon_min': 0,
            'Enz_per_taxon_max': 40,
            'Enz_Prod_min': 0.00001,
            'Enz_Prod_max': 0.0001,
            'Constit_Prod_min': 0.00001,
            'Constit_Prod_max': 0.0001,
            'Osmo_per_taxon_min': 1,
            'Osmo_per_taxon_max': 1,
            'Osmo_Consti_Prod_min': 0.0000001,
            'Osmo_Consti_Prod_max': 0.000001,
            'Osmo_Induci_Prod_min': 0.01,
            'Osmo_Induci_Prod_max': 0.1,
            'CUE_ref': 0.5,
            'CUE_temp': -0.005,
            'death_rate_bac': 0.001,
            'death_rate_fun': 0.001,
            'beta_bac': 10,
            'beta_fun': 10,
            'wp_fc': -2.0,
            'wp_th': -6.0,
            'alpha': 0.01,
            'Enz_C_cost': 1,
            'Enz_N_cost': 0.3,
            'Enz_P_cost': 0,
            'Enz_Maint_cost': 5,
            'Vmax0_min': 5.0,
            'Vmax0_max': 50.0,
            'Uptake_Vmax0_min': 1.0,
            'Uptake_Vmax0_max': 10.0,
            'Uptake_Ea_min': 35.0,
            'Uptake_Ea_max': 35.0,
            'Km_min': 0.01,
            'Uptake_Km_min': 0.001,
            'Vmax_Km': 1,
            'Vmax_Km_int': 0,
            'Uptake_Vmax_Km': 0.2,
            'Uptake_Vmax_Km_int': 0,
            'Km_error': 0,
            'Specif_factor': 1
        }
        params.update(parameters_dict)

        # 1. Write parameters.csv
        params_df = pd.DataFrame(list(params.items()))
        params_df.to_csv(os.path.join(temp_dir, 'parameters.csv'), header=False, index=False)

        # 2. Write initial_substrates.csv
        # IMPORTANT: Must include DeadMic and DeadEnz for recycling dead biomass
        special_substrates = ['DeadMic', 'DeadEnz']
        organic_substrate_names = [f'Sub{i+1}' for i in range(n_substrates)]
        all_substrate_names = special_substrates + organic_substrate_names

        c_per_sub = substrate_c / n_substrates if n_substrates > 0 else 0
        n_per_sub = substrate_n / n_substrates if n_substrates > 0 else 0
        p_per_sub = substrate_p / n_substrates if n_substrates > 0 else 0

        # DeadMic and DeadEnz start at 0; organic substrates get the initial carbon
        substrate_c_values = [0, 0] + [c_per_sub] * n_substrates
        substrate_n_values = [0, 0] + [n_per_sub] * n_substrates
        substrate_p_values = [0, 0] + [p_per_sub] * n_substrates

        substrates_df = pd.DataFrame({
            '': all_substrate_names,
            'C': substrate_c_values,
            'N': substrate_n_values,
            'P': substrate_p_values
        })
        substrates_df.to_csv(os.path.join(temp_dir, 'initial_substrates.csv'), index=False)

        # 3. Write sub_mon_inputs.csv (substrate and monomer inputs)
        # Format: first column is substrate names, then 'Sub' and 'Mon' columns
        # Only organic substrates get external inputs (not DeadMic/DeadEnz)
        inputs_df = pd.DataFrame({
            '': organic_substrate_names,
            'Sub': [0.1] * n_substrates,  # Substrate input rates
            'Mon': [0.01] * n_substrates  # Monomer input rates
        })
        inputs_df.to_csv(os.path.join(temp_dir, 'sub_mon_inputs.csv'), index=False)

        # 4. Write enzyme_ea.csv (activation energies)
        # Must include all substrates to match substrate count
        # DeadMic and DeadEnz don't require enzymatic degradation but need entries
        ea_values_min = [0.0, 0.0] + [37.0] * n_substrates  # DeadMic/DeadEnz get 0
        ea_values_max = [0.0, 0.0] + [37.0] * n_substrates
        ea_df = pd.DataFrame({
            '': all_substrate_names,
            'Ea_min': ea_values_min,
            'Ea_max': ea_values_max
        })
        ea_df.to_csv(os.path.join(temp_dir, 'enzyme_ea.csv'), index=False)

        # 5. Write climate.csv (365 days of temperature and moisture)
        climate_df = pd.DataFrame({
            'Day': range(1, 366),
            'Temp': [temperature] * 365,
            'Psi': [moisture_psi] * 365
        })
        climate_df.to_csv(os.path.join(temp_dir, 'climate.csv'), index=False)

    @staticmethod
    def run_timestep(grid: Grid, day: int) -> Dict[str, float]:
        """
        Run one daily timestep of DEMENTpy.

        Args:
            grid: Grid instance
            day: Day number (0-364)

        Returns:
            Dictionary with key outputs from this timestep
        """
        # Run all grid methods for one day
        grid.degradation(day)
        grid.uptake(day)
        grid.metabolism(day)
        grid.mortality(day)
        grid.reproduction(day)

        # Extract outputs
        outputs = DEMENTLibrary.extract_outputs(grid)

        return outputs

    @staticmethod
    def extract_outputs(grid: Grid) -> Dict[str, float]:
        """
        Extract coupling-relevant outputs from Grid.

        Args:
            grid: Grid instance

        Returns:
            Dictionary with:
                - available_N: Available N for plant uptake (g/m²)
                - respiration: CO2 respiration (g C/m²)
                - microbial_biomass_c: Total microbial C (g/m²)
                - microbial_biomass_n: Total microbial N (g/m²)
                - substrate_c: Remaining substrate C (g/m²)
                - cue_system: Microbial carbon use efficiency
        """
        # Aggregate across grid cells
        total_substrate_c = grid.Substrates.values.sum()

        # Get NH4 (ammonium) which is plant-available mineral nitrogen
        # Monomers DataFrame has rows: NH4, PO4, DeadMic, DeadEnz, Mon3, etc.
        # and columns: C, N, P
        is_NH4 = grid.Monomers.index == 'NH4'
        total_monomer_n = grid.Monomers.loc[is_NH4, 'N'].sum() if is_NH4.any() else 0.0

        total_microbe_c = grid.Microbes['C'].sum()
        total_microbe_n = grid.Microbes['N'].sum()

        # Respiration (if available from last metabolism call)
        respiration = float(grid.Respiration) if not np.isnan(grid.Respiration) else 0.0

        # CUE (if available)
        cue = float(grid.CUE_system) if not np.isnan(grid.CUE_system) else 0.5

        # Available N = NH4 ammonium (plant-available mineral N)
        available_n = total_monomer_n

        outputs = {
            'available_N': available_n,
            'respiration': respiration,
            'microbial_biomass_c': total_microbe_c,
            'microbial_biomass_n': total_microbe_n,
            'substrate_c': total_substrate_c,
            'cue_system': cue
        }

        return outputs

    @staticmethod
    def get_state(grid: Grid) -> Dict[str, Any]:
        """
        Extract complete state from Grid for persistence.

        Args:
            grid: Grid instance

        Returns:
            State dictionary that can be saved and restored
        """
        state = {
            'Substrates': grid.Substrates.copy(),
            'Enzymes': grid.Enzymes.copy(),
            'Monomers': grid.Monomers.copy(),
            'Microbes': grid.Microbes.copy(),
            'Respiration': grid.Respiration,
            'CUE_system': grid.CUE_system,
            # Add other state variables as needed
        }

        return state

    @staticmethod
    def set_state(grid: Grid, state: Dict[str, Any]):
        """
        Restore Grid state from saved state dictionary.

        Args:
            grid: Grid instance to update
            state: State dictionary from get_state()
        """
        grid.Substrates = state['Substrates'].copy()
        grid.Enzymes = state['Enzymes'].copy()
        grid.Monomers = state['Monomers'].copy()
        grid.Microbes = state['Microbes'].copy()
        grid.Respiration = state['Respiration']
        grid.CUE_system = state['CUE_system']


# Convenience function for quick testing
def test_library_interface():
    """Test the library interface."""
    print("Testing DEMENTpy Library Interface\n")

    print("1. Creating Grid...")
    grid = DEMENTLibrary.create_grid(
        end_time=10,
        gridsize=5,
        n_taxa=3,
        substrate_c=100.0,
        microbial_c=10.0
    )
    print(f"   ✓ Grid created: {grid.gridsize}x{grid.gridsize}, {grid.n_taxa} taxa")

    print("\n2. Running 10 timesteps...")
    for day in range(10):
        outputs = DEMENTLibrary.run_timestep(grid, day)
        if day == 0 or day == 9:
            print(f"   Day {day}: Resp={outputs['respiration']:.4f}, "
                  f"AvailN={outputs['available_N']:.4f}, "
                  f"MicC={outputs['microbial_biomass_c']:.2f}")

    print("\n3. Testing state save/restore...")
    state = DEMENTLibrary.get_state(grid)
    print(f"   ✓ State saved ({len(state)} components)")

    # Create new grid and restore state
    grid2 = DEMENTLibrary.create_grid(end_time=10, gridsize=5, n_taxa=3)
    DEMENTLibrary.set_state(grid2, state)
    print(f"   ✓ State restored to new grid")

    print("\n✅ Library interface test complete!")


if __name__ == "__main__":
    test_library_interface()
