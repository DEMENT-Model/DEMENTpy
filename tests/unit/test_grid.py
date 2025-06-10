"""Tests for the Grid class."""

import pandas as pd
import pytest

from initialization import initialize_data
from grid import Grid


@pytest.fixture
def grid():
    """Initialize a Grid object using the initialize_data function."""
    input_dir = 'grassland'
    runtime = pd.read_csv(input_dir+'/runtime.txt', header=None, index_col=0, sep='\t')
    data = initialize_data(runtime, input_dir)
    return Grid(runtime, data)


def test_initialization_runs(grid):
    """Test that the initialization of Grid works without an error."""
    pass


def test_degradation_runs(grid):
    """Test that Grid.degradation works without an error."""
    grid.degradation(0)


def test_uptake_runs(grid):
    """Test that Grid.uptake works without an error."""
    grid.degradation(0) # degradation needs to run to initialize some DataFrames
    grid.uptake(0)


def test_metabolism_runs(grid):
    """Test that Grid.metabolism works without an error."""
    grid.metabolism(0)


def test_mortality_runs(grid):
    """Test that Grid.mortality works without an error."""
    grid.mortality(0)


def test_reproduction_runs(grid):
    """Test that Grid.reproduction works without an error."""
    grid.reproduction(0)
