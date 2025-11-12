"""Tests for the initialization of the data."""

import pandas as pd
import pytest

from initialization import initialize_data


def test_initialize_data():
    """Test that initialize_data works without an error."""
    input_dir = 'grassland'
    runtime = pd.read_csv(input_dir+'/runtime.txt', header=None, index_col=0, sep='\t')
    data = initialize_data(runtime, input_dir)
    assert isinstance(data, dict)
