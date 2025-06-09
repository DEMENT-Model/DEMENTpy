"""Test(s) to use with pytest-profiling to identify bottlenecks in the code."""
import pytest

@pytest.mark.profiling
def test_profiling(monkeypatch):

    # Define the command line arguments
    import sys
    argv = ['dementpy.py', 'grassland', 'output', '20250402', 'scrubland']
    monkeypatch.setattr(sys, 'argv', argv)

    # Move to subfolder so input and output folders will be correct
    import os
    os.chdir('src')

    # Run dementpy
    import dementpy
    dementpy.main()
