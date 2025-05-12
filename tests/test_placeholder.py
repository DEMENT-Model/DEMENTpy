"""Tests to show how to use pytest."""

import numpy as np
import pytest


def test_basic():
    """Illustrates how to write a simple test."""
    assert 1 + 1 == 2


def test_raises_exception():
    """Illustrates how to write a test that checks for an exception."""
    with pytest.raises(ZeroDivisionError):
        1 / 0


@pytest.mark.parametrize("x,y,r", [(3, 4, 5), (-5, -12, 13)])
def test_multiple_parameters(x, y, r):
    """Illustrates how to write a test to run with multiple input arguments."""
    assert np.sqrt(x**2 + y**2) == r
