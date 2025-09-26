import pytest
import numpy as np

from a5py.data.marker.state import MarkerState, Structure


def test_init():
    state = MarkerState()
    with pytest.raises(TypeError):
        state.getval("r")

    n = 10
    mrk = (n * Structure)()
    print(getattr(mrk[2], "r"))
    MarkerState.from_params(mrk)
    assert isinstance(state, MarkerState)