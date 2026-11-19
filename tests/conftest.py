import pytest
import numpy as np

@pytest.fixture(autouse=True)
def reseed_numpy():
    np.random.seed(123)
