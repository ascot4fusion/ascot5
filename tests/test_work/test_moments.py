import numpy as np

from a5py.data.hist import Dist
from a5py.work.moments import density

def test_density():
    dist = Dist(
        np.full((3,4,5), 1),
        {
            "r": np.arange(4),
            "z": np.arange(5),
            "ekin": np.arange(6),
        }
    )
    density(dist, np.full((3,4), 1.))