import pytest

import numpy as np

from a5py.composite.run import Run

from a5py.data import AscotData
from a5py.data.marker.fieldline import FieldlineMarker

def test_run():
    a5 = AscotData()
    mrk = a5.create_fieldlinemarker(r=np.full(3, 6.2), z=np.full(3, 0.0))
    run = Run()
    run._setup(mrk._cdata)
