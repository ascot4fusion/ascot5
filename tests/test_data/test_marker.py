import unyt
import numpy as np
import pytest

from a5py import Ascot
from a5py.data.marker.state import MarkerState

@pytest.mark.parametrize(
        "marker",
        ["FieldLine", "Guidingcenter", "Particle"],
        )
def test_solve_setup_markers(marker):
    ascot = Ascot()
    arr = np.ones(100,)
    if marker == "FieldLine":
        mrk = ascot.data.create_fieldlinemarker(r=arr*6.3*unyt.m, z=0.0*unyt.m)
    elif marker == "Particle":
        mrk = ascot.data.create_particlemarker(
            species="proton", charge=1*unyt.e, r=arr*6.3*unyt.m, z=0.0*unyt.m,
            vr=0.0*unyt.m/unyt.s, vphi=1.0e6*unyt.m/unyt.s, vz=0.0*unyt.m/unyt.s,
            )
    elif marker == "Guidingcenter":
        mrk = ascot.data.create_guidingcentermarker(
            species="proton", charge=1*unyt.e, r=arr*6.3*unyt.m, z=0.0*unyt.m,
            ekin=1.0*unyt.eV, pitch=0.0,
            )
    bjacb = np.full((mrk.n, 12), 1.0)
    obj = MarkerState.from_markerinput(mrk, bjacb, idx=np.array([0]))
    assert obj.ids == 1
