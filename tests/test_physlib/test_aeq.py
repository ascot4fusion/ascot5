import pytest
import numpy as np

import a5py.physlib as physlib

@pytest.mark.parametrize(
        ("A", "epsilon", "kappa", "delta", "xpoint", "symmetric", "axisr"),
        [
            (1.0, 0.32, 1.0, 0.0, None, True, 1.01),
            (-0.155, 0.32, 1.7, 0.33, (0.88, -0.6), False, 1.05),
            # The X-point is from the paper: (1-1.1*delta*eps, 1.1*kappa*eps)
            (0.0, 0.78, 2.0, 0.35, (1.0 - 1.1*0.35*0.78, 1.1*2.0*0.78), True,
             1.27),
        ],
        ids=("iter-circular", "iter-baseline", "nstx-doublenull"),
    )
def test_analyticalequilibrium(
    A, epsilon, kappa, delta, xpoint, symmetric, axisr,
    ):
    """Test the analytical equilibrium model.

    Here we make use of the fact that psi = 0 at the LCFS.
    """
    c = physlib.aeq.parameters2coefficients(
        A, epsilon, kappa, delta, xpoint, symmetric,
        )
    rmin, rmax, zmin, zmax = (
        np.amax((1.0 - epsilon - 0.1, 1e-8)),
        1.0 + epsilon + 0.1,
        -1.2*kappa*epsilon - 0.1,
        1.2*kappa*epsilon + 0.1,
    )
    def psi(x, y):
        return physlib.aeq.psi(x, y, c, A)

    if symmetric:
        X1, Y1 = np.meshgrid(np.linspace(rmin, rmax, 200),
                           np.linspace(0.0, zmax, 200))
        X2, Y2 = np.meshgrid(np.linspace(rmin, rmax, 200),
                           np.linspace(0.0, zmin, 200))
        assert np.isclose(psi(X1, Y1), psi(X2, Y2)).all()

    high_point_inside = (1.0 - epsilon*delta, kappa*epsilon - 1e-3)
    high_point_outside = (1.0 - epsilon*delta, kappa*epsilon + 1e-3)
    assert psi(*high_point_inside) < 0.0
    if symmetric and xpoint is not None:
        # Double-null
        assert psi(*high_point_outside) < 0.0
    else:
        assert psi(*high_point_outside) > 0.0

    assert np.isclose(
        physlib.aeq.find_axis(np.append(c, A))[0],
        axisr,
        atol=1e-2
        )
