"""Test magnetic field module."""
import time

import unyt
import pytest
import numpy as np
from scipy.interpolate import make_interp_spline

from a5py import Ascot
from a5py.physlib.formulas import cart2pol, cart2pol_vec
from a5py.engine.interpolate import evaluate, qnt_mapping
from a5py.libascot import LIBASCOT
from a5py.templates import (
    PremadeMagneticField, TransmuteAnalyticalBfieldToSplines,
)
from a5py.physlib import aeq

from .conf import test_data, create_input

@pytest.fixture
def ascot():
    return Ascot()


def test_bfieldcartesian_init(ascot):
    obj = ascot.data.create_bfieldcartesian(
        bxyz=test_data["BfieldCartesian"]["bxyz"],
        jacobian=test_data["BfieldCartesian"]["jacobian"],
        axisrz=test_data["BfieldCartesian"]["axisrz"],
        rhoval=test_data["BfieldCartesian"]["rhoval"],
        )
    assert obj.psival.v == obj.rhoval.v


def test_interpolate_cartesian():
    jac = np.array([
        [1., 2., 3.],
        [4., 5., 6.],
        [7., 8., 9.],
        ]) * unyt.T/unyt.m
    obj = create_input("BfieldCartesian", jacobian=jac)
    x, y, z = 2.*unyt.m, 0.*unyt.m, 0.*unyt.m
    r, phi, z = cart2pol(x, y, z)
    br0, bphi0, bz0, rho, psi = evaluate(
        r, phi, z, 0.*unyt.s, "br", "bphi", "bz", "rho", "psi", bfield=obj
        )

    jac = obj.jacobian
    bx = obj.bxyz[0] + jac[0,0]*x + jac[0,1]*y + jac[0,2]*z
    by = obj.bxyz[1] + jac[1,0]*x + jac[1,1]*y + jac[1,2]*z
    bz = obj.bxyz[2] + jac[2,0]*x + jac[2,1]*y + jac[2,2]*z

    br, bphi, bz = cart2pol_vec(x, bx, y, by, z, bz)
    assert np.isclose(br, br0)
    assert np.isclose(bz, bz0)
    assert np.isclose(bphi, bphi0)
    assert np.isclose(rho, obj.rhoval)
    assert np.isclose(psi, obj.psival)

    jac0 = evaluate(
        r, phi, z, 0.*unyt.s, "brdr", "brdphi", "brdz", "bphidr", "bphidphi",
        "bphidz", "bzdr", "bzdphi", "bzdz", bfield=obj
        )
    dr = 1.e-1 * unyt.m
    br, bphi, bz = evaluate(
        r+dr, phi, z, 0.*unyt.s, "br", "bphi", "bz", bfield=obj
        )
    assert np.isclose(br - br0, jac0[0] * dr)
    assert np.isclose(bphi - bphi0, jac0[3] * dr)
    assert np.isclose(bz - bz0, jac0[6] * dr)

    dz = 1.e-1 * unyt.m
    br, bphi, bz = evaluate(
        r, phi, z+dz, 0.*unyt.s, "br", "bphi", "bz", bfield=obj
        )
    assert np.isclose(br - br0, jac0[2] * dz)
    assert np.isclose(bphi - bphi0, jac0[5] * dz)
    assert np.isclose(bz - bz0, jac0[8] * dz)

    dphi = 1.e-3 * unyt.deg
    br, bphi, bz = evaluate(
        r, phi+dphi, z, 0.*unyt.s, "br", "bphi", "bz", bfield=obj
        )
    dphi = dphi.to("rad").v
    assert np.isclose(br - br0, jac0[1] * dphi)
    assert np.isclose(bphi - bphi0, jac0[4] * dphi)
    assert np.isclose(bz - bz0, jac0[7] * dphi)


def test_interpolate_analytical(ascot):
    template = PremadeMagneticField(ascot, field="iter-baseline")
    obj = template.create_input()

    n = 100
    rng = np.random.default_rng(123)
    r = (4 + (8.2 - 4) * rng.random(n))*unyt.m
    z = (-2 + (2 - -2) * rng.random(n))*unyt.m

    rho0, psi0 = evaluate(
        r, 0*unyt.deg, z, 0.*unyt.s, "rho", "psi", bfield=obj
        )
    r0 = obj.rmajor
    psi = aeq.psi(r/r0, z/r0, obj.coefficients[:-1], obj.coefficients[-1]) * obj.psiscaling
    assert np.isclose(psi, psi0).all()
    rho = np.sqrt( (psi - obj.psilimits[0]) / (obj.psilimits[1] - obj.psilimits[0]) )
    assert np.isclose(rho, rho0).all()

    r0, z0, psi0 = r[0], z[0], psi[0]
    dr, dz = 1.e-4 * unyt.m, 1.e-4 * unyt.m
    psi = evaluate(
        unyt.unyt_array([r0+dr, r0]), 0*unyt.deg, unyt.unyt_array([z0, z0+dz]),
        0.*unyt.s, "psi", bfield=obj
        )

    bz = (unyt.rad*(psi[0] - psi0) / (r0*dr)).to("T")
    br = -(unyt.rad*(psi[1] - psi0) / (r0*dz)).to("T")
    br0, bz0 = evaluate(r0, 0*unyt.deg, z0, 0.*unyt.s, "br", "bz", bfield=obj)
    assert np.isclose(br, br0, atol=1e-3)
    assert np.isclose(bz, bz0, atol=1e-3)

    r, z, phi = r0, z0, 0.0*unyt.deg
    br0, bphi0, bz0 = evaluate(
        r, phi, z, 0.*unyt.s, "br", "bphi", "bz", bfield=obj
        )

    jac0 = evaluate(
        r, phi, z, 0.*unyt.s, "brdr", "brdphi", "brdz", "bphidr", "bphidphi",
        "bphidz", "bzdr", "bzdphi", "bzdz", bfield=obj
        )
    dr = 1.e-4 * unyt.m
    br, bphi, bz = evaluate(
        r+dr, phi, z, 0.*unyt.s, "br", "bphi", "bz", bfield=obj
        )
    assert np.isclose(br - br0, jac0[0] * dr)
    assert np.isclose(bphi - bphi0, jac0[3] * dr)
    assert np.isclose(bz - bz0, jac0[6] * dr)

    dz = 1.e-4 * unyt.m
    br, bphi, bz = evaluate(
        r, phi, z+dz, 0.*unyt.s, "br", "bphi", "bz", bfield=obj
        )
    assert np.isclose(br - br0, jac0[2] * dz)
    assert np.isclose(bphi - bphi0, jac0[5] * dz)
    assert np.isclose(bz - bz0, jac0[8] * dz)

    dphi = 1.e-4 * unyt.deg
    br, bphi, bz = evaluate(
        r, phi+dphi, z, 0.*unyt.s, "br", "bphi", "bz", bfield=obj
        )
    dphi = dphi.to("rad").v
    assert np.isclose(br - br0, jac0[1] * dphi)
    assert np.isclose(bphi - bphi0, jac0[4] * dphi)
    assert np.isclose(bz - bz0, jac0[7] * dphi)


def _phi_interpolant(variant, nodes, values, period):
    """Return the scheme's periodic 1D interpolant of values tabulated at nodes
    along phi.

    The linear scheme interpolates piecewise linearly and the spline one with
    the (unique) periodic cubic spline. The returned function takes phi in
    degrees and returns the nu:th derivative per radian.
    """
    k = 1 if variant == "BfieldLinear3D" else 3
    x = np.deg2rad(np.append(nodes, nodes[0] + period))
    spline = make_interp_spline(
        x, np.append(values, values[0]), k=k, bc_type="periodic",
        )
    # bc_type doesn't make the k=1 spline wrap around, so set it explicitly
    spline.extrapolate = "periodic"
    return lambda phi, nu=0: spline(np.deg2rad(phi), nu)


@pytest.mark.parametrize("variant", ["BfieldLinear3D", "BfieldSpline3D"])
@pytest.mark.parametrize("period", [360, 20], ids=["full", "partial"])
def test_interpolate_3d(variant, period):
    """Test 3D fields against data that their interpolation reproduces exactly.

    The poloidal flux and the tabulated field components are bilinear in
    (R, z), which both schemes reproduce exactly, and the field components are
    also scaled by arbitrary values tabulated along phi, so the result is known
    exactly in terms of the 1D interpolant along phi. This tests the periodic
    wrap-around with the period being the full torus or a fraction of it (e.g.
    18 TF coils). The R and z grids of psi differ from those of B to verify
    that the two are handled separately.
    """
    rgrid = np.linspace(4.0, 8.0, 21) * unyt.m
    zgrid = np.linspace(-2.0, 2.0, 41) * unyt.m
    phigrid = np.linspace(0, period, 7)[:-1] * unyt.deg
    rgridpsi = np.linspace(3.5, 8.5, 31) * unyt.m
    zgridpsi = np.linspace(-2.5, 2.5, 61) * unyt.m

    p = np.array([1.0, 0.3, -0.2, 0.05]) # psi = p0 + p1 R + p2 z + p3 R z
    a = np.array([0.1, 0.02, -0.03, 0.01]) # br perturbation
    b = np.array([5.0, -0.4, 0.1, 0.02]) # bphi
    c = np.array([-0.2, 0.01, 0.04, -0.03]) # bz perturbation
    bilinear = lambda k, r, z: k[0] + k[1]*r + k[2]*z + k[3]*r*z
    ddr = lambda k, r, z: k[1] + k[3]*z
    ddz = lambda k, r, z: k[2] + k[3]*r
    g = np.random.default_rng(3).uniform(0.5, 1.5, phigrid.size)

    rp, zp = np.meshgrid(rgridpsi.v, zgridpsi.v, indexing="ij")
    rb, gb, zb = np.meshgrid(rgrid.v, g, zgrid.v, indexing="ij")
    psilimits = unyt.unyt_array([-1.0, 4.0], "Wb/rad")
    obj = create_input(
        variant, rgrid=rgrid, phigrid=phigrid, zgrid=zgrid,
        rgridpsi=rgridpsi, zgridpsi=zgridpsi,
        psi=bilinear(p, rp, zp) * unyt.Wb/unyt.rad,
        br=bilinear(a, rb, zb) * gb * unyt.T,
        bphi=bilinear(b, rb, zb) * gb * unyt.T,
        bz=bilinear(c, rb, zb) * gb * unyt.T, psilimits=psilimits,
        )
    assert np.allclose(obj.phigrid.to("deg").v, phigrid.v)

    # Points inside the (R, z) grid and up to two periods outside the
    # tabulated one along phi
    rng = np.random.default_rng(1)
    n = 200
    r = 4.1 + 3.8 * rng.random(n)
    z = -1.9 + 3.8 * rng.random(n)
    phi = period * rng.uniform(-2, 3, n)
    br, bphi, bz, psi, rho, *jac = evaluate(
        r*unyt.m, phi*unyt.deg, z*unyt.m, 0.*unyt.s, "br", "bphi", "bz", "psi",
        "rho", *qnt_mapping["bjac"], bfield=obj,
        )
    interpolant = _phi_interpolant(variant, phigrid.v, g, period)
    G, dG = interpolant(phi), interpolant(phi, nu=1)

    # B = B_tabulated + grad(psi) x e_phi / R
    psi0 = bilinear(p, r, z)
    psidr, psidz, psidrdz = ddr(p, r, z), ddz(p, r, z), p[3]
    assert np.allclose(psi.v, psi0)
    assert np.allclose(
        rho.v, np.sqrt((psi0 - psilimits[0].v) / np.diff(psilimits.v)),
        )
    assert np.allclose(br.v, bilinear(a, r, z) * G - psidz / r)
    assert np.allclose(bphi.v, bilinear(b, r, z) * G)
    assert np.allclose(bz.v, bilinear(c, r, z) * G + psidr / r)

    # Second derivatives of psi other than the mixed one vanish
    expected = {
        "brdr": ddr(a, r, z) * G + psidz / r**2 - psidrdz / r,
        "brdphi": bilinear(a, r, z) * dG,
        "brdz": ddz(a, r, z) * G,
        "bphidr": ddr(b, r, z) * G,
        "bphidphi": bilinear(b, r, z) * dG,
        "bphidz": ddz(b, r, z) * G,
        "bzdr": ddr(c, r, z) * G - psidr / r**2,
        "bzdphi": bilinear(c, r, z) * dG,
        "bzdz": ddz(c, r, z) * G + psidrdz / r,
        }
    for name, value in zip(qnt_mapping["bjac"], jac):
        assert np.allclose(value.v, expected[name], atol=1e-12), name


@pytest.mark.parametrize("variant", ["BfieldLinear3D", "BfieldSpline3D"])
def test_interpolate_3d_outside_grid(variant):
    """Test that the evaluation outside the grid fails and NaN is returned."""
    obj = create_input(variant)
    br = evaluate(
        9.0*unyt.m, 0*unyt.deg, 0.0*unyt.m, 0.*unyt.s, "br", bfield=obj,
        )
    assert np.isnan(br)


@pytest.fixture(scope="module")
def tabulated_fields():
    """ITER-like circular field tabulated with each scheme.

    Returns the analytical field and {name: field tabulated from it}. The 3D
    fields are axisymmetric copies of the 2D one.
    """
    ascot = Ascot()
    reference = PremadeMagneticField(
        ascot, field="iter-circular", axisymmetric=True,
        ).create_input()
    spline2d = TransmuteAnalyticalBfieldToSplines(
        ascot, field=reference, nr=60, rlim=[4.0, 8.5] * unyt.m,
        nz=100, zlim=[-4.0, 4.0] * unyt.m,
        ).create_input()
    phigrid = np.linspace(0, 360, 5)[:-1] * unyt.deg
    bphi = np.repeat(spline2d.bphi[:, None, :], phigrid.size, axis=1)
    data3d = {
        "rgrid": spline2d.rgrid, "phigrid": phigrid, "zgrid": spline2d.zgrid,
        "psi": spline2d.psi, "bphi": bphi, "br": 0*bphi, "bz": 0*bphi,
        "axisrz": reference.axisrz, "psilimits": reference.psilimits,
        "preview": True,
        }
    return reference, {
        "Spline2D": spline2d,
        "Spline3D": ascot.data.create_bfieldspline3d(**data3d),
        "Linear3D": ascot.data.create_bfieldlinear3d(**data3d),
        }


def _plasma_points(n, seed):
    """Random points uniformly distributed in the poloidal plane inside
    r/a = 0.9 of the ITER-like circular plasma.

    The plasma boundary is left out, since on the inboard side it is only
    about 3 cells from the R grid edge, where the natural boundary condition of
    the splines degrades their accuracy.
    """
    # Major radius and inverse aspect ratio of "iter-circular" in
    # PremadeMagneticField
    rmajor, epsilon = 6.2, 0.32
    rng = np.random.default_rng(seed)
    rminor = rmajor * epsilon * 0.9 * np.sqrt(rng.random(n))
    angle = rng.uniform(-np.pi, np.pi, n)
    r = (rmajor + rminor * np.cos(angle)) * unyt.m
    phi = 360 * rng.random(n) * unyt.deg
    z = rminor * np.sin(angle) * unyt.m
    return r, phi, z


def _relative_errors(field, reference, r, phi, z):
    """Return the maximum errors of B, psi and each Jacobian component relative
    to |B|, psi1 - psi0 and the Frobenius norm of the Jacobian of the reference
    field.

    The errors of B and psi are under "b" and "psi", and those of the Jacobian
    components under their names (qnt_mapping["bjac"]).
    """
    bcomp = ("br", "bphi", "bz")
    qnt = (*bcomp, "psi", *qnt_mapping["bjac"])
    val, ref = (
        dict(zip(qnt, evaluate(r, phi, z, 0*unyt.s, *qnt, bfield=f)))
        for f in (field, reference)
        )

    def jacobian(v):
        # (1/R) dB/dphi so that all components are in T/m
        return np.array([
            v[name].v / r.v if name.endswith("dphi") else v[name].v
            for name in qnt_mapping["bjac"]
            ])

    bnorm = np.sqrt(sum(ref[name].v**2 for name in bcomp))
    jac, jacref = jacobian(val), jacobian(ref)
    jnorm = np.sqrt(np.sum(jacref**2, axis=0))
    psinorm = np.abs(np.diff(reference.psilimits.v))[0]
    return {
        "b": max(np.max(np.abs((val[name] - ref[name]).v) / bnorm)
                 for name in bcomp),
        "psi": np.max(np.abs((val["psi"] - ref["psi"]).v)) / psinorm,
        **dict(zip(
            qnt_mapping["bjac"], np.max(np.abs(jac - jacref) / jnorm, axis=1),
            )),
        }


TOLERANCES = {
    "Spline2D": {"b": 2e-6, "psi": 3e-7, "jac": 2e-4},
    "Spline3D": {"b": 2e-6, "psi": 3e-7, "jac": 2e-4},
    "Linear3D": {"b": 1e-2, "psi": 2e-3, "jac": 5e-2},
    }
"""Tolerances of the errors returned by _relative_errors.

The "jac" tolerance applies to each Jacobian component. They are about three
times the errors of the fields.
"""


@pytest.mark.parametrize("variant", list(TOLERANCES))
def test_interpolate_tabulated(tabulated_fields, variant):
    """Test tabulated fields against the analytical field they were tabulated
    from.
    """
    reference, fields = tabulated_fields
    errors = _relative_errors(
        fields[variant], reference, *_plasma_points(1000, seed=2),
        )
    tolerances = TOLERANCES[variant]
    for name, error in errors.items():
        if variant == "Linear3D" and name in ("brdz", "bzdr"):
            # Bilinear psi has no second derivatives d2psi/dz2 and d2psi/dR2
            # (see test_interpolate_3d), which the poloidal field contributes
            # to these
            continue
        assert error < tolerances.get(name, tolerances["jac"]), (name, error)


def test_benchmark_linear3d_vs_spline(tabulated_fields, monkeypatch):
    """Compare accuracy and evaluation speed of the linear interpolation
    against splines using an analytical ITER-like circular field as reference.

    Only the C evaluation is timed, and on a single OpenMP thread, since that
    is how simulations use the field (each thread evaluates it for its own
    markers). The Python side of evaluate (unit conversions, allocating the
    results) runs on one thread and would otherwise dominate the timing.

    Run with ``pytest -s`` to see the results.
    """
    reference, tabulated = tabulated_fields
    spline2d = tabulated["Spline2D"]
    fields = {"Analytical": reference} | tabulated
    n = 100_000
    r, phi, z = _plasma_points(n, seed=3)

    ctimes = []
    interpolate_c = LIBASCOT.ascot_interpolate
    def timed_interpolate_c(*args):
        start = time.perf_counter()
        interpolate_c(*args)
        ctimes.append(time.perf_counter() - start)
    monkeypatch.setattr(LIBASCOT, "ascot_interpolate", timed_interpolate_c)

    print(
        f"\nBenchmark: iter-circular, {spline2d.rgrid.size}x"
        f"{spline2d.zgrid.size} (R,z) grid, {n} points, B + psi + Jacobian, "
        "C evaluation on one thread"
        )
    print(f"{'scheme':<11}{'err B':>11}{'err psi':>11}{'err jac':>11}"
          f"{'time [ms]':>11}{'ns/point':>10}")
    nthreads = LIBASCOT.omp_get_max_threads()
    LIBASCOT.omp_set_num_threads(1)
    try:
        for name, field in fields.items():
            errors = _relative_errors(field, reference, r, phi, z)
            jacerror = max(errors[name] for name in qnt_mapping["bjac"])
            ctimes.clear()
            for _ in range(3):
                evaluate(r, phi, z, 0*unyt.s, "br", "bphi", "bz", "psi",
                         *qnt_mapping["bjac"], bfield=field)
            best = min(ctimes)
            print(f"{name:<11}{errors['b']:>11.2e}{errors['psi']:>11.2e}"
                  f"{jacerror:>11.2e}{1e3*best:>11.1f}"
                  f"{1e9*best/n:>10.0f}")
    finally:
        LIBASCOT.omp_set_num_threads(nthreads)
