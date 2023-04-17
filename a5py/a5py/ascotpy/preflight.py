"""
Checks to help user to ensure the inputs are ok before running the simulation.
"""

import numpy as np
import matplotlib.pyplot as plt
from a5py.ascotpy import Ascotpy
from a5py.misc import openfigureifnoaxes

def check_inputs_present(ascotpy):
    """
    Check required inputs are present for this run.
    """

    # Determine what we are simulating.
    runtype = "ascot"

    msg = []
    if runtype == "ascot":
        if not hasattr(ascotpy.hdf5, "bfield"):
            msg += ["Error: bfield data is missing"]
        if not hasattr(ascotpy.hdf5, "efield"):
            msg += ["Error: efield data is missing"]
        if not hasattr(ascotpy.hdf5, "plasma"):
            msg += ["Error: plasma data is missing"]
        if not hasattr(ascotpy.hdf5, "marker"):
            msg += ["Error: marker data is missing"]
        if not hasattr(ascotpy.hdf5, "options"):
            msg += ["Error: options are missing"]
        if not hasattr(ascotpy.hdf5, "boozer"):
            msg += ["Error: boozer data is missing"]
        if not hasattr(ascotpy.hdf5, "wall"):
            msg += ["Error: wall data is missing"]
        if not hasattr(ascotpy.hdf5, "mhd"):
            msg += ["Error: mhd data is missing"]


    return msg

def check_options_consistent(ascotpy):
    """
    Check that options are consistent with inputs and each other.
    """
    msg = []
    opt = ascotpy.hdf5.options.active

    msg0 = opt.validatevalues()
    if len(msg0) > 0:
        msg += ["Error: following options had invalid parameters:"]
        msg += msg0

    opt = opt.read()
    if opt["ENABLE_MHD"] == 1 and \
       ascotpy.hdf5.bfield.active.get_type() == "B_STS":
        msg += ["Error: cannot enable MHD for stellarators"]

    if opt["SIM_MODE"] in [1,2,3] and \
       ascotpy.hdf5.marker.active.get_type() == "mrk_fl":
        msg += ["Error: SIM_MODE invalid for field-line markers"]

    # Make rough estimates on how much memory is consumed by diagnostics
    high_memory_consumption = 1e9 # In bits
    # Orbit = Npoint * Nfields (~10) * Nmrk * 8 bit
    orb_mem = opt["ORBITWRITE_NPOINT"] * 10 \
              * ascotpy.hdf5.marker.active.read()["n"] * 8
    # Distributions = NR * Nz * ... * 8 bit
    rzp = opt["DIST_NBIN_R"] * opt["DIST_NBIN_Z"] * opt["DIST_NBIN_PHI"]
    rtp = opt["DIST_NBIN_RHO"] * opt["DIST_NBIN_THETA"] * opt["DIST_NBIN_PHI"]
    p2d = opt["DIST_NBIN_PPA"] * opt["DIST_NBIN_PPE"]
    p3d = opt["DIST_NBIN_PR"] * opt["DIST_NBIN_PZ"] * opt["DIST_NBIN_PPHI"]

    if opt["ENABLE_ORBITWRITE"] == 1 and orb_mem > high_memory_consumption:
        msg += ["Warning: orbit diagnostic memory consumption high (~" +
                str(int(orb_mem / 1e9)) + "Gb)"]

    if opt["ENABLE_DIST_5D"] == 1 and rzp * p2d * 8 > high_memory_consumption:
        msg += ["Warning: 5D distribution memory consumption high (~" +
                str(int(rzp * p2d * 8 / 1e9)) + "Gb)"]

    if opt["ENABLE_DIST_6D"] == 1 and rzp * p3d * 8 > high_memory_consumption:
        msg += ["Warning: 6D distribution memory consumption high (~" +
                str(int(rzp * p3d * 8 / 1e9)) + "Gb)"]

    if opt["ENABLE_DIST_RHO5D"] == 1 and rtp * p2d * 8 > high_memory_consumption:
        msg += ["Warning: rho5D distribution memory consumption high (~" +
                str(int(rtp * p2d * 8 / 1e9)) + "Gb)"]

    if opt["ENABLE_DIST_RHO6D"] == 1 and rtp * p3d * 8 > high_memory_consumption:
        msg += ["Warning: rho6D distribution memory consumption high (~" +
                str(int(rtp * p3d * 8 / 1e9)) + "Gb)"]

    return msg


def check_bfield_psi0(ascotpy):
    """
    Checks whether psi0 given in input is actually extreme value.

    Because psi is interpolated with splines, there might be numerical error
    that causes psi (near the axis) to have more extreme value than psi0
    which is given in input. This leads to imaginary rho and termination
    of the simulation if marker ends up there.

    This check uses Monte Carlo method to 1. Draw phi 2. Evaluate axis (R,z)
    3. Draw random (R,z) coordinates within 10 cm of the axis. 4. Evaluate psi
    at that point and compare to psi0. Process repeats N times. Check passes
    if all evaluations are valid.
    """

    data = ascotpy.hdf5.bfield.active.read()
    if "psi0" not in data:
        return []

    psi0 = data["psi0"]
    psi1 = data["psi1"]
    psi0 = -7.1

    N     = 10000
    phi   = np.random.rand(N,) * 360
    theta = np.random.rand(N,) * 2 * np.pi

    axis = ascotpy.evaluate(1, phi, 0, 0, "axis")
    z0 = axis["axisz"]
    r0 = axis["axisr"]

    R = 0.1 # 10 cm
    r = R * np.cos(theta) + r0
    z = R * np.sin(theta) + z0
    psi = ascotpy.evaluate(r, phi, z, 0, "psi")

    if psi0 < psi1 and any(psi < psi0):
        return ["Error: psi0 = %.2e but we found near axis that psi = %.2e" \
                % (psi0, np.amin(psi)) ]
    elif psi0 > psi1 and any(psi > psi0):
        return ["Error: psi0 = %.2e but we found near axis that psi = %.2e" \
                % (psi0, np.amax(psi)) ]

    return []


def plot_top_view(ascotpy, axes=None):
    """
    Plot top view of the machine showing Ip, Bphi, and markers.

    Assumes bfield is initialized in ascotpy.
    """
    showfig = False
    if axes == None:
        fig  = plt.figure()
        axes = fig.add_subplot(1,1,1)
        showfig = True

    axis = ascotpy.evaluate(1, 0, 0, 0, "axis")
    z0 = axis["axisz"]
    r0 = axis["axisr"]

    rmin = r0-r0/2
    rmax = r0+r0/2
    dphi = 10 * np.pi/180 # So that b and j quivers dont overlap

    r   = np.linspace(rmin, rmax, 10)
    phi = np.linspace(0, 360, 18, endpoint=False) * np.pi/180
    phi, r = np.meshgrid(phi, r)
    r   = r.ravel()
    phi = phi.ravel()

    br   = np.squeeze(ascotpy.evaluate(r, phi, z0, 0, "br"))
    bphi = np.squeeze(ascotpy.evaluate(r, phi, z0, 0, "bphi"))
    jr   = np.squeeze(ascotpy.evaluate(r, phi + dphi, z0, 0, "jr"))
    jphi = np.squeeze(ascotpy.evaluate(r, phi + dphi, z0, 0, "jphi"))

    x  = np.cos(phi) * r
    y  = np.sin(phi) * r
    bx = np.cos(phi) * br - np.sin(phi) * bphi
    by = np.sin(phi) * br + np.cos(phi) * bphi
    bnorm = np.sqrt(bx**2 + by**2)

    xj = np.cos(phi+dphi) * r
    yj = np.sin(phi+dphi) * r
    jx = np.cos(phi+dphi) * jr - np.sin(phi+dphi) * jphi
    jy = np.sin(phi+dphi) * jr + np.cos(phi+dphi) * jphi
    jnorm = np.sqrt(jx**2 + jy**2)

    axes.quiver(x,y, bx/bnorm, by/bnorm,
                color="blue", scale=40)
    axes.quiver(xj,yj, jx/jnorm, jy/jnorm,
                color="red", scale=40)

    axes.set_aspect("equal", adjustable="box")

    marker = ascotpy.hdf5.marker.active.read()
    x  = np.cos(marker["phi"] * np.pi/180) * marker["r"]
    y  = np.sin(marker["phi"] * np.pi/180) * marker["r"]

    axes.scatter(x,y, s=1, c="black", zorder=-2)

    axes.legend(("Magnetic field", "Plasma current", "Markers"))

    if showfig:
        plt.show()


@openfigureifnoaxes
def plot_energypitch(ascotpy, axes=None):
    """
    Plot marker energy-pitch histogram
    """
    ascotpy.hdf5.marker.active.plot_hist_energypitch(ascotpy, axes=axes)


@openfigureifnoaxes
def plot_rhophi(ascotpy, axes=None):
    """
    Plot marker rho-phi histogram
    """
    ascotpy.hdf5.marker.active.plot_hist_rhophi(ascotpy, axes=axes)


if __name__ == '__main__':
    ascotpy = Ascotpy("ascot.h5")
    ascotpy.init(bfield=True)

    msg = []
    msg += check_inputs_present(ascotpy)
    msg += check_options_consistent(ascotpy)
    msg += check_bfield_psi0(ascotpy)

    for s in msg:
        print(s)

    if len(msg) == 0:
        print("Preflight checks ok!")

    fig  = plt.figure()
    plot_top_view(ascotpy,    axes=fig.subplot(2,2,[1,3]))
    plot_rhophi(ascotpy,      axes=fig.subplot(2,2,2))
    plot_energypitch(ascotpy, axes=fig.subplot(2,2,4))
    plt.show()

    ascotpy.free(bfield=True)
