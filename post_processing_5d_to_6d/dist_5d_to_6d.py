import numpy as np
import matplotlib.pyplot as plt
from a5py import Ascot
import unyt
from a5py.ascot5io.dist import DistData
import math
import os
import bisect
import time
import sys


"""Helper functions to convert 5D distribution (r,phi,z,ppar,pperp) to 6D Cartesian distribution"""

def unit_vector(v):
    """Calculate the unit vector of a given vector."""
    if (v.ndim == 1):
        magnitude = np.linalg.norm(v)
        return v / magnitude
    elif (v.ndim == 2):
        magnitude = np.linalg.norm(v, axis = 1)
        return v / magnitude[:,np.newaxis]
    else:
        raise ValueError(f"Input vector must be either 1D or 2D, but got {v.ndim}D")

def cylindrical_to_cartesian(r, phi, z):
    phi_rad = np.radians(phi.value)*unyt.dimensionless
    x = r * np.cos(phi_rad)
    y = r * np.sin(phi_rad)
    return x, y, z

def cylindrical_to_spherical(r, phi, z):
    rho = np.sqrt(r**2 + z**2)
    theta = np.arccos(z/(np.sqrt(r**2 + z**2)))
    return rho, phi, theta

def spherical_to_cartesian(r, phi, theta):
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return x, y, z


def magnetic_field_cylindrical_to_cartesian(B_r, B_phi, B_z, phi):
    """
    Convert a magnetic field vector from cylindrical coordinates to Cartesian coordinates.

    Parameters:
        B_r (float or np.ndarray): Radial component of the magnetic field in cylindrical coordinates.
        B_phi (float or np.ndarray): Azimuthal component of the magnetic field in cylindrical coordinates.
        B_z (float or np.ndarray): Axial (z) component of the magnetic field in cylindrical coordinates.
        phi (float or np.ndarray): Azimuthal angle in radians.

    Returns:
        tuple: A tuple (B_x, B_y, B_z) representing the magnetic field components in Cartesian coordinates:
            - B_x (float or np.ndarray): x-component of the magnetic field.
            - B_y (float or np.ndarray): y-component of the magnetic field.
            - B_z (float or np.ndarray): z-component of the magnetic field (unchanged from input).
    """
    B_x = B_r * np.cos(phi) - B_phi * np.sin(phi)
    B_y = B_r * np.sin(phi) + B_phi * np.cos(phi)
    return B_x, B_y, B_z

# Not used in the code
def cartesian_to_spherical2(x, y, z):
    # Calculate spherical coordinates
    rho = np.sqrt(x**2 + y**2 + z**2)
    # [-pi,pi]
    phi = np.arctan2(y, x)
    # [0,pi]
    theta = np.arccos(z/rho)
    return rho, phi, theta

def cartesian_to_spherical(r, phi, z, Ax, Ay, Az):
    """Converting cartesian momentum vector to spherical using cylindrical position coordinates

    Args:
        r (np.array): radial coordinate
        phi (np.array): phi angle coordinate
        z (np.array): z coordinate
        Ax (np.array): Cartesian vector x-component
        Ay (np.array): Cartesian vector y-component
        Az (np.array): Cartesian vector z-component

    Returns:
        np.array: Spherical vector components
    """
    # [0,pi]
    theta = np.arctan2(r,z)
    Ar = Ax*np.sin(theta)*np.cos(phi) + Ay*np.sin(theta)*np.sin(phi) + Az*np.cos(theta)
    Atheta = Ax*np.cos(theta)*np.cos(phi) + Ay*np.cos(theta)*np.sin(phi) - Az*np.sin(theta)
    Aphi = -Ax*np.sin(phi) + Ay*np.cos(phi)    
    return Ar, Atheta, Aphi

def cartesian_to_cylindrical2(x, y, z):
    # Calculate spherical coordinates
    r = np.sqrt(x**2 + y**2)
    # [-pi,pi]
    phi = np.arctan2(y, x)
    return r, phi, z

def cartesian_to_cylindrical(Ax, Ay, Az, phi):
    Ar = Ax*np.cos(phi) + Ay*np.sin(phi)
    Aphi = -Ax*np.sin(phi)+Ay*np.cos(phi)
    return Ar, Aphi, Az

def bfield_momentum_to_MeV(ppar, pperp):
    momentum = np.sqrt(ppar**2 + pperp**2)# * unyt.kg * unyt.m / unyt.s
    e = momentum**2/(2*1.674927471e-27*unyt.kg)
    return e.to("MeV")


def calculate_3d_momentum(ppar, pperp, r, phi, z, b_vec, n_samples, pcoord):
    """
    Compute 3D momentum vectors from parallel and perpendicular components
    with respect to the local magnetic field direction.

    Parameters:
        ppar (ndarray): Parallel momentum (n_positions,)
        pperp (ndarray): Perpendicular momentum (n_positions,)
        r, phi, z (ndarray): Position arrays (n_positions,)
        b_vec (ndarray): Magnetic field vectors (3, n_positions)
        n_samples (int): Number of gyro-angle samples per position
        pcoord (str): Output coordinate system ('cartesian', 'spherical', or 'cylindrical')

    Returns:
        Momentum components in the selected coordinate system.
    """
    bhat = unit_vector(b_vec.T)
    # Arbitrary basis vector along the z-axis
    e1 = np.array([0, 0, 1])
    #  Basis vectors e1 and e2 perpendicular to the B_field
    e2 = np.cross(bhat, e1)
    e1 = unit_vector(e2)
    e2 = np.cross(bhat, e1)

    # Generate a random number for zeta in the range from 0 to 2*pi
    zeta = np.random.uniform(0, 2*np.pi, n_samples)

    # Compute cos(zeta)
    c = np.cos(zeta)
    s = np.sin(zeta)
    n_positions = b_vec.shape[1]
    perphat = np.empty((n_positions, n_samples, 3))
    e1 = e1[:,:,np.newaxis]
    e2 = e2[:,:,np.newaxis]
    # Does the sign matter (charge thing with alphas)?
    perphat[:,:,0] = s * e1[:,0,:] + c * e2[:,0,:];
    perphat[:,:,1] = s * e1[:,1,:] + c * e2[:,1,:];
    perphat[:,:,2] = s * e1[:,2,:] + c * e2[:,2,:];

    px = (ppar * bhat[:,0])[:, np.newaxis] + pperp * perphat[:,:,0];
    py = (ppar * bhat[:,1])[:, np.newaxis] + pperp * perphat[:,:,1];
    pz = (ppar * bhat[:,2])[:, np.newaxis] + pperp * perphat[:,:,2];
    if (pcoord == "cartesian"):
        return px, py, pz
    elif (pcoord == "spherical"):
        r_tiled = np.reshape(np.tile(r, n_samples),(r.shape[0], n_samples))
        phi_tiled = np.reshape(np.tile(phi, n_samples),(phi.shape[0], n_samples))
        z_tiled = np.reshape(np.tile(z, n_samples),(z.shape[0], n_samples))
        pr, pphi, ptheta = cartesian_to_spherical(r_tiled, phi_tiled, z_tiled, px, py, pz)
        return pr, pphi, ptheta
    elif (pcoord == "cylindrical"):
        prho, pphi, pz = cartesian_to_cylindrical(px, py, pz, np.reshape(np.tile(phi, n_samples),(phi.shape[0], n_samples)))
        return prho, pphi, pz
    else:
        raise ValueError(f"Unknown coordinate: {pcoord}")


def dist_5d_to_6d_momentumloop_general(dist_5d, n_samples, a5, pcoord, np1, np2, np3):
    print(f"Shape of 5D distribution {dist_5d.distribution().shape}, abscissae {dist_5d.abscissae}")
    r_edges = dist_5d.abscissa_edges("r")
    r_values = dist_5d.abscissa("r")
    phi_edges = dist_5d.abscissa_edges("phi")
    phi_values = dist_5d.abscissa("phi")
    z_edges = dist_5d.abscissa_edges("z")
    z_values = dist_5d.abscissa("z")
    ppar_values = dist_5d.abscissa("ppar")
    pperp_values = dist_5d.abscissa("pperp")
    ppar_edges = dist_5d.abscissa_edges("ppar")
    #pperp_edges = dist_5d.abscissa_edges("pperp")
    time_values = dist_5d.abscissa("time")

    # Do p edges need to be changed in some situations?
    p1_min, p1_max = ppar_edges.min()*2, ppar_edges.max()*2
    p2_min, p2_max = ppar_edges.min()*2, ppar_edges.max()*2
    p3_min, p3_max = ppar_edges.min()*2, ppar_edges.max()*2

    p1_edges = np.linspace(p1_min, p1_max, np1)# * ppar_edges.units
    p2_edges = np.linspace(p2_min, p2_max, np2)# * ppar_edges.units
    p3_edges = np.linspace(p3_min, p3_max, np3)# * ppar_edges.units
    hist_6d_arr = np.zeros((len(r_edges)-1, len(phi_edges)-1, len(z_edges)-1, len(p1_edges)-1, len(p2_edges)-1, len(p3_edges)-1), dtype=np.float64)

    r_grid, phi_grid, z_grid = np.meshgrid(r_values, phi_values, z_values, indexing='ij')
    # Let's ravel the grid
    r_grid = r_grid.ravel()
    phi_grid = phi_grid.ravel()
    z_grid = z_grid.ravel()
    ir_grid = np.digitize(r_grid, r_edges)-1
    iphi_grid = np.digitize(phi_grid, phi_edges)-1
    iz_grid = np.digitize(z_grid, z_edges)-1 

    a5.input_init(bfield=True)
    br, bphi, bz = a5.input_eval(r_grid, phi_grid, z_grid, time_values, 'br', 'bphi', 'bz')
    a5.input_free()
    bx, by, bz = magnetic_field_cylindrical_to_cartesian(br, bphi, bz, phi_grid)
    print(f"Iterations: {len(ppar_values)*len(pperp_values)}")
    itn = 0

    ir_grid = np.tile(ir_grid, n_samples)
    iphi_grid = np.tile(iphi_grid, n_samples)
    iz_grid = np.tile(iz_grid, n_samples)
    dist_5d  = dist_5d.integrate(copy=True,  charge=np.s_[:], time=np.s_[:] )
    for ippar, ppar in enumerate(ppar_values):
        for ipperp, pperp in enumerate(pperp_values):
            particles_weight = np.tile(dist_5d.histogram()[:,:,:,ippar,ipperp].ravel(), n_samples)
            p1, p2, p3 = calculate_3d_momentum(ppar, pperp,r_grid, phi_grid, z_grid, np.array([bx, by, bz]), n_samples, pcoord)
            ip1_grid = (np.digitize(p1, p1_edges)-1).ravel()
            ip2_grid = (np.digitize(p2, p2_edges)-1).ravel()
            ip3_grid = (np.digitize(p3, p3_edges)-1).ravel()
            np.add.at(hist_6d_arr, (ir_grid, iphi_grid, iz_grid, ip1_grid, ip2_grid, ip3_grid), (particles_weight/n_samples).v )
            if (itn % 500 == 0):
                    print(f"Iteration {itn}")
            itn += 1
    if (pcoord == "spherical"):
        dist_6d = DistData(hist_6d_arr, r=r_edges, phi = phi_edges, z=z_edges, pr = p1_edges, pphi = p2_edges, ptheta = p3_edges)
    elif (pcoord == "cartesian"):
        dist_6d = DistData(hist_6d_arr, r=r_edges, phi = phi_edges, z=z_edges, px = p1_edges, py = p2_edges, pz = p3_edges)
    elif (pcoord == "cylindrical"):
        dist_6d = DistData(hist_6d_arr, r=r_edges, phi = phi_edges, z=z_edges, prho = p1_edges, pphi = p2_edges, pz = p3_edges)
    else:
        raise ValueError("Unsupported coordinate system")

    print(f"6D distribution shape:, {dist_6d._distribution.shape}")
    return dist_6d



def process_momentum(dist_6d, plot_title, pcoord, save_path=None):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 5), sharey=True)

    if (pcoord == "cartesian"):
        px_dist_6d  = dist_6d.integrate(copy=True,  r=np.s_[:], z=np.s_[:], phi=np.s_[:], pz=np.s_[:], py=np.s_[:] )
        py_dist_6d  = dist_6d.integrate(copy=True,  r=np.s_[:], z=np.s_[:], phi=np.s_[:], px=np.s_[:], pz=np.s_[:] )
        pz_dist_6d  = dist_6d.integrate(copy=True,  r=np.s_[:], z=np.s_[:], phi=np.s_[:], px=np.s_[:], py=np.s_[:] )

        px_dist_6d.plot(axes=ax1)
        py_dist_6d.plot(axes=ax2)
        pz_dist_6d.plot(axes=ax3)
        fig.suptitle("Momentum distribution components", fontsize = 13)
        plt.tight_layout()
        plt.show()
    elif (pcoord == "cylindrical"):
        prho_dist_6d  = dist_6d.integrate(copy=True,  r=np.s_[:], z=np.s_[:], phi=np.s_[:], pz=np.s_[:], pphi=np.s_[:] )
        pphi_dist_6d  = dist_6d.integrate(copy=True,  r=np.s_[:], z=np.s_[:], phi=np.s_[:], prho=np.s_[:], pz=np.s_[:] )
        pz_dist_6d  = dist_6d.integrate(copy=True,  r=np.s_[:], z=np.s_[:], phi=np.s_[:], pphi=np.s_[:], prho=np.s_[:] )

        prho_dist_6d.plot(axes=ax1)
        pphi_dist_6d.plot(axes=ax2)
        pz_dist_6d.plot(axes=ax3)

    elif (pcoord == "spherical"):
        pr_dist_6d  = dist_6d.integrate(copy=True,  r=np.s_[:], z=np.s_[:], phi=np.s_[:], ptheta=np.s_[:], pphi=np.s_[:] )
        pphi_dist_6d  = dist_6d.integrate(copy=True,  r=np.s_[:], z=np.s_[:], phi=np.s_[:], pr=np.s_[:], ptheta=np.s_[:] )
        ptheta_dist_6d  = dist_6d.integrate(copy=True,  r=np.s_[:], z=np.s_[:], phi=np.s_[:], pphi=np.s_[:], pr=np.s_[:] )

        pr_dist_6d.plot(axes=ax1)
        pphi_dist_6d.plot(axes=ax2)
        ptheta_dist_6d.plot(axes=ax3)
    
    else:
        raise ValueError("Unknown coordinate")
    ax1.set_xlim(-1.3e-19, 1.3e-19)
    ax2.set_xlim(-1.3e-19, 1.3e-19)
    ax3.set_xlim(-1.3e-19, 1.3e-19)
    if save_path != None:
        fig.savefig(f"{save_path}/{plot_title}.png")
        print(f"Figure saved {save_path}/{plot_title}.png")



def generate_markers(a5, dist, n_markers, pcoordinate, bins, save_path=None):
    # Generate markers
    anum   = 1
    znum   = 0
    mass   = 1.0087*unyt.amu
    charge = 0*unyt.e
    mrk = a5.markergen.generate_6d(n_markers, mass, charge, anum, znum, dist, mode="prt", pcoord=pcoordinate)
    if (pcoordinate == "cartesian"):
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
        ax1.hist(mrk["px"].v, bins = bins)
        ax1.set_xlabel("px", fontsize = 13)
        ax2.hist(mrk["py"].v, bins = bins)
        ax2.set_xlabel("py",fontsize = 13)
        ax3.hist(mrk["pz"].v, bins = bins)
        ax3.set_xlabel("pz",fontsize = 13)
        fig.suptitle("Components of sampled marker momentum distributions", fontsize = 13)
        plt.tight_layout()
        plt.show()
        if save_path != None:
            fig.savefig(f"{save_path}/sampled_markers_cartesian.png")
            print(f"Figure saved {save_path}/sampled_markers_cartesian.png")   
    return mrk

def mrk_to_serpent(mrk, bins, save_path = None):
    px, py, pz = mrk["px"], mrk["py"], mrk["pz"]
    r, phi, z = mrk["r"], mrk["phi"], mrk["z"]
    x, y, z = cylindrical_to_cartesian(r, phi, z)
    
    m_n = 1.674927471e-27 * unyt.kg
    ekin = ((px**2 + py**2 + pz**2) / (2 * m_n)).to("MeV")
    print(f"Mean energy {np.mean(ekin.v):.3f} MeV, std {np.std(ekin.v):.3f} MeV")
    dir_vec = unit_vector(np.column_stack((px, py, pz)))
    fig = plt.figure(figsize=(18, 8))  # Increase height to fit second row

    # First row
    ax1 = fig.add_subplot(2, 3, 5)
    ax1.hist(ekin.v, bins=bins)
    ax1.set_xlabel("Energy [MeV]",fontsize = 13)
    plt.axvline(x=14.1, c="black")

    ax2 = fig.add_subplot(2, 3, 1)
    ax2.hist(dir_vec.v[:, 0], bins=bins)
    ax2.set_xlabel("u",fontsize = 13)

    ax3 = fig.add_subplot(2, 3, 2, sharey=ax2)
    ax3.hist(dir_vec.v[:, 1], bins=bins)
    ax3.set_xlabel("v",fontsize = 13)

    ax4 = fig.add_subplot(2, 3, 3, sharey=ax2)
    ax4.hist(dir_vec.v[:, 2], bins=bins)
    ax4.set_xlabel("w", fontsize = 13)

    # Second row – example additional plot
    ax5 = fig.add_subplot(2, 3, 4)
    ax5.scatter(x, z)  # Replace weight.v with the desired data
    ax5.set_xlabel("x",  fontsize = 13)
    ax5.set_ylabel("z",  fontsize = 13)
    fig.suptitle("Serpent marker distributions", fontsize = 13)
    
    hist1, _, _ = np.histogram2d(r, z)
    ax6 = fig.add_subplot(2, 3, 6)
    im1 = ax6.imshow(hist1.T, origin="lower", extent=[r.min(), r.max(), z.min(), z.max()], cmap="hot", aspect="auto")
    plt.colorbar(im1, ax=ax6)
    ax6.set_xlabel("x",  fontsize = 13)
    ax6.set_ylabel("z",  fontsize = 13)
    plt.tight_layout()
    plt.show()
    
    isotropic_dir = isotropic_directions_rejection(len(px))
    plot_direction_vectors(isotropic_dir, dir_vec, bins = 40)
    
    
    if save_path != None:
        fig.savefig(f"{save_path}/sampled_serpent_markers.png")
        print(f"Figure saved {save_path}/sampled_serpent_markers.png")
    return ekin, dir_vec,x,y,z

# TODO: Check if this isotropic sampling is ok
def sample_isotropic_directions(n_samples):
    phi = np.random.uniform(0, 2 * np.pi, n_samples)       # Azimuthal angle
    cos_theta = np.random.uniform(-1, 1, n_samples)        # Cosine of polar angle
    sin_theta = np.sqrt(1 - cos_theta**2)

    x = sin_theta * np.cos(phi)
    y = sin_theta * np.sin(phi)
    z = cos_theta

    directions = np.vstack((x, y, z)).T  # Shape (n_samples, 3)
    
    return directions



def plot_direction_vectors(analytical_directions, sampled_directions, bins=50):
    fig, axes = plt.subplots(3, 2, figsize=(10, 8), sharey='row')

    labels = ['u', 'v', 'w']

    for i in range(3):
        # Analytical in left column
        axes[i, 0].hist(analytical_directions[:, i], bins=bins, alpha=0.7, color='tab:blue')
        axes[i, 0].set_xlabel(labels[i], fontsize=12)
        if i == 0:
            axes[i, 0].set_title("Analytical", fontsize=13)

        # Sampled in right column
        axes[i, 1].hist(sampled_directions[:, i], bins=bins, alpha=0.7, color='tab:orange')
        axes[i, 1].set_xlabel(labels[i], fontsize=12)
        if i == 0:
            axes[i, 1].set_title("Sampled from 6D dist.", fontsize=13)
    plt.suptitle("Components of unit direction vectors")
    plt.tight_layout()
    plt.show()


def isotropic_directions_rejection(n, rng=np.random.default_rng()):
    """
    Generate `n` isotropic direction vectors uniformly distributed over the unit sphere.

    Parameters:
        n (int): Number of direction vectors to generate.
        rng: numpy random generator instance (default: np.random.default_rng())

    Returns:
        directions (np.ndarray): Array of shape (n, 3) with unit direction vectors.
    """
    directions = np.empty((n, 3))
    i = 0
    while i < n:
        rand1 = 2.0 * rng.random() - 1.0
        rand2 = 2.0 * rng.random() - 1.0
        C1 = rand1**2 + rand2**2
        if C1 > 1.0:
            continue

        rand3 = 2.0 * rng.random() - 1.0
        C2 = np.sqrt(1.0 - rand3**2)

        u = C2 * (rand1**2 - rand2**2) / C1
        v = C2 * 2.0 * rand1 * rand2 / C1
        w = rand3

        directions[i] = [u, v, w]
        i += 1
    return directions

# TODO: Add possibility to make afsi run before 6D conversion
def iter_analytical_5D(a5, rmin,rmax,nr, phimin,phimax, nphi, zmin, zmax, nz, nppar, npperp, nsamples = 1000, pparmin=-1.1e-19, pparmax = 1.1e-19, pperpmin = 0, pperpmax = 1.1e-19, desc = None):
    # Feed the parameters correctly
    r = np.linspace(rmin, rmax, nr)
    z = np.linspace(zmin, zmax, nz)
    phi = np.linspace(phimin, phimax, nphi)
    ppar = np.linspace(pparmin, pparmax, nppar)
    pperp = np.linspace(pperpmin,pperpmax, npperp)
    a5.afsi.thermal(
    "DT_He4n", nmc=nsamples, r=r, z=z,phi=phi,
    ppar1=ppar, pperp1=pperp, ppar2=ppar, pperp2=pperp
    )
    dist_5d = a5.data.active.getdist("prod2")
    if desc == None:
        desc = f"R{nr}_Z{nz}_ppar{nppar}_pperp{npperp}"
        a5.data.active.set_desc(desc)
    return dist_5d, a5
    

