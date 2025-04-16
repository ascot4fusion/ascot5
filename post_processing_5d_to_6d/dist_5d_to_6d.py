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


def unit_vector(v):
    """Calculate the unit vector of a given vector."""
    if (v.ndim == 1):
        magnitude = np.linalg.norm(v)
        return v / magnitude
    elif (v.ndim == 2):
        magnitude = np.linalg.norm(v, axis = 1)
        return v / magnitude[:,np.newaxis]
    elif (v.ndim == 4):
        magnitude = np.linalg.norm(v, axis = 0)
        return v / magnitude
    else:
        raise ValueError("Vector must be 1D or 2D")

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
    B_x = B_r * np.cos(phi) - B_phi * np.sin(phi)
    B_y = B_r * np.sin(phi) + B_phi * np.cos(phi)
    return B_x, B_y, B_z

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


def calculate_3d_momentum_grid_samples(ppar, pperp, r, phi, z, b_vec, n_samples, pcoord):
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
        #e_kin = (p**2/(2*1.674927471e-27*unyt.kg)).to("MeV")
        return pr, pphi, ptheta
    elif (pcoord == "cylindrical"):
        prho, pphi, pz = cartesian_to_cylindrical(px, py, pz, np.reshape(np.tile(phi, n_samples),(phi.shape[0], n_samples)))
        return prho, pphi, pz
    else:
        raise ValueError("Unknown coordinate")


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
    p1_edges = np.linspace(-2e-19, 2e-19, np1)*ppar_edges.units
    p2_edges = np.linspace(-2e-19, 2e-19, np2)*ppar_edges.units
    p3_edges = np.linspace(-2e-19, 2e-19, np3)*ppar_edges.units
    hist_6d_arr = np.empty((len(r_edges)-1, len(phi_edges)-1, len(z_edges)-1, len(p1_edges)-1, len(p2_edges)-1, len(p3_edges)-1), dtype=np.float64)

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
    for ippar, ppar in enumerate(ppar_values):
        for ipperp, pperp in enumerate(pperp_values):
            particles_weight = np.tile(dist_5d.histogram()[:,:,:,ippar,ipperp,0,0].ravel(), n_samples)
            p1, p2, p3 = calculate_3d_momentum_grid_samples(ppar, pperp,r_grid, phi_grid, z_grid, np.array([bx, by, bz]), n_samples, pcoord)
            ip1_grid = (np.digitize(p1, p1_edges)-1).ravel()
            ip2_grid = (np.digitize(p2, p2_edges)-1).ravel()
            ip3_grid = (np.digitize(p3, p3_edges)-1).ravel()
            np.add.at(hist_6d_arr, (ir_grid, iphi_grid, iz_grid, ip1_grid, ip2_grid, ip3_grid), (particles_weight/n_samples).v )
            if (itn % 100 == 0):
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




def process_momentum(dist_6d, dist_5d, plot_title, pcoord):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 5), sharey=True)

    if (pcoord == "cartesian"):
        px_dist_6d  = dist_6d.integrate(copy=True,  r=np.s_[:], z=np.s_[:], phi=np.s_[:], pz=np.s_[:], py=np.s_[:] )
        py_dist_6d  = dist_6d.integrate(copy=True,  r=np.s_[:], z=np.s_[:], phi=np.s_[:], px=np.s_[:], pz=np.s_[:] )
        pz_dist_6d  = dist_6d.integrate(copy=True,  r=np.s_[:], z=np.s_[:], phi=np.s_[:], px=np.s_[:], py=np.s_[:] )

        px_dist_6d.plot(axes=ax1)
        py_dist_6d.plot(axes=ax2)
        pz_dist_6d.plot(axes=ax3)

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


    ax1.set_xlim(-1.5e-19, 1.5e-19)
    ax2.set_xlim(-1.5e-19, 1.5e-19)
    ax3.set_xlim(-1.5e-19, 1.5e-19)
    plt.savefig(f"plots/updated_momentum/{plot_title}.png")
    print(f"Figure saved plots/updated_momentum/{plot_title}.png")



def generate_markers(a5, dist, n_markers, pcoordinate, bins):
    # Generate markers
    anum   = 1
    znum   = 0
    mass   = 1.0087*unyt.amu
    charge = 0*unyt.e
    mrk = a5.markergen.generate_6d(n_markers, mass, charge, anum, znum, dist, mode="prt", pcoord=pcoordinate)
    if (pcoordinate == "cartesian"):
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
        ax1.hist(mrk["px"].v, bins = bins)
        ax1.set_xlabel("px")
        ax2.hist(mrk["py"].v, bins = bins)
        ax2.set_xlabel("py")
        ax3.hist(mrk["pz"].v, bins = bins)
        ax3.set_xlabel("pz")
        
    return mrk

def mrk_to_serpent(mrk, bins):
    px = mrk["px"]
    py = mrk["py"]
    pz = mrk["pz"]
    r = mrk["r"]
    phi = mrk["phi"]
    z = mrk["z"]
    x,y,z = cylindrical_to_cartesian(r, phi, z)
    ekin = ((px**2+py**2+pz**2)/(2*1.674927471e-27*unyt.kg)).to("MeV")
    dir_vec = unit_vector(np.column_stack((px, py, pz)))
    fig = plt.figure(figsize = (18,4))
    ax1 = fig.add_subplot(1,4,1)
    ax1.hist(ekin.v, bins = bins)
    ax1.set_xlabel("Energy [MeV]")
    plt.axvline(x = 14.1, c = "black")

    ax2 = fig.add_subplot(1, 4, 2)
    ax2.hist(dir_vec.v[:, 0], bins=bins)
    ax2.set_xlabel("x-dir")

    # Third subplot sharing y-axis with ax2
    ax3 = fig.add_subplot(1, 4, 3, sharey=ax2)
    ax3.hist(dir_vec.v[:, 1], bins=bins)
    ax3.set_xlabel("y-dir")

    # Fourth subplot sharing y-axis with ax2
    ax4 = fig.add_subplot(1, 4, 4, sharey=ax2)
    ax4.hist(dir_vec.v[:, 2], bins=bins)
    ax4.set_xlabel("z-dir")

    plt.show()
    return ekin, dir_vec,x,y,z


def main(args):
    input_file = args[1]
    n_samples = int(args[2])
    run_title = args[3]
    a5 = Ascot(input_file)

    dist_5d = a5.data.active.getdist("prod2")
    for coord in ["cartesian", "cylindrical", "spherical"]:
        dist_6d_file = f"hist6d_{run_title}_{coord}.npy"
        start_time = time.time()
        print(f"N samples {n_samples}")
        dist_6d = dist_5d_to_6d_momentumloop_general(dist_5d, n_samples, a5, pcoord=coord, np1 = 101, np2=101, np3 = 101)
        print(f"Saving distribution to a file {dist_6d_file}")
        np.save(dist_6d_file, dist_6d.histogram())
        print("Processing distribution data")
        plot_title = f"{run_title}_{coord}"
        process_momentum(dist_6d, dist_5d, plot_title, coord)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time: {elapsed_time/(60):.2f} minutes")
        print(f"6D distribution shape: {dist_6d._distribution.shape}, abscissae {dist_6d.abscissae}")
        print("Distribution saved")

if __name__ == "__main__":
    main(sys.argv)


