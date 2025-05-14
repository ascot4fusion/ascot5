import numpy as np
import matplotlib.pyplot as plt
from a5py import Ascot
import pandas as pd

def fill_output_file(file, dist_5d, a5, samples_per_bin):

    r_edges = dist_5d.abscissa_edges("r")
    z_edges = dist_5d.abscissa_edges("z")
    phi_edges = dist_5d.abscissa_edges("phi")

    it = 0
    for ir in range(len(r_edges)-1):
        # Get the r bin edges
        rmin_bin = r_edges[ir]
        rmax_bin = r_edges[ir+1]
        for iphi in range(len(phi_edges)-1):
            # Get the phi bin edges
            phimin_bin = phi_edges[iphi]
            phimax_bin = phi_edges[iphi+1]
            for iz in range(len(z_edges)-1):
                # Get the number of samples for each bin
                n_samples_bin = samples_per_bin[ir,iphi,iz]
                # Get the z bin edges
                zmin_bin = z_edges[iz]
                zmax_bin = z_edges[iz+1]
                # call afsi for each bin
                prod2  = iter_6D_markers(a5, rmin_bin,rmax_bin, phimin_bin, phimax_bin, zmin_bin, zmax_bin, nsamples=n_samples_bin)
                it = it+1
                print(it)
                # Save the results
                np.savetxt(file, prod2, fmt='%1.6f')  # write the array
                file.write('\n')
    
    return

def samples_in_bin(dist_5d, n_samples):

    dist_5d.integrate(time=np.s_[:], charge=np.s_[:])
    spat_dist = dist_5d.integrate(copy=True, ppar=np.s_[:], pperp=np.s_[:])
    spat_hist = spat_dist.histogram()
    prob_3d = spat_hist/np.sum(spat_hist)
    prob_flat = prob_3d.flatten()

    samples_per_bin_flat = np.random.multinomial(n_samples, prob_flat)
    samples_per_bin = samples_per_bin_flat.reshape(prob_3d.shape)

    return samples_per_bin

def spatial_dist(a5, rmin,rmax,nr, phimin,phimax, nphi, zmin, zmax, nz, nsamples, pparmin=-1.1e-19, pparmax = 1.1e-19, pperpmin = 0, pperpmax = 1.1e-19, desc = None):
    # Feed the parameters correctly
    r = np.linspace(rmin, rmax, nr)
    z = np.linspace(zmin, zmax, nz)
    phi = np.linspace(phimin, phimax, nphi)
    ppar = np.linspace(pparmin, pparmax, 2)
    pperp = np.linspace(pperpmin,pperpmax, 2)
    a5.afsi.thermal(
    "DT_He4n", nmc=nsamples, r=r, z=z,phi=phi,
    ppar1=ppar, pperp1=pperp, ppar2=ppar, pperp2=pperp)
    dist_5d = a5.data.active.getdist("prod2")
    if desc == None:
            desc = f"R{nr}_Z{nz}_ppar{2}_pperp{2}"
            a5.data.active.set_desc(desc)
    return dist_5d

def iter_6D_markers(a5, rmin,rmax, phimin,phimax, zmin, zmax, nsamples):
    # Feed the parameters correctly
    r = np.linspace(rmin, rmax, 2)
    z = np.linspace(zmin, zmax, 2)
    phi = np.linspace(phimin, phimax, 2)

    prod2 = a5.afsi.thermal_3D(
    "DT_He4n", nmc=nsamples, r=r, z=z,phi=phi)
    return prod2

def plot_markers(output_path):
     
    df = pd.read_csv(output_path, sep=" ", header=None)
    r = df[0].to_numpy()
    phi = df[1].to_numpy()
    z = df[2].to_numpy()
    u = df[3].to_numpy()
    v = df[4].to_numpy()
    w = df[5].to_numpy()
    en = df[6].to_numpy()

    plt.subplot(2,3,1)
    plt.hist(r, bins=50)
    plt.title("r")
    plt.subplot(2,3,2)
    plt.hist(phi, bins=50)
    plt.title("phi")
    plt.subplot(2,3,3)
    plt.hist(z, bins=50)
    plt.title("z")
    plt.subplot(2,3,4)
    plt.hist(u, bins=50)
    plt.title("u")
    plt.subplot(2,3,5)
    plt.hist(v, bins=50)
    plt.title("v")
    plt.subplot(2,3,6)
    plt.hist(w, bins=50)
    plt.title("w")
    plt.show()

    plt.hist(en, bins=50)
    plt.axvline(x=14.1, color='black')
    plt.title("en")
    plt.xlabel("Energy [MeV]")
    plt.show()

    return