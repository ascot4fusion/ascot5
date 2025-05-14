import numpy as np
import matplotlib.pyplot as plt
from serpent_markers_functions import *


def main():

    # Specify ascot file path
    input_file = "/home/giacomo/ascot5/post_processing_5d_to_6d/ascot.h5"
    output_path = "/home/giacomo/ascot5/post_processing_5d_to_6d/neutron_dist.txt"
    a5 = Ascot(input_file)

    #Set geometric parameters to run afsi
    rmin = 5.2
    rmax = 7.2
    nr = 51
    phimin = -0.1
    phimax = 0.1
    nphi = 2
    zmin = -1.0
    zmax = 1.0
    nz = 51
    n_samples = 100000

    # dist_5d = spatial_dist(a5, rmin,rmax,nr, phimin,phimax, nphi, zmin, zmax, nz, n_samples)
    # dist_5d = a5.data.R21Z21PPAR2PPERP2.getdist("prod2")
    dist_5d = a5.data.R51Z51PPAR2PPERP2.getdist("prod2")
    # dist_5d = a5.data.active.getdist("prod2")

    samples_per_bin = samples_in_bin(dist_5d, n_samples)

    output_file = open(output_path, "w")
    fill_output_file(output_file, dist_5d, a5, samples_per_bin)
    output_file.close()

    plot_markers(output_path)

    return

main()