import os 
# Change to same folder as dist_5d_to_6d.py
os.chdir("ascot5/post_processing_5d_to_6d")

from dist_5d_to_6d import *


def main():
    # Specify ascot file path
    input_file = "source_iter_test/ascot.h5" 
    n_samples = 1000 
    run_title = "iter_test" 
    a5 = Ascot(input_file)
    np1 = 51
    np2 = 51
    np3 = 51
    pbins = 40
    plot_save_path = None # Specify save path for plots, if you want to save them. 
    n_markers = 10000
    dist_5d = a5.data.active.getdist("prod2")
    
    #TODO: Adding afsi run in this script 
    # dist_5d, a5 = iter_analytical_5D()
    #coordinates = ["cartesian"]#, "cylindrical", "spherical"]
    coord = "cartesian"
    
    dist_6d_file = f"hist6d_{run_title}_{coord}.npy"
    start_time = time.time()
    print(f"N samples {n_samples}")
    dist_6d = dist_5d_to_6d_momentumloop_general(dist_5d, n_samples, a5, pcoord=coord, np1 = np1, np2=np2, np3 = np3)
    # Saving dist histogram is commented out
    #print(f"Saving distribution to a file {dist_6d_file}")
    #np.save(dist_6d_file, dist_6d.histogram())
    print("Processing distribution data")
    plot_title = f"{run_title}_{coord}"
    process_momentum(dist_6d, dist_5d, plot_title, coord, plot_save_path)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time/(60):.2f} minutes")
    print(f"6D distribution shape: {dist_6d._distribution.shape}, abscissae {dist_6d.abscissae}")
    print("Distribution saved")
    markers = generate_markers(a5, dist_6d,n_markers, "cartesian", pbins, plot_save_path)
    serp_markers = mrk_to_serpent(markers, pbins, plot_save_path)
     
if __name__ == "__main__":
    main()

print(process_momentum)