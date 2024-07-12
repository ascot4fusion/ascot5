/**
 * @file bbnbi5_main.c
 * @brief BBNBI5 main program
 */
#include <getopt.h>
#include <math.h>
#ifdef MPI
  #include <mpi.h>
#endif
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ascot5.h"
#include "gitver.h"
#include "math.h"
#include "hdf5_interface.h"
#include "hdf5io/hdf5_helpers.h"
#include "print.h"
#include "simulate.h"
#include "particle.h"
#include "B_field.h"
#include "plasma.h"
#include "neutral.h"
#include "wall.h"
#include "asigma.h"
#include "nbi.h"
#include "diag.h"
#include "bbnbi5.h"

/**
 * @brief Main function for BBNBI5
 *
 * @param  argc argument count of the command line arguments
 * @param  argv argument vector of the command line arguments
 *
 * @return Zero if simulation was completed
 */
int main(int argc, char** argv) {

    /* Read and parse command line arguments */
    int nprt;    /* Number of markers to be generated in total */
    real t1, t2; /* Markers are initialized in this time-spawn */
    sim_offload_data sim;
    if(bbnbi_read_arguments(argc, argv, &sim, &nprt, &t1, &t2) != 0) {
        abort();
        return 1;
    }

    /* Init MPI (but BBNBI 5 does not yet support MPI) */
    sim.mpi_rank = 0;
    sim.mpi_root = 0;
    print_out0(VERBOSE_MINIMAL, sim.mpi_rank, sim.mpi_root, "BBNBI5\n");
    print_out0(VERBOSE_MINIMAL, sim.mpi_rank, sim.mpi_root,
               "Tag %s\nBranch %s\n\n", GIT_VERSION, GIT_BRANCH);

    /* Read data needed for bbnbi simulation */
    real* nbi_offload_array;
    real* B_offload_array;
    real* plasma_offload_array;
    real* neutral_offload_array;
    real* wall_offload_array;
    int*  wall_int_offload_array;
    real* asigma_offload_array;
    real* diag_offload_array;
    if( hdf5_interface_read_input(&sim, hdf5_input_bfield | hdf5_input_plasma |
                                  hdf5_input_neutral | hdf5_input_wall |
                                  hdf5_input_asigma | hdf5_input_nbi |
                                  hdf5_input_options,
                                  &B_offload_array, NULL, &plasma_offload_array,
                                  &neutral_offload_array, &wall_offload_array,
                                  &wall_int_offload_array, NULL, NULL,
                                  &asigma_offload_array, &nbi_offload_array,
                                  NULL, NULL) ) {
        print_out0(VERBOSE_MINIMAL, sim.mpi_rank, sim.mpi_root,
                   "Input initialization failed\n");
        abort();
        return 1;
    }

    /* Disable diagnostics that are not supported */
    sim.diag_offload_data.diagorb_collect   = 0;
    sim.diag_offload_data.diagtrcof_collect = 0;

    /* Initialize diagnostics */
    if( diag_init_offload(&sim.diag_offload_data, &diag_offload_array, 0) ) {
        print_out0(VERBOSE_MINIMAL, sim.mpi_rank, sim.mpi_root,
                       "\nFailed to initialize diagnostics.\n"
                       "See stderr for details.\n");
            abort();
            return 1;
    }
    real diag_offload_array_size = sim.diag_offload_data.offload_array_length
        * sizeof(real) / (1024.0*1024.0);
    print_out0(VERBOSE_IO, sim.mpi_rank, sim.mpi_root,
               "Initialized diagnostics, %.1f MB.\n", diag_offload_array_size);
    simulate_init_offload(&sim);

    /* QID for this run */
    char qid[11];
    hdf5_generate_qid(qid);

    /* Write bbnbi run group to HDF5 */
    if(sim.mpi_rank == sim.mpi_root) {
        print_out0(VERBOSE_IO, sim.mpi_rank, sim.mpi_root,
                   "\nPreparing output.\n");
        if( hdf5_interface_init_results(&sim, qid, "bbnbi") ) {
            print_out0(VERBOSE_MINIMAL, sim.mpi_rank, sim.mpi_root,
                       "\nInitializing output failed.\n"
                       "See stderr for details.\n");
            /* Free offload data and terminate */
            abort();
            return 1;
        }
        strcpy(sim.qid, qid);
    }

    /* Inject markers from the injectors and trace them */
    particle_state* p;
    bbnbi_simulate(
        &sim, nprt, t1, t2, B_offload_array, plasma_offload_array,
        neutral_offload_array, wall_offload_array, wall_int_offload_array,
        asigma_offload_array, nbi_offload_array, &p, diag_offload_array);

    /* Write output */
    if(sim.mpi_rank == sim.mpi_root) {
        if( hdf5_interface_write_state(sim.hdf5_out, "state", nprt, p) ) {
            print_out0(VERBOSE_MINIMAL, sim.mpi_rank, sim.mpi_root,
                       "\n"
                       "Writing marker state failed.\n"
                       "See stderr for details.\n"
                       "\n");
        }
        free(p);
        print_out0(VERBOSE_NORMAL, sim.mpi_rank, sim.mpi_root,
                   "\nMarker state written.\n");

        hdf5_interface_write_diagnostics(
            &sim, diag_offload_array, sim.hdf5_out);
    }
    print_out0(VERBOSE_MINIMAL, sim.mpi_rank, sim.mpi_root, "\nDone\n");

    return 0;
}