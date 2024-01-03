/**
 * @file bmc_main.c
 * @brief Backward Monte Carlo simulation program
 */
#define _XOPEN_SOURCE 500 /**< drand48 requires POSIX 1995 standard */
#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ascot5.h"
#include "consts.h"
#include "math.h"
#include "wall.h"
#include "diag.h"
#include "B_field.h"
#include "plasma.h"
#include "print.h"
#include "simulate.h"
#include "hdf5_interface.h"
#include "offload.h"
#include "mpi_interface.h"
#include "bmc_mesh.h"
#include "simulate/simulate_bmc.h"

int read_arguments(int argc, char** argv, sim_offload_data* sim);

/**
 * @brief Main function for bmc_main
 *
 * @param  argc argument count of the command line arguments
 * @param  argv argument vector of the command line arguments
 *
 * @return Zero if simulation was completed
 */
int main(int argc, char** argv) {

    /* Read and parse command line arguments */
    sim_offload_data sim;
    if( read_arguments(argc, argv, &sim) ) {
        abort();
        return 1;
    }

    char qid[11];
    hdf5_generate_qid(qid);

    /* Get MPI rank and set qid for the run*/
    int mpi_rank, mpi_size, mpi_root;
    mpi_interface_init(argc, argv, &sim, &mpi_rank, &mpi_size, &mpi_root);

    print_out0(VERBOSE_MINIMAL, mpi_rank, "BMC_MAIN\n");

    print_out0(VERBOSE_NORMAL, mpi_rank,
               "Initialized MPI, rank %d, size %d.\n", mpi_rank, mpi_size);

    /* Offload data arrays that are allocated when input is read */
    real* B_offload_array;
    real* E_offload_array;
    real* plasma_offload_array;
    real* wall_offload_array;
    int* wall_int_offload_array;
    real* diag_offload_array;

    /* Read input from the HDF5 file */
    if( hdf5_interface_read_input(&sim,
                                  hdf5_input_options | hdf5_input_bfield |
                                  hdf5_input_efield  | hdf5_input_plasma |
                                  hdf5_input_wall,
                                  &B_offload_array, &E_offload_array,
                                  &plasma_offload_array, NULL,
                                  &wall_offload_array,  &wall_int_offload_array,
                                  NULL, NULL, NULL, NULL,
                                  NULL, NULL) ) {
        print_out0(VERBOSE_MINIMAL, mpi_rank,
                   "\nInput reading or initializing failed.\n"
                   "See stderr for details.\n");
        abort();
        return 1;
    };

    /* No offloading implemented so inputs are initialized now */
    sim_data sim_data;
    sim_init(&sim_data, &sim);
    random_init(&sim_data.random_data, time(NULL));
    B_field_init(&sim_data.B_data, &sim.B_offload_data, B_offload_array);
    E_field_init(&sim_data.E_data, &sim.E_offload_data, E_offload_array);
    plasma_init(&sim_data.plasma_data, &sim.plasma_offload_data,
                plasma_offload_array);
    wall_init(&sim_data.wall_data, &sim.wall_offload_data, wall_offload_array,
              wall_int_offload_array);
    diag_init_offload(&sim.diag_offload_data, &diag_offload_array, 0);

    print_out0(VERBOSE_NORMAL, mpi_rank, "\nInitializing probability mesh.\n");

    bmc_mesh mesh;
    if(bmc_mesh_init(
            sim.diag_offload_data.dist5D.min_r,
            sim.diag_offload_data.dist5D.max_r,
            sim.diag_offload_data.dist5D.n_r,
            sim.diag_offload_data.dist5D.min_phi,
            sim.diag_offload_data.dist5D.max_phi,
            sim.diag_offload_data.dist5D.n_phi,
            sim.diag_offload_data.dist5D.min_z,
            sim.diag_offload_data.dist5D.max_z,
            sim.diag_offload_data.dist5D.n_z,
            sim.diag_offload_data.dist5D.min_ppara,
            sim.diag_offload_data.dist5D.max_ppara,
            sim.diag_offload_data.dist5D.n_ppara,
            sim.diag_offload_data.dist5D.min_pperp,
            sim.diag_offload_data.dist5D.max_pperp,
            sim.diag_offload_data.dist5D.n_pperp,
            &mesh)) {
        goto CLEANUP_FAILURE;
    }

    /* Initialize results group in the output file */
    if(mpi_rank == mpi_root) {
        print_out0(VERBOSE_IO, mpi_rank, "\nPreparing output.\n")
	  if( hdf5_interface_init_results(&sim, qid, "bmc") ) {
            print_out0(VERBOSE_MINIMAL, mpi_rank,
                    "\nInitializing output failed.\n"
                    "See stderr for details.\n");
            /* Free offload data and terminate */
            goto CLEANUP_FAILURE;
        };
        strcpy(sim.qid, qid);
    }

    /* Run simulation in either time-dependent or -independent mode */
    if(!sim_data.bmc_timedependent) {
        /* Evaluate the push matrix once */
        size_t start = 0, stop = mesh.size;
        real* r     = (real*) malloc( mesh.size * HERMITE_KNOTS * sizeof(real));
        real* phi   = (real*) malloc( mesh.size * HERMITE_KNOTS * sizeof(real));
        real* z     = (real*) malloc( mesh.size * HERMITE_KNOTS * sizeof(real));
        real* ppara = (real*) malloc( mesh.size * HERMITE_KNOTS * sizeof(real));
        real* pperp = (real*) malloc( mesh.size * HERMITE_KNOTS * sizeof(real));
        int* fate = (int*) malloc( mesh.size * HERMITE_KNOTS * sizeof(int));

        simulate_bmc_gc(&sim_data, &mesh, sim_data.bmc_timestep,
                        sim_data.bmc_tstart, start, stop,
                        r, phi, z, ppara, pperp, fate);
        for(real t=sim_data.bmc_tstart; t <= sim_data.bmc_tstop;
            t += sim_data.bmc_timestep) {
            /* Update the probability repeatedly until the simulation is
             * complete */
            bmc_mesh_update(&mesh, start, stop, r, phi, z, ppara, pperp, fate);
            bmc_mesh_finishstep(&mesh);
            real tot = 0;
            for(int i=0; i<mesh.size; i++) {
                //tot = fmax(mesh.val_prev[i], tot);
                tot += mesh.val_prev[i];
            }
            //printf("%g\n", tot);
        }
        free(r);
        free(phi);
        free(z);
        free(ppara);
        free(pperp);
        free(fate);
    }
    else {
        for(size_t i = 0; i < 10; i++) {
            /* Evaluate the push matrix */

            /* Update the probability */
        }
    }
    for(size_t i=0; i<mesh.size; i++) {
        diag_offload_array[i] = mesh.val_prev[i];
    }

    mpi_interface_finalize();

    /* Write output */
    hdf5_interface_write_diagnostics(
            &sim, diag_offload_array, sim.hdf5_out);

    print_out0(VERBOSE_MINIMAL, mpi_rank, "\nDone.\n");

    return 0;


/* Free resources in case simulation crashes */
CLEANUP_FAILURE:

    mpi_interface_finalize();

    abort();
    return 1;
}


/**
 * @brief Read command line arguments and modify sim struct accordingly
 *
 * The command line arguments are in, out, mpi_size, and mpi_rank which
 * correspond to input file, output file, number of MPI processes and
 * the rank of this process. These are stored in simulation offload data struct
 * as (default values, used if the specific argument was not given, are in
 * parenthesis):
 *
 * - sim->hdf5_in     = "in.h5" ("ascot.h5")
 * - sim->hdf5_out    = "out" (sim->hdf5_in is copied here)
 * - sim->mpi_rank    = 0
 * - sim->mpi_size    = 0
 * - sim->desc        = "No description"
 *
 * If the arguments could not be parsed, this function returns a non-zero exit
 * value.
 *
 * @param argc argument count as given to main()
 * @param argv argument vector as given to main()
 * @param sim pointer to offload data struct
 *
 * @return Zero if success
 */
int read_arguments(int argc, char** argv, sim_offload_data* sim) {
    struct option longopts[] = {
        {"in", required_argument, 0, 1},
        {"out", required_argument, 0, 2},
        {"mpi_size", required_argument, 0, 3},
        {"mpi_rank", required_argument, 0, 4},
        {"d", required_argument, 0, 5},
        {"options", required_argument, 0, 6},
        {"bfield",  required_argument, 0, 7},
        {"efield",  required_argument, 0, 8},
        {"marker",  required_argument, 0, 9},
        {"wall",    required_argument, 0, 10},
        {"plasma",  required_argument, 0, 11},
        {"neutral", required_argument, 0, 12},
        {"boozer",  required_argument, 0, 13},
        {"mhd",     required_argument, 0, 14},
        {0, 0, 0, 0}
    };

    // Initialize default values
    sim->hdf5_in[0]     = '\0';
    sim->hdf5_out[0]    = '\0';
    sim->mpi_rank       = 0;
    sim->mpi_size       = 0;
    strcpy(sim->description, "No description.");
    sim->qid_options[0] = '\0';
    sim->qid_bfield[0]  = '\0';
    sim->qid_efield[0]  = '\0';
    sim->qid_marker[0]  = '\0';
    sim->qid_wall[0]    = '\0';
    sim->qid_plasma[0]  = '\0';
    sim->qid_neutral[0] = '\0';
    sim->qid_boozer[0]  = '\0';
    sim->qid_mhd[0]     = '\0';

    // Read user input
    int c;
    int slen;  // String length
    while((c = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
        switch(c) {
            case 1:
                /* The .hdf5 filename can be specified with or without the
                 * trailing .h5 */
                slen = strlen(optarg);
                if ( slen > 3 && !strcmp(optarg+slen-3, ".h5") ) {
                    strcpy(sim->hdf5_in, optarg);
                    (sim->hdf5_in)[slen-3]='\0';
                }
                else {
                    strcpy(sim->hdf5_in, optarg);
                }
                break;
            case 2:
                /* The .hdf5 filename can be specified with or without the
                 * trailing .h5 */
                slen = strlen(optarg);
                if ( slen > 3 && !strcmp(optarg+slen-3, ".h5") ) {
                    strcpy(sim->hdf5_out, optarg);
                    (sim->hdf5_out)[slen-3]='\0';
                }
                else {
                    strcpy(sim->hdf5_out, optarg);
                }
                break;
            case 3:
                sim->mpi_size = atoi(optarg);
                break;
            case 4:
                sim->mpi_rank = atoi(optarg);
                break;
            case 5:
                strcpy(sim->description, optarg);
                break;
            case 6:
                strcpy(sim->qid_options, optarg);
                break;
            case 7:
                strcpy(sim->qid_bfield, optarg);
                break;
            case 8:
                strcpy(sim->qid_efield, optarg);
                break;
            case 9:
                strcpy(sim->qid_marker, optarg);
                break;
            case 10:
                strcpy(sim->qid_wall, optarg);
                break;
            case 11:
                strcpy(sim->qid_plasma, optarg);
                break;
            case 12:
                strcpy(sim->qid_neutral, optarg);
                break;
            case 13:
                strcpy(sim->qid_boozer, optarg);
                break;
            case 14:
                strcpy(sim->qid_mhd, optarg);
                break;
            default:
                // Unregonizable argument(s). Tell user how to run ascot5_main
                print_out(VERBOSE_MINIMAL,
                          "\nUnrecognized argument. The valid arguments are:\n");
                print_out(VERBOSE_MINIMAL,
                          "--in input file without .h5 extension (default: ascot)\n");
                print_out(VERBOSE_MINIMAL,
                          "--out output file without .h5 extension (default: same as input)\n");
                print_out(VERBOSE_MINIMAL,
                          "--mpi_size number of independent processes\n");
                print_out(VERBOSE_MINIMAL,
                          "--mpi_rank rank of independent process\n");
                print_out(VERBOSE_MINIMAL,
                          "--d run description maximum of 250 characters\n");
                return 1;
        }
    }

    /* Default value for input file is ascot.h5, and for output same as input
     * file. Adujust hdf5_in and hdf5_out accordingly. For output file, we don't
     * add the .h5 extension here. */
    if(sim->hdf5_in[0] == '\0' && sim->hdf5_out[0] == '\0') {
        // No input, use default values for both
        strcpy(sim->hdf5_in, "ascot.h5");
        strcpy(sim->hdf5_out, "ascot.h5");
    }
    else if(sim->hdf5_in[0] == '\0' && sim->hdf5_out[0] != '\0') {
        // Output file is given but the input file is not
        strcpy(sim->hdf5_in, "ascot.h5");
    }
    else if(sim->hdf5_in[0] != '\0' && sim->hdf5_out[0] == '\0') {
        // Input file is given but the output is not
        strcat(sim->hdf5_in, ".h5");
        strcpy(sim->hdf5_out, sim->hdf5_in);
    }
    else {
        // Both input and output files are given
        strcat(sim->hdf5_in, ".h5");
        strcat(sim->hdf5_out, ".h5");
    }
    return 0;
}