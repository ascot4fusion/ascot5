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

int bbnbi_read_arguments(int argc, char** argv, sim_offload_data* sim,
                         int* nprt, real* t1, real* t2);

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

/**
 * @brief Read command line arguments
 *
 * Read in command line arguments, input and output names and mpi parameters
 * are stored in sim structure as with ascot5, number of markers is passed
 * as an argument.
 *
 * If the arguments could not be parsed, this function returns a non-zero exit
 * value.
 *
 * @param argc argument count as given to main()
 * @param argv argument vector as given to main()
 * @param sim pointer to offload data struct
 * @param nprt pointer to integer where number of markers is stored
 * @param t1 pointer to store beginning of time interval
 * @param t2 pointer to store end of the time interval
 *
 * @return Zero if success
 */
int bbnbi_read_arguments(int argc, char** argv, sim_offload_data* sim,
                         int* nprt, real* t1, real* t2) {
    struct option longopts[] = {
        {"in",       required_argument, 0,  1},
        {"out",      required_argument, 0,  2},
        {"mpi_size", required_argument, 0,  3},
        {"mpi_rank", required_argument, 0,  4},
        {"d",        required_argument, 0,  5},
        {"bfield",   required_argument, 0,  6},
        {"wall",     required_argument, 0,  7},
        {"plasma",   required_argument, 0,  8},
        {"nbi",      required_argument, 0,  9},
        {"n",        required_argument, 0, 10},
        {"t1",       required_argument, 0, 11},
        {"t2",       required_argument, 0, 12},
        {0, 0, 0, 0}
    };

    /* Initialize default values */
    sim->hdf5_in[0]     = '\0';
    sim->hdf5_out[0]    = '\0';
    sim->qid_options[0] = '\0';
    sim->qid_bfield[0]  = '\0';
    sim->qid_efield[0]  = '\0';
    sim->qid_marker[0]  = '\0';
    sim->qid_wall[0]    = '\0';
    sim->qid_plasma[0]  = '\0';
    sim->qid_neutral[0] = '\0';
    sim->qid_boozer[0]  = '\0';
    sim->qid_mhd[0]     = '\0';
    sim->qid_asigma[0]  = '\0';
    sim->qid_nbi[0]     = '\0';
    sim->mpi_rank       = 0;
    sim->mpi_size       = 0;
    *nprt               = 10000;
    *t1                 = 0.0;
    *t2                 = 0.0;
    strcpy(sim->description, "TAG");

    /* Read user input */
    int c;
    int slen;
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
                strcpy(sim->qid_bfield, optarg);
                break;
            case 7:
                strcpy(sim->qid_wall, optarg);
                break;
            case 8:
                strcpy(sim->qid_plasma, optarg);
                break;
            case 9:
                strcpy(sim->qid_nbi, optarg);
                break;
            case 10:
                *nprt = atoi(optarg);
                break;
            case 11:
                *t1 = atof(optarg);
                break;
            case 12:
                *t2 = atof(optarg);
                break;
            default:
                // Unregonizable argument(s). Tell user how to run ascot5_main
                print_out(VERBOSE_MINIMAL,
                          "\nUnrecognized argument. The valid arguments are:\n");
                print_out(VERBOSE_MINIMAL,
                          "--in input file (default: ascot.h5)\n");
                print_out(VERBOSE_MINIMAL,
                          "--out output file (default: same as input)\n");
                print_out(VERBOSE_MINIMAL,
                          "--mpi_size number of independent processes\n");
                print_out(VERBOSE_MINIMAL,
                          "--mpi_rank rank of independent process\n");
                print_out(VERBOSE_MINIMAL,
                          "--d run description maximum of 250 characters\n");
                print_out(VERBOSE_MINIMAL,
                          "--n number of markers to generate (default: 10000)\n");
                print_out(VERBOSE_MINIMAL,
                          "--t1 time when injectors are turned on (default: 0.0 s)\n");
                print_out(VERBOSE_MINIMAL,
                          "--t2 time when injectors are turned off (default: 0.0 s)\n");
                return 1;
        }
    }

    /* Default value for input file is ascot.h5, and for output same as input
     * file. Adujust hdf5_in and hdf5_out accordingly. For output file, we don't
     * add the .h5 extension here. */
    if(sim->hdf5_in[0] == '\0' && sim->hdf5_out[0] == '\0') {
        /* No input, use default values for both */
        strcpy(sim->hdf5_in,  "ascot.h5");
        strcpy(sim->hdf5_out, "ascot.h5");
    }
    else if(sim->hdf5_in[0] == '\0' && sim->hdf5_out[0] != '\0') {
        /* Output file is given but the input file is not */
        strcpy(sim->hdf5_in,  "ascot.h5");
        strcat(sim->hdf5_out, ".h5");
    }
    else if(sim->hdf5_in[0] != '\0' && sim->hdf5_out[0] == '\0') {
        /* Input file is given but the output is not */
        strcat(sim->hdf5_in, ".h5");
        strcpy(sim->hdf5_out, sim->hdf5_in);
    }
    else {
        /* Both input and output files are given */
        strcat(sim->hdf5_in,  ".h5");
        strcat(sim->hdf5_out, ".h5");
    }
    return 0;
}
