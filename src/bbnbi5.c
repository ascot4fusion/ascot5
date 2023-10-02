/**
 * @file bbnbi5.c
 * @brief BBNBI5 main program
 *
 * BBNBI5 models neutral beam injectors and is used to evaluate shine-through
 * and beam birth-profile. Neutral markers are generated from injector geometry
 * and traced until they ionize or hit the wall. Several injectors can be
 * modelled simultaneously keeping in mind that in this case the output
 * the injector from which a particle originated is lost.
 */
#define _XOPEN_SOURCE 500 /**< drand48 requires POSIX 1995 standard */
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
#include "consts.h"
#include "gitver.h"
#include "math.h"
#include "hdf5_interface.h"
#include "hdf5io/hdf5_helpers.h"
#include "hdf5io/hdf5_marker.h"
#include "print.h"
#include "simulate.h"
#include "random.h"
#include "particle.h"
#include "B_field.h"
#include "plasma.h"
#include "wall.h"
#include "nbi.h"
#include "diag.h"

int read_arguments(int argc, char** argv, sim_offload_data* sim, int* nprt,
                   real* t1, real* t2, int* writemarker);

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
    int nprt; /* Number of markers to be generated excluding shinetrough */
    real t1, t2; /* Markers are initialized in this time-spawn */
    int writemarker; /* Store and write marker input */
    sim_offload_data sim;
    if(read_arguments(argc, argv, &sim, &nprt, &t1, &t2, &writemarker) != 0) {
        abort();
        return 1;
    }

    /* QID for this run */
    char qid[11];
    hdf5_generate_qid(qid);

    int mpi_rank = 0; /* BBNBI 5 does not yet support MPI */
    int mpi_root = 0;
    print_out0(VERBOSE_MINIMAL, mpi_rank, "BBNBI5\n");

#ifdef GIT_VERSION
    print_out0(VERBOSE_MINIMAL, mpi_rank,
               "Tag %s\nBranch %s\n\n", GIT_VERSION, GIT_BRANCH);
#else
    print_out0(VERBOSE_MINIMAL, mpi_rank, "Not under version control\n\n");
#endif

    /* Read data needed for nbi simulation */
    real* nbi_offload_array;
    real* B_offload_array;
    real* plasma_offload_array;
    real* wall_offload_array;
    int*  wall_int_offload_array;
    real* diag_offload_array;
    if( hdf5_interface_read_input(&sim, hdf5_input_bfield | hdf5_input_plasma |
                                  hdf5_input_wall | hdf5_input_nbi |
                                  hdf5_input_options,
                                  &B_offload_array, NULL, &plasma_offload_array,
                                  NULL, &wall_offload_array,
                                  &wall_int_offload_array, NULL, NULL, NULL,
                                  &nbi_offload_array, NULL, NULL) ) {
        print_out0(VERBOSE_MINIMAL, mpi_rank, "Input initialization failed\n");
        abort();
        return 1;
    }
    /* These should be taken from options but keeping them here for now */
    sim.diag_offload_data.diagorb_collect   = 0;
    sim.diag_offload_data.dist5D_collect    = 1;
    sim.diag_offload_data.dist6D_collect    = 0;
    sim.diag_offload_data.distrho5D_collect = 0;
    sim.diag_offload_data.distrho6D_collect = 0;
    sim.diag_offload_data.distCOM_collect   = 0;
    sim.diag_offload_data.diagtrcof_collect = 0;
    diag_init_offload(&sim.diag_offload_data, &diag_offload_array, 0);
    real diag_offload_array_size = sim.diag_offload_data.offload_array_length
        * sizeof(real) / (1024.0*1024.0);
    print_out0(VERBOSE_IO, mpi_rank,
               "Initialized diagnostics, %.1f MB.\n", diag_offload_array_size);
    simulate_init_offload(&sim);

    /* Write bbnbi run group to HDF5 */
    if(mpi_rank == mpi_root) {
        print_out0(VERBOSE_IO, mpi_rank, "\nPreparing output.\n");
        if( hdf5_interface_init_results(&sim, qid, "bbnbi") ) {
            print_out0(VERBOSE_MINIMAL, mpi_rank,
                       "\nInitializing output failed.\n"
                       "See stderr for details.\n");
            /* Free offload data and terminate */
            abort();
            return 1;
        }
        strcpy(sim.qid, qid);
    }

    /* Initialize input data */
    sim_data sim_data;
    sim_init(&sim_data, &sim);
    B_field_init(&sim_data.B_data, &sim.B_offload_data, B_offload_array);
    plasma_init(&sim_data.plasma_data, &sim.plasma_offload_data,
                plasma_offload_array);
    wall_init(&sim_data.wall_data, &sim.wall_offload_data, wall_offload_array,
              wall_int_offload_array);
    random_init(&sim_data.random_data, time(NULL));
    nbi_init(&sim_data.nbi_data, &sim.nbi_offload_data, nbi_offload_array);
    diag_init(&sim_data.diag_data, &sim.diag_offload_data, diag_offload_array);

    /* The number of markers generated per injector is proportional to injector
     * power */
    real total_power = 0;
    for(int i=0; i < sim_data.nbi_data.ninj; i++) {
        total_power += sim_data.nbi_data.inj[i].power;
    }

    /* Simulate requested number of markers into array of particle structs */
    particle* p = NULL;
    if(writemarker) {
        (particle*) malloc(nprt*sizeof(particle));
    }
    int nprt_generated = 0;

    for(int i = 0; i < sim_data.nbi_data.ninj; i++) {
        int nprt_inj;
        if(i == sim_data.nbi_data.ninj-1) {
            nprt_inj = nprt - nprt_generated;
        }
        else {
            nprt_inj = sim_data.nbi_data.inj[i].power/total_power * nprt;
        }
        nbi_generate(&p[nprt_generated], nprt_inj, t1, t2,
                     &(sim_data.nbi_data.inj[i]), &sim_data.B_data,
                     &sim_data.plasma_data, &sim_data.wall_data,
                     &sim_data.random_data, &sim_data.diag_data);

        nprt_generated += nprt_inj;
        print_out0(VERBOSE_NORMAL, mpi_rank,
                   "Generated %d markers for injector %d.\n",
                   nprt_inj, i+1);
    }
    print_out0(VERBOSE_IO, mpi_rank, "\nWriting %d markers.\n",
               nprt_generated);

    if(mpi_rank == mpi_root) {
        hdf5_interface_write_diagnostics(
            &sim, diag_offload_array, sim.hdf5_out);
    }

    /* Copy markers from particle structs into input_particle structs to be
     * written into the h5 file */

    input_particle* ip = (input_particle*) malloc(nprt*sizeof(input_particle));
    for(int i=0; i < nprt; i++) {
        ip[i].type = input_particle_type_p;
        ip[i].p = p[i];
        ip[i].p.id = i+1;
    }

    /* Write marker output */
    if(writemarker) {
        char qid[11];
        hdf5_generate_qid(qid);

        hid_t of = hdf5_create(sim.hdf5_out);
        hdf5_close(of);
        of = hdf5_open(sim.hdf5_out);
        hdf5_marker_write_particle(of, nprt, ip, qid);

        char path[256];
        hdf5_gen_path("/marker/prt_XXXXXXXXXX", qid, path);

        char desc[256];
        sprintf(desc, "BBNBIRUN%s", sim.qid);
        hdf5_write_string_attribute(of, path, "description", desc);

        time_t t = time(NULL);
        struct tm tm = *localtime(&t);
        char date[21];
        sprintf(date, "%04d-%02d-%02d %02d:%02d:%02d.", tm.tm_year + 1900,
                tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
        hdf5_write_string_attribute(of, path, "date",  date);

        /* Set this run as active. */
        hdf5_write_string_attribute(of, "/marker", "active",  qid);
        hdf5_close(of);
    }

    print_out0(VERBOSE_MINIMAL, mpi_rank, "\nDone\n");

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
 * @param writemarker pointer to flag indicating whether a marker input
 *        is created
 *
 * @return Zero if success
 */
int read_arguments(int argc, char** argv, sim_offload_data* sim, int* nprt,
                   real* t1, real* t2, int* writemarker) {
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
        {"markers",  required_argument, 0, 13},
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
    *writemarker        = 1;
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
            case 13:
                *writemarker = atof(optarg);
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
                print_out(VERBOSE_MINIMAL,
                          "--marker flag indicating whether a marker input is written (default: 1)\n");
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
