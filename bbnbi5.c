/**
 * @file bbnbi5.c
 * @brief BBNBI5 main program
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
#include "B_field.h"
#include "consts.h"
#include "hdf5_interface.h"
#include "hdf5io/hdf5_helpers.h"
#include "hdf5io/hdf5_nbi.h"
#include "hdf5io/hdf5_marker.h"
#include "math.h"
#include "nbi.h"
#include "particle.h"
#include "plasma.h"
#include "print.h"
#include "random.h"
#include "wall.h"

int read_arguments(int argc, char** argv, sim_offload_data* sim, int* nprt);

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
    int nprt;
    sim_offload_data sim;
    if(read_arguments(argc, argv, &sim, &nprt) != 0) {
        return 1;
    }

    /* Initialize data needed for nbi simulation */
    real* B_offload_array;
    real* plasma_offload_array;
    real* wall_offload_array;
    hdf5_interface_read_input(&sim, hdf5_input_bfield | hdf5_input_plasma |
                              hdf5_input_wall, &B_offload_array, NULL,
                              &plasma_offload_array, NULL, &wall_offload_array,
                              NULL, NULL, NULL, NULL);

    B_field_data B_data;
    B_field_init(&B_data, &sim.B_offload_data, B_offload_array);

    plasma_data plasma_data;
    plasma_init(&plasma_data, &sim.plasma_offload_data, plasma_offload_array);

    wall_data wall_data;
    wall_init(&wall_data, &sim.wall_offload_data, wall_offload_array);

    random_data rng;
    random_init(&rng, time(NULL));

    /* NBI data read and initialized separately for now */
    int n_inj;
    nbi_injector* inj;
    hid_t f = hdf5_open(sim.hdf5_in);
    hdf5_nbi_read(f, &n_inj, &inj);
    hdf5_close(f);

    real total_power = 0;

    for(int i=0; i < n_inj; i++) {
        printf("Injector %d:\n", i+1);
        printf("id: %d\n", inj[i].id);
        printf("n_beamlet: %d\n", inj[i].n_beamlet);
        printf("power: %le\n", inj[i].power);
        printf("energy: %le\n", inj[i].energy);
        printf("efrac: %le %le %le\n", inj[i].efrac[0], inj[i].efrac[1], inj[i].efrac[2]);
        printf("divergence: %le %le %le %le %le\n", inj[i].div_h, inj[i].div_v, inj[i].div_halo_frac, inj[i].div_halo_h, inj[i].div_halo_v);
        printf("anum: %d\n", inj[i].anum);
        printf("znum: %d\n", inj[i].znum);
        printf("mass: %le\n", inj[i].mass);
        printf("\n");

        total_power += inj[i].power;
    }


    /* Simulate requested number of markers into array of particle structs */
    particle* p = (particle*) malloc(nprt*sizeof(particle));
    int nprt_generated = 0;

    for(int i = 0; i < n_inj; i++) {
        int nprt_inj;

        if(i == n_inj-1) {
            nprt_inj = nprt - nprt_generated;
        }
        else {
            nprt_inj = inj[i].power/total_power * nprt;
        }

        nbi_generate(nprt_inj, &p[nprt_generated], &inj[i], &B_data,
                     &plasma_data, &wall_data, &rng);

        nprt_generated += nprt_inj;
        printf("Generated %d markers for injector %d.\n", nprt_inj, i+1);
    }

    printf("\nWriting %d markers.\n", nprt_generated);

    /* Copy markers from particle structs into input_particle structs to be
     * written into the h5 file */
    input_particle* ip = (input_particle*) malloc(nprt*sizeof(input_particle));
    for(int i=0; i < nprt; i++) {
        ip[i].type = input_particle_type_p;
        ip[i].p = p[i];
    }

    char qid[11];
    hdf5_generate_qid(qid);

    hid_t of = hdf5_create(sim.hdf5_out);
    hdf5_close(of);
    of = hdf5_open(sim.hdf5_out);
    hdf5_marker_write_particle(of, nprt, ip, qid);

    /* Write metadata */
    char path[256];
    hdf5_gen_path("/marker/prt_XXXXXXXXXX", qid, path);

    hdf5_write_string_attribute(f, path, "description",  sim.description);

    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    char date[21];
    sprintf(date, "%04d-%02d-%02d %02d:%02d:%02d.", tm.tm_year + 1900,
            tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    hdf5_write_string_attribute(f, path, "date",  date);

    /* Set this run as active. */
    hdf5_write_string_attribute(f, "/marker", "active",  qid);

    hdf5_close(of);

    printf("\nDone.\n");

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
 *
 * @return Zero if success
 */
int read_arguments(int argc, char** argv, sim_offload_data* sim, int* nprt) {
    struct option longopts[] = {
        {"in", required_argument, 0, 1},
        {"out", required_argument, 0, 2},
        {"mpi_size", required_argument, 0, 3},
        {"mpi_rank", required_argument, 0, 4},
        {"d", required_argument, 0, 5},
        {"bfield",  required_argument, 0, 6},
        {"wall",    required_argument, 0, 7},
        {"plasma",  required_argument, 0, 8},
        {"n",       required_argument, 0, 9},
        {0, 0, 0, 0}
    };

    // Initialize default values
    sim->hdf5_in[0]     = '\0';
    sim->hdf5_out[0]    = '\0';
    sim->mpi_rank       = 0;
    sim->mpi_size       = 0;
    strcpy(sim->description, "No description.");
    sim->qid_bfield[0]  = '\0';
    sim->qid_wall[0]    = '\0';
    sim->qid_plasma[0]  = '\0';
    *nprt               = 10000;

    // Read user input
    int c;
    int slen;  // String length
    while((c = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
        switch(c) {
            case 1:
                // The .hdf5 filename can be specified with or without the trailing .h5
                slen = strlen(optarg);
                if ( slen > 3 && !strcmp(optarg+slen-3,".h5") ) {
                    strncpy(sim->hdf5_in,optarg,slen-3);
                    (sim->hdf5_in)[slen-3]='\0';
                }
                else
                    strcpy(sim->hdf5_in, optarg);
                break;
            case 2:
                // The .hdf5 filename can be specified with or without the trailing .h5
                slen = strlen(optarg);
                if ( slen > 3 && !strcmp(optarg+slen-3,".h5") ) {
                    strncpy(sim->hdf5_out,optarg,slen-3);
                    (sim->hdf5_out)[slen-3]='\0';
                }
                else
                    strcpy(sim->hdf5_out, optarg);
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
                *nprt = atoi(optarg);
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
                print_out(VERBOSE_MINIMAL,
                    "--n number of markers to generate, (default: 10000)\n");
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
