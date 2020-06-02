/**
 * @file ascot5_main.c
 * @brief ASCOT5 stand-alone program
 *
 * This program reads data from input HDF5 file, simulates the given markers,
 * and writes the output data to a HDF5 file.
 *
 * The input and output files can be separate.
 *
 * Example:
 *
 *     ascot5_main --in=in --out=out
 *
 * Here "in" refers to in.h5 where input data is located and "out" to out.h5
 * where results will be stored. If no input argument is given the data is read
 * from ascot.h5. If not output argument is given the results are stored in the
 * input file.
 *
 * This program assumes that the input file contains magnetic field, electric
 * field, plasma, wall, and neutral data along with markers and options. See
 * hdf5_input.c for details. This program uses the input fields that are marked
 * as active (the HDF5 file can contain multiple instances of same input types
 * but only the active one is used here).
 *
 * The results are stored under /results/ group in output HDF5 file. The group
 * is created if one does not exists. For each run a specific "run" group is
 * created, which has the format run-XXXXXXXXXX, where "XXXXXXXXXX" is randomly
 * generated identification number (QID). The run group holds information when
 * the run was started, which input fields were used (referenced by their QIDs),
 * and at least the marker initial and end states if the simulation succeeded.
 * Also any other diagnostic data that was used is stored there.
 *
 * This program uses MPI by dividing the number of markers equally to all MPI
 * processes. The markers are not suffled so user is advised to do it beforehand
 * to ensure work is evenly distributed. A single MPI process can be simulated
 * with:
 *
 *     ascot5_main --mpi_size=size --mpi_rank=rank
 *
 * where size refers to number of MPI processes and rank is the process being
 * run (between [0, size-1]). Running the program this way does not use MPI.
 * This is intended to be used in Condor-like environments.
 *
 * You can add a description of the simulation as:
 *
 * ascot5_main --d="This is a test run"
 *
 * which is written in HDF5 file at the run group specific to this simulation.
 *
 * In addition to output data, the simulation progress may be written in
 * *.stdout files with each MPI process having dedicated file. See ascot5.h for
 * details.
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
#include "particle.h"
#include "endcond.h"
#include "hdf5_interface.h"
#include "offload.h"
#include "gitver.h"
#include "mpi_interface.h"

#ifdef TRAP_FPE
#include <fenv.h>
#endif

int read_arguments(int argc, char** argv, sim_offload_data* sim);
void marker_summary(particle_state* p, int n);

/**
 * @brief Main function for ascot5_main
 *
 * This function calls functions that read input data from the disk, and
 * functions that initialize the offload data structs and offload arrays.
 * Actual simulation is done by calling simulate(). Once the simulation has been
 * completed, offload arrays are deallocated and the results are written to the
 * disk.
 *
 * MPI level parallelisation is done here as well as the offloading.
 *
 * @param  argc argument count of the command line arguments
 * @param  argv argument vector of the command line arguments
 *
 * @return Zero if simulation was completed
 */
int main(int argc, char** argv) {

#ifdef TRAP_FPE
	/* This will raise floating point exceptions */
    feenableexcept(FE_DIVBYZERO| FE_INVALID | FE_OVERFLOW);
	/*
	 * If you are hunting a specific exception, you can disable the exceptions in other parts
	 * of the code by surrounding it with: */
    /*
     * fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
     *  --- your  code here ---
     * feenableexcept(FE_DIVBYZERO  | FE_INVALID | FE_OVERFLOW);
     *
     * */

#endif


    /* Read and parse command line arguments */
    sim_offload_data sim;
    if( read_arguments(argc, argv, &sim) ) {
        abort();
        return 1;
    }

    /* Get MPI rank and set qid for the run*/
    char qid[11];
    hdf5_generate_qid(qid);

    int mpi_rank, mpi_size, mpi_root;
    mpi_interface_init(argc, argv, &sim, &mpi_rank, &mpi_size, &mpi_root);

    print_out0(VERBOSE_MINIMAL, mpi_rank,
               "ASCOT5_MAIN\n");

#ifdef GIT_VERSION
    print_out0(VERBOSE_MINIMAL, mpi_rank,
               "Tag %s\nBranch %s\n\n", GIT_VERSION, GIT_BRANCH);
#else
    print_out0(VERBOSE_MINIMAL, mpi_rank,
               "Not under version control\n\n");
#endif

    print_out0(VERBOSE_NORMAL, mpi_rank,
               "Initialized MPI, rank %d, size %d.\n", mpi_rank, mpi_size);

    /* Total number of markers to be simulated */
    int n_tot;

    /* Marker input struct */
    input_particle* p;

    /* Offload data arrays that are allocated when input is read */
    real* B_offload_array;
    real* E_offload_array;
    real* plasma_offload_array;
    real* neutral_offload_array;
    real* wall_offload_array;
    real* boozer_offload_array;
    real* mhd_offload_array;

    /* Read input from the HDF5 file */
    if( hdf5_interface_read_input(&sim,
                                  hdf5_input_options | hdf5_input_bfield |
                                  hdf5_input_efield  | hdf5_input_plasma |
                                  hdf5_input_neutral | hdf5_input_wall |
                                  hdf5_input_marker | hdf5_input_boozer |
                                  hdf5_input_mhd,
                                  &B_offload_array, &E_offload_array,
                                  &plasma_offload_array, &neutral_offload_array,
                                  &wall_offload_array, &boozer_offload_array,
                                  &mhd_offload_array, &p, &n_tot) ) {
        print_out0(VERBOSE_MINIMAL, mpi_rank,
                   "\nInput reading or initializing failed.\n"
                   "See stderr for details.\n");
        abort();
        return 1;
    };
    simulate_init_offload(&sim);

    /* Pack offload data into single array and free individual offload arrays */
    /* B_offload_array is needed for marker evaluation and is freed later */
    real* offload_array;
    offload_package offload_data;
    offload_init_offload(&offload_data, &offload_array);
    offload_pack(&offload_data, &offload_array, B_offload_array,
                 sim.B_offload_data.offload_array_length);

    offload_pack(&offload_data, &offload_array, E_offload_array,
                 sim.E_offload_data.offload_array_length);
    E_field_free_offload(&sim.E_offload_data, &E_offload_array);

    offload_pack(&offload_data, &offload_array, plasma_offload_array,
                 sim.plasma_offload_data.offload_array_length);
    plasma_free_offload(&sim.plasma_offload_data, &plasma_offload_array);

    offload_pack(&offload_data, &offload_array, neutral_offload_array,
                 sim.neutral_offload_data.offload_array_length);
    neutral_free_offload(&sim.neutral_offload_data, &neutral_offload_array);

    offload_pack(&offload_data, &offload_array, wall_offload_array,
                 sim.wall_offload_data.offload_array_length);
    wall_free_offload(&sim.wall_offload_data, &wall_offload_array);

    offload_pack(&offload_data, &offload_array, boozer_offload_array,
                 sim.boozer_offload_data.offload_array_length);
    boozer_free_offload(&sim.boozer_offload_data, &boozer_offload_array);

    offload_pack(&offload_data, &offload_array, mhd_offload_array,
                 sim.mhd_offload_data.offload_array_length);
    mhd_free_offload(&sim.mhd_offload_data, &mhd_offload_array);

    /* Initialize diagnostics offload data.
     * Separate arrays for host and target */
#ifdef TARGET
    real* diag_offload_array_mic0;
    real* diag_offload_array_mic1;
    diag_init_offload(&sim.diag_offload_data, &diag_offload_array_mic0, n_tot);
    diag_init_offload(&sim.diag_offload_data, &diag_offload_array_mic1, n_tot);
#else
    real* diag_offload_array_host;
    diag_init_offload(&sim.diag_offload_data, &diag_offload_array_host, n_tot);
#endif

    real diag_offload_array_size = sim.diag_offload_data.offload_array_length
        * sizeof(real) / (1024.0*1024.0);
    print_out0(VERBOSE_NORMAL, mpi_rank,
               "Initialized diagnostics, %.1f MB.\n", diag_offload_array_size);

    /* Choose which markers are used in this MPI process. Simply put, markers
     * are divided into mpi_size sequential blocks and the mpi_rank:th block
     * is chosen for this simulation. */
    int start_index, n;
    mpi_my_particles(&start_index, &n, n_tot, mpi_rank, mpi_size);
    p += start_index;

    /* Set up particlestates on host, needs magnetic field evaluation */
    print_out0(VERBOSE_NORMAL, mpi_rank,
               "\nInitializing marker states.\n");
    B_field_data Bdata;
    B_field_init(&Bdata, &sim.B_offload_data, B_offload_array);
    particle_state* ps = (particle_state*) malloc(n * sizeof(particle_state));
    for(int i = 0; i < n; i++) {
        particle_input_to_state(&p[i], &ps[i], &Bdata);
    }
    /* We can now free the Bfield offload array */
    B_field_free_offload(&sim.B_offload_data, &B_offload_array);
    free(p-start_index); // Input markers are no longer required
    print_out0(VERBOSE_NORMAL, mpi_rank,
               "Estimated memory usage %.1f MB.\n",
               (sizeof(real) * n) / (1024.0*1024.0));
    print_out0(VERBOSE_NORMAL, mpi_rank,
               "Marker states initialized.\n");

    if(mpi_rank == mpi_root) {
        /* Initialize results group in the output file */
        print_out0(VERBOSE_IO, mpi_rank, "\nPreparing output.\n")
        if( hdf5_interface_init_results(&sim, qid) ) {
            print_out0(VERBOSE_MINIMAL, mpi_rank,
                       "\nInitializing output failed.\n"
                       "See stderr for details.\n");
            /* Free offload data and terminate */
            goto CLEANUP_FAILURE;
        };
        strcpy(sim.qid, qid);
    }

    particle_state* ps_gathered;
    int n_gathered;

    mpi_gather_particlestate(ps, &ps_gathered, &n_gathered, n_tot,
                             mpi_rank, mpi_size, mpi_root);

    if(mpi_rank == mpi_root) {
        /* Write inistate */
        if(hdf5_interface_write_state(sim.hdf5_out, "inistate", n_gathered,
                                      ps_gathered)) {
            print_out0(VERBOSE_MINIMAL, mpi_rank,
                       "\n"
                       "Writing inistate failed.\n"
                       "See stderr for details.\n"
                       "\n");
            /* Free offload data and terminate */
            goto CLEANUP_FAILURE;
        }
        print_out0(VERBOSE_NORMAL, mpi_rank,
                   "\nInistate written.\n");
    }

    free(ps_gathered);

    /* Divide markers among host and target */
#ifdef TARGET
    int n_mic = n / TARGET;
    int n_host = 0;
#else
    int n_mic = 0;
    int n_host = n;
#endif

    double mic0_start = 0, mic0_end=0,
        mic1_start=0, mic1_end=0,
        host_start=0, host_end=0;

    fflush(stdout);

    /* Allow threads to spawn threads */
    omp_set_nested(1);

    /* Actual marker simulation happens here. Threads are spawned which
     * distribute the execution between target(s) and host. Both input and
     * diagnostic offload arrays are mapped to target. Simulation is initialized
     * at the target and completed within the simulate() function.*/
    #pragma omp parallel sections num_threads(3)
    {
        /* Run simulation on first target */
#if TARGET >= 1
        #pragma omp section
        {
            mic0_start = omp_get_wtime();

            #pragma omp target device(0) map( \
                ps[0:n_mic], \
                offload_array[0:offload_data.offload_array_length], \
                diag_offload_array_mic0[0:sim.diag_offload_data.offload_array_length] \
            )
            simulate(1, n_mic, ps, &sim, &offload_data, offload_array,
                diag_offload_array_mic0);

            mic0_end = omp_get_wtime();
        }
#endif

        /* Run simulation on second target */
#if TARGET >= 2
        #pragma omp section
        {
            mic1_start = omp_get_wtime();

            #pragma omp target device(1) map( \
                ps[n_mic:2*n_mic], \
                offload_array[0:offload_data.offload_array_length], \
                diag_offload_array_mic1[0:sim.diag_offload_data.offload_array_length] \
            )
            simulate(2, n_mic, ps+n_mic, &sim, &offload_data, offload_array,
                diag_offload_array_mic1);

            mic1_end = omp_get_wtime();
        }
#endif

        /* No target, marker simulation happens where the code execution began.
         * Offloading is only emulated. */
#ifndef TARGET
        #pragma omp section
        {
            host_start = omp_get_wtime();
            simulate(0, n_host, ps+2*n_mic, &sim, &offload_data,
                offload_array, diag_offload_array_host);
            host_end = omp_get_wtime();
        }
#endif
    }

    /* Code execution returns to host. */
    print_out0(VERBOSE_NORMAL, mpi_rank, "mic0 %lf s, mic1 %lf s, host %lf s\n",
        mic0_end-mic0_start, mic1_end-mic1_start, host_end-host_start);

    /* Free input data */
    offload_free_offload(&offload_data, &offload_array);

    mpi_gather_particlestate(ps, &ps_gathered, &n_gathered, n_tot,
                             mpi_rank, mpi_size, mpi_root);

    if(mpi_rank == mpi_root) {

        /* Write endstate */
        if( hdf5_interface_write_state(sim.hdf5_out, "endstate", n_gathered,
                                       ps_gathered)) {
            print_out0(VERBOSE_MINIMAL, mpi_rank,
                   "\nWriting endstate failed.\n"
                   "See stderr for details.\n");
            /* Free offload data and terminate */
            goto CLEANUP_FAILURE;
        }
        print_out0(VERBOSE_NORMAL, mpi_rank,
                   "Endstate written.\n");
    }

    /* Combine diagnostic data and write it to HDF5 file */
    print_out0(VERBOSE_MINIMAL, mpi_rank,
                   "\nCombining and writing diagnostics.\n");
    int err_writediag = 0;
#ifdef TARGET
    diag_sum(&sim.diag_offload_data,
             diag_offload_array_mic0, diag_offload_array_mic1);
    mpi_gather_diag(&sim.diag_offload_data, diag_offload_array_mic0, n_tot,
                    mpi_rank, mpi_size, mpi_root);

    if(mpi_rank == mpi_root) {
        err_writediag = hdf5_interface_write_diagnostics(&sim,
            diag_offload_array_mic0, sim.hdf5_out);
    }
#else
    mpi_gather_diag(&sim.diag_offload_data, diag_offload_array_host, n_tot,
                    mpi_rank, mpi_size, mpi_root);

    if(mpi_rank == mpi_root) {
        err_writediag = hdf5_interface_write_diagnostics(&sim,
            diag_offload_array_host, sim.hdf5_out);
    }
#endif
    if(err_writediag) {
        print_out0(VERBOSE_MINIMAL, mpi_rank,
                   "\nWriting diagnostics failed.\n"
                   "See stderr for details.\n");
        /* Free offload data and terminate */
        goto CLEANUP_FAILURE;
    }
    else {
        print_out0(VERBOSE_MINIMAL, mpi_rank,
                   "Diagnostics written.\n");
    }

    /* Free offload data */

    mpi_interface_finalize();

#ifdef TARGET
    diag_free_offload(&sim.diag_offload_data, &diag_offload_array_mic0);
    diag_free_offload(&sim.diag_offload_data, &diag_offload_array_mic1);
#else
    diag_free_offload(&sim.diag_offload_data, &diag_offload_array_host);
#endif

    if(mpi_rank == mpi_root) {
        marker_summary(ps_gathered, n_gathered);
    }

    free(ps);
    free(ps_gathered);

    print_out0(VERBOSE_MINIMAL, mpi_rank, "\nDone.\n");

    return 0;


/* Free resources in case simulation crashes */
CLEANUP_FAILURE:

    mpi_interface_finalize();

#ifdef TARGET
    diag_free_offload(&sim.diag_offload_data, &diag_offload_array_mic0);
    diag_free_offload(&sim.diag_offload_data, &diag_offload_array_mic1);
#else
    diag_free_offload(&sim.diag_offload_data, &diag_offload_array_host);
#endif

    offload_free_offload(&offload_data, &offload_array);

    free(ps);
    free(ps_gathered);

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
                // The .hdf5 filename can be specified with or without the trailing .h5
                slen = strlen(optarg);
                if ( slen > 3 && !strcmp(optarg+slen-3,".h5") ) {
                    strncpy(sim->hdf5_in,optarg,slen-3);
                    (sim->hdf5_in)[slen-3]=NULL;
                }
                else
                    strcpy(sim->hdf5_in, optarg);
                break;
            case 2:
                // The .hdf5 filename can be specified with or without the trailing .h5
                slen = strlen(optarg);
                if ( slen > 3 && !strcmp(optarg+slen-3,".h5") ) {
                    strncpy(sim->hdf5_out,optarg,slen-3);
                    (sim->hdf5_out)[slen-3]=NULL;
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

/**
 * @brief Writes a summary of what happened to the markers during simulation
 *
 * This function writes a summary of marker end conditions and possible
 * simulation-time errors. Since simulation can have billions and billions of
 * markers, we only show how many markers had specific error or end condition.
 *
 * End conditions and errors are plotted in human-readable format.
 *
 * @param ps array of marker states after simulation has finished
 * @param n number of markers in the array
 */
void marker_summary(particle_state* ps, int n) {

    print_out(VERBOSE_MINIMAL, "\nSummary of results:\n");

    /* Temporary arrays that are needed to store unique end conditions and
     * errors. We can have at most n different values. */
    int* temp = (int*)malloc(n*sizeof(int));
    int* unique = (int*)malloc(n*sizeof(int));
    int* count = (int*)malloc((n+1)*sizeof(int));

    /* First we find out what the end conditions are */
    for(int i=0; i<n; i++) {
        temp[i] = ps[i].endcond;
    }
    math_uniquecount(temp, unique, count, n);
    count[n] = 0; // This ensures the following while loops are terminated.

    /* Print the end conditions, if marker has multiple end conditions, separate
     * those with "and". */
    int i = 0;
    while(count[i] > 0) {
        // Initialize
        int endconds[32];
        for(int j=0; j<32;j++) {
            endconds[j] = 0;
        }
        char endcondstr[256];
        endcondstr[0] = '\0';

        // Represent all end conditions with a single string and print it
        int j = 0;
        endcond_parse(unique[i], endconds);
        while(endconds[j]) {
            if(j>0) {
                strcat(endcondstr, " and ");
            }
            char temp[256];
            endcond_parse2str(endconds[j], temp);
            strcat(endcondstr, temp);
            j++;
        }
        if(j == 0) {
            sprintf(endcondstr, "Aborted");
        }
        print_out(VERBOSE_MINIMAL, "%9d markers had end condition %s\n",
                  count[i], endcondstr);
        i++;
    }

    // Empty line between end conditions and errors
    print_out(VERBOSE_MINIMAL, "\n");

    // Find all errors
    for(int i=0; i<n; i++) {
        temp[i] = (int)(ps[i].err);
    }
    math_uniquecount(temp, unique, count, n);

    // Go through each unique error and represent is a string
    i = 0;
    while(count[i] > 0) {
        if(unique[i] == 0) {
            // Value 0 indicates no error occurred so skip that
            i++;
            continue;
        }
        char msg[256];
        char line[256];
        char file[256];
        error_parse2str(unique[i], msg, line, file);
        print_out(VERBOSE_MINIMAL,
                  "%9d markers were aborted with an error message:\n"
                  "          %s\n"
                  "          at line %s in %s\n",
                  count[i], msg, line, file);
        i++;
    }

    // If count[0] equals to number of markers and their error field is zero,
    // we have no markers that were aborted.
    if(count[0] == n && unique[0] == 0) {
        print_out(VERBOSE_MINIMAL,
                  "          No markers were aborted.\n");
    }

    // Free temporary arrays
    free(temp);
    free(unique);
    free(count);
}
