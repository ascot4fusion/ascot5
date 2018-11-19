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
 * In addition to output data, the simulation progress may be written in
 * *.stdout files with each MPI process having dedicated file. See ascot5.h for
 * details.
 */
#define _XOPEN_SOURCE 500
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

int read_arguments(int argc, char** argv, sim_offload_data* sim);
void generate_qid(char* qid);
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

    /* Read and parse command line arguments */
    sim_offload_data sim;
    if( read_arguments(argc, argv, &sim) ) {
        abort();
        return 1;
    }

    /* Get MPI rank and set qid for the run.
     * qid rules: The actual random unique qid is used in MPI or single-process
     * runs. If this is a multi-process run (user-defined MPI rank and size),
     * e.g. condor run, we set qid = 5 000 000 000 since 32 bit integers don't
     * go that high. The actual qid is assigned when results are combined.*/
    int mpi_rank, mpi_size;
    char qid[] = "5000000000";

    if(sim.mpi_size == 0) {
#ifdef MPI
        /* MPI run */
        int provided;
        MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        sim.mpi_rank = mpi_rank;
        sim.mpi_size = mpi_size;
        generate_qid(qid);
#else
        /* MPI was not included while compiling       */
        /* Give warning  and run a single process run */
        mpi_rank = 0;
        mpi_size = 1;
        generate_qid(qid);
#endif
    }
    else {
        /* Emulate MPI run (Condor-like run) */
        /* Use user-defined size and rank    */
        mpi_rank = sim.mpi_rank;
        mpi_size = sim.mpi_size;
    }

    print_out0(VERBOSE_NORMAL, mpi_rank,
               "Initialized MPI, rank %d, size %d.\n", mpi_rank, mpi_size);

    /* Number of markers to be simulated */
    int n;
    /* Marker input struct */
    input_particle* p;

    /* Offload data arrays that are allocated when input is read */
    real* B_offload_array;
    real* E_offload_array;
    real* plasma_offload_array;
    real* neutral_offload_array;
    real* wall_offload_array;

    /* Read input from the HDF5 file */
    if( hdf5_interface_read_input(&sim, &B_offload_array, &E_offload_array,
                                  &plasma_offload_array, &neutral_offload_array,
                                  &wall_offload_array, &p, &n) ) {
        print_out0(VERBOSE_MINIMAL, mpi_rank,
                   "\nInput reading or initializing failed.\n"
                   "See stderr for details.\n");
        abort();
        return 1;
    };

    /* Pack offload data into single array */
    real* offload_array;
    offload_package offload_data;
    offload_init_offload(&offload_data, &offload_array);
    offload_pack(&offload_data, &offload_array, B_offload_array,
                 sim.B_offload_data.offload_array_length);
    offload_pack(&offload_data, &offload_array, E_offload_array,
                 sim.E_offload_data.offload_array_length);
    offload_pack(&offload_data, &offload_array, plasma_offload_array,
                 sim.plasma_offload_data.offload_array_length);
    offload_pack(&offload_data, &offload_array, neutral_offload_array,
                 sim.neutral_offload_data.offload_array_length);
    offload_pack(&offload_data, &offload_array, wall_offload_array,
                 sim.wall_offload_data.offload_array_length);

    /* Initialize diagnostics offload data.
     * Separate arrays for host and target */
#ifdef TARGET
    real* diag_offload_array_mic0;
    real* diag_offload_array_mic1;
    diag_init_offload(&sim.diag_offload_data, &diag_offload_array_mic0);
    diag_init_offload(&sim.diag_offload_data, &diag_offload_array_mic1);
#else
    real* diag_offload_array_host;
    diag_init_offload(&sim.diag_offload_data, &diag_offload_array_host);
#endif

    real diag_offload_array_size = sim.diag_offload_data.offload_array_length
        * sizeof(real) / (1024.0*1024.0);
    print_out0(VERBOSE_NORMAL, mpi_rank,
               "Initialized diagnostics, %.1f MB.\n", diag_offload_array_size);

    /* Set output filename for this MPI process. */
    if(mpi_size == 1) {
        strcat(sim.hdf5_out, ".h5");
    }
    else {
        char temp[256];
        sprintf(temp, "_%06d.h5", mpi_rank);
        strcat(sim.hdf5_out, temp);
    }

    /* Choose which markers are used in this MPI process. Simply put, markers
     * are divided into mpi_size sequential blocks and the mpi_rank:th block
     * is chosen for this simulation. */
    int start_index = mpi_rank * (n / mpi_size);
    p += start_index;

    if(mpi_rank == mpi_size-1) {
        n = n - mpi_rank * (n / mpi_size);
    }
    else {
        n = n / mpi_size;
    }

    /* Set up particlestates on host, needs magnetic field evaluation */
    print_out0(VERBOSE_NORMAL, mpi_rank,
               "\nInitializing marker states.\n");
    B_field_data Bdata;
    B_field_init(&Bdata, &sim.B_offload_data, B_offload_array);
    particle_state* ps = (particle_state*) malloc(n * sizeof(particle_state));
    for(int i = 0; i < n; i++) {
        particle_input_to_state(&p[i], &ps[i], &Bdata);
    }
    free(p-start_index); // Input markers are no longer required
    print_out0(VERBOSE_NORMAL, mpi_rank,
               "Estimated memory usage %.1f MB.\n",
               (sizeof(real) * n) / (1024.0*1024.0));
    print_out0(VERBOSE_NORMAL, mpi_rank,
               "Marker states initialized.\n");

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

    /* Write inistate */
    if( hdf5_interface_write_state(sim.hdf5_out, "inistate", n, ps) ) {
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

    /* Write endstate */
    if( hdf5_interface_write_state(sim.hdf5_out, "endstate", n, ps) ) {
        print_out0(VERBOSE_MINIMAL, mpi_rank,
                   "\nWriting endstate failed.\n"
                   "See stderr for details.\n");
        /* Free offload data and terminate */
        goto CLEANUP_FAILURE;
    }
    print_out0(VERBOSE_NORMAL, mpi_rank,
               "Endstate written.\n");

    /* Combine diagnostic data and write it to HDF5 file */
    print_out0(VERBOSE_MINIMAL, mpi_rank,
                   "\nCombining and writing diagnostics.\n");
    int err_writediag = 0;
#ifdef TARGET
    diag_sum(&sim.diag_offload_data,
             diag_offload_array_mic0, diag_offload_array_mic1);
    err_writediag = hdf5_interface_write_diagnostics(
        &sim, diag_offload_array_mic0, sim.hdf5_out);
#else
    err_writediag = hdf5_interface_write_diagnostics(
        &sim, diag_offload_array_host, sim.hdf5_out);
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
    goto CLEANUP_SUCCESS;

CLEANUP_SUCCESS:

#ifdef MPI
    MPI_Finalize();
#endif

    B_field_free_offload(&sim.B_offload_data, &B_offload_array);
    E_field_free_offload(&sim.E_offload_data, &E_offload_array);
    plasma_free_offload(&sim.plasma_offload_data, &plasma_offload_array);
    wall_free_offload(&sim.wall_offload_data, &wall_offload_array);
    neutral_free_offload(&sim.neutral_offload_data, &neutral_offload_array);

#ifdef TARGET
    diag_free_offload(&sim.diag_offload_data, &diag_offload_array_mic0);
    diag_free_offload(&sim.diag_offload_data, &diag_offload_array_mic1);
#else
    diag_free_offload(&sim.diag_offload_data, &diag_offload_array_host);
#endif

    offload_free_offload(&offload_data, &offload_array);

    marker_summary(ps, n);
    free(ps);

    print_out0(VERBOSE_MINIMAL, mpi_rank, "Done.\n");

    return 0;

CLEANUP_FAILURE:

#ifdef MPI
    MPI_Finalize();
#endif

    B_field_free_offload(&sim.B_offload_data, &B_offload_array);
    E_field_free_offload(&sim.E_offload_data, &E_offload_array);
    plasma_free_offload(&sim.plasma_offload_data, &plasma_offload_array);
    wall_free_offload(&sim.wall_offload_data, &wall_offload_array);
    neutral_free_offload(&sim.neutral_offload_data, &neutral_offload_array);

#ifdef TARGET
    diag_free_offload(&sim.diag_offload_data, &diag_offload_array_mic0);
    diag_free_offload(&sim.diag_offload_data, &diag_offload_array_mic1);
#else
    diag_free_offload(&sim.diag_offload_data, &diag_offload_array_host);
#endif

    offload_free_offload(&offload_data, &offload_array);

    free(ps);

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
 * - sim->hdf5_in  = "in.h5" ("ascot.h5")
 * - sim->hdf5_out = "out.h5" (sim->hdf5_in is copied here)
 * - sim->mpi_rank = 0
 * - sim->mpi_size = 0
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
        {0, 0, 0, 0}
    };

    sim->hdf5_in[0]  = '\0';
    sim->hdf5_out[0] = '\0';
    sim->mpi_rank = 0;
    sim->mpi_size = 0;

    int c;
    while((c = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
        switch(c) {
        case 1:
            strcpy(sim->hdf5_in, optarg);
            break;
        case 2:
            strcpy(sim->hdf5_out, optarg);
            break;
        case 3:
            sim->mpi_size = atoi(optarg);
            break;
        case 4:
            sim->mpi_rank = atoi(optarg);
            break;
        default:
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
            return 1;
        }
    }

    if(sim->hdf5_in[0] == '\0' && sim->hdf5_out[0] == '\0') {
        strcpy(sim->hdf5_in, "ascot.h5");
        strcpy(sim->hdf5_out, "ascot");
    }
    else if(sim->hdf5_in[0] == '\0' && sim->hdf5_out[0] != '\0') {
        strcpy(sim->hdf5_in, "ascot.h5");
    }
    else if(sim->hdf5_in[0] != '\0' && sim->hdf5_out[0] == '\0') {
        strcpy(sim->hdf5_out, sim->hdf5_in);
        strcat(sim->hdf5_in, ".h5");
    }
    else {
        strcat(sim->hdf5_in, ".h5");
    }
    strcpy(sim->outfn, sim->hdf5_out);
    return 0;
}

/**
 * @brief Generate an identification number for a run
 *
 * The identification number (QID) is a 32 bit unsigned integer represented in a
 * string format, i.e., by ten characters. QID is a random integer between 0 and
 * 4 294 967 295, and it is padded with leading zeroes in string representation.
 *
 * @param a pointer to 11 chars wide array where generated QID is stored
 */
void generate_qid(char* qid) {

    /* Seed random number generator with current time */
    srand48( time(NULL) );

    /* Generate a 32 bit random integer by generating signed 32 bit random
     * integers with mrand48() and choosing the first one that is positive */
    long int qint = -1;
    while(qint < 0) {
        qint = mrand48();
    }

    /* Convert the random number to a string format */
    sprintf(qid, "%010lu", (long unsigned int)qint);
}

void marker_summary(particle_state* ps, int n) {

    print_out(VERBOSE_MINIMAL, "\nSummary of results:\n")

    int* temp = (int*)malloc(n*sizeof(int));
    int* unique = (int*)malloc(n*sizeof(int));
    int* count = (int*)malloc(n*sizeof(int));

    for(int i=0; i<n; i++) {
        temp[i] = ps[i].endcond;
    }
    math_uniquecount(temp, unique, count, n);

    int i = 0;
    while(count[i]) {
        print_out(VERBOSE_MINIMAL, "%9d markers had end condition %d\n",
                  count[i], unique[i]);
        i++;
    }

    for(int i=0; i<n; i++) {
        temp[i] = (int)(ps[i].err);
    }
    math_uniquecount(temp, unique, count, n);

    print_out(VERBOSE_MINIMAL, "\n");
    i = 0;
    while(unique[i] > 0) {
        print_out(VERBOSE_MINIMAL, "%9d markers were aborted with error %d\n",
                  count[i], unique[i]);
        i++;
    }
    if(count[0] == n) {
        print_out(VERBOSE_MINIMAL, "No markers were aborted.\n");
    }

    free(temp);
    free(unique);
    free(count);
}
