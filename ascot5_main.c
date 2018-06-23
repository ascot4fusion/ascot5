/**
 * @file ascot5_main.c
 * @brief ASCOT5
 */
#define _XOPEN_SOURCE 500
#include <getopt.h>
#include <math.h>
#ifdef MPI
  #include <mpi.h>
#endif
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ascot5.h"
#include "consts.h"
#include "wall.h"
#include "diag.h"
#include "B_field.h"
#include "plasma.h"
#include "print.h"
#include "simulate.h"
#include "particle.h"
#include "endcond.h"
#include "hdf5io/hdf5_diag.h"
#include "hdf5io/hdf5_input.h"
#include "hdf5io/hdf5_orbits.h"
#include "hdf5io/hdf5_particlestate.h"
#include "offload.h"

int read_options(int argc, char** argv, sim_offload_data* sim);
void generate_qid(char* qid);

int main(int argc, char** argv) {
    /* Prepare simulation parameters and data for offload */
    sim_offload_data sim;
    sim.mpi_rank = 0;
    sim.mpi_size = 1;

    read_options(argc, argv, &sim);

    /* Get MPI rank and set qid for the run.
     * qid rules: The actual random unique qid is used in MPI or single-process
     * runs. If this is a multi-process run (user-defined MPI rank and size),
     * e.g. condor run, we set qid = 5 000 000 000 since 32 bit integers don't
     * go that high. The actual qid is assigned when results are combined. */
    int mpi_rank, mpi_size;
    char qid[] = "5000000000";
#ifdef MPI
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if(sim.mpi_size == 0) {
        /* Let MPI determine size and rank */
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        sim.mpi_rank = mpi_rank;
        sim.mpi_size = mpi_size;
        generate_qid(qid);
    }
    else {
        /* Use user-defined size and rank */
        mpi_rank = sim.mpi_rank;
        mpi_size = sim.mpi_size;
        if(mpi_size == 1) {
            generate_qid(qid);
        }
    }
#else
    if(sim.mpi_size == 0) {
        /* Use default values (single process) */
        mpi_rank = 0;
        mpi_size = 1;
        generate_qid(qid);
    }
    else {
        /* Use user-defined size and rank */
        mpi_rank = sim.mpi_rank;
        mpi_size = sim.mpi_size;
        if(mpi_size == 1) {
            generate_qid(qid);
        }
    }
#endif

    print_out0(VERBOSE_NORMAL, mpi_rank, "Initialized MPI, rank %d, size %d.\n", mpi_rank, mpi_size);

    int err = 0;
    int n;
    input_particle* p;
    real* B_offload_array;
    real* E_offload_array;
    real* plasma_offload_array;
    real* neutral_offload_array;
    real* wall_offload_array;
    err = hdf5_input(&sim, &B_offload_array, &E_offload_array, &plasma_offload_array, 
                     &neutral_offload_array, &wall_offload_array, &p, &n);
    if(err) {return 0;};

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

#ifdef TARGET
    real* diag_offload_array_mic0;
    real* diag_offload_array_mic1;
    diag_init_offload(&sim.diag_offload_data, &diag_offload_array_mic0);
    diag_init_offload(&sim.diag_offload_data, &diag_offload_array_mic1);
#else
    real* diag_offload_array_host;
    diag_init_offload(&sim.diag_offload_data, &diag_offload_array_host);
#endif
    
    print_out0(VERBOSE_NORMAL, mpi_rank, "Initialized diagnostics, %.1f MB.\n", sim.diag_offload_data.offload_array_length * sizeof(real) / (1024.0*1024.0));

    /* Set output filename for this MPI process. */
    if(mpi_size == 1) {
	strcat(sim.hdf5_out, ".h5");
    }
    else {
        char temp[256];
        sprintf(temp, "_%06d.h5", mpi_rank);
	strcat(sim.hdf5_out, temp);
    }

    err = hdf5_initoutput(&sim, qid);
    if(err) {return 0;};
    strcpy(sim.qid, qid);

    int start_index = mpi_rank * (n / mpi_size);
    p += start_index;

    if(mpi_rank == mpi_size-1) {
        n = n - mpi_rank * (n / mpi_size);
    }
    else {
        n = n / mpi_size;
    }

    /* Set up particlestates on host, needs magnetic field evaluation */
    B_field_data Bdata;
    B_field_init(&Bdata, &sim.B_offload_data, B_offload_array);

    print_out0(VERBOSE_NORMAL, mpi_rank, "Magnetic field initialization complete.\n");

    particle_state* ps = (particle_state*) malloc(n * sizeof(particle_state));
    for(int i = 0; i < n; i++) {
        particle_input_to_state(&p[i], &ps[i], &Bdata);
    }
    
    hdf5_particlestate_write(sim.hdf5_out, qid, "inistate", n, ps);

    print_out0(VERBOSE_NORMAL, mpi_rank, "Markers initialized and inistate written.\n");

    #ifdef TARGET
        int n_mic = n / TARGET;
        int n_host = 0;
    #else
        int n_mic = 0;
        int n_host = n;
    #endif

    double mic0_start = 0, mic0_end=0, mic1_start=0, mic1_end=0, host_start=0, host_end=0;
    
    fflush(stdout);

    omp_set_nested(1);
    
    #pragma omp parallel sections num_threads(3)
    {
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
    /* Code excution returns to host. */

    print_out0(VERBOSE_NORMAL, mpi_rank, "Writing endstate.");

    hdf5_particlestate_write(sim.hdf5_out, qid, "endstate", n, ps);

    print_out0(VERBOSE_NORMAL, mpi_rank, "mic0 %lf s, mic1 %lf s, host %lf s\n",
        mic0_end-mic0_start, mic1_end-mic1_start, host_end-host_start);
    
    /* Combine histograms */
    #ifdef TARGET
        diag_sum(&sim.diag_offload_data, diag_offload_array_mic0,diag_offload_array_mic1);
        hdf5_diag_write(&sim, diag_offload_array_mic0, sim.hdf5_out, qid);
    #else
        hdf5_diag_write(&sim, diag_offload_array_host, sim.hdf5_out, qid);
    #endif
    

    #ifdef MPI
        MPI_Finalize();
    #endif

    B_field_free_offload(&sim.B_offload_data, &B_offload_array);
    plasma_free_offload(&sim.plasma_offload_data, &plasma_offload_array);
    wall_free_offload(&sim.wall_offload_data, &wall_offload_array);
    #ifdef TARGET
        diag_free_offload(&sim.diag_offload_data, &diag_offload_array_mic0);
	diag_free_offload(&sim.diag_offload_data, &diag_offload_array_mic1);
    #else
        diag_free_offload(&sim.diag_offload_data, &diag_offload_array_host);
    #endif
    offload_free_offload(&offload_data, &offload_array);
    
    free(p-start_index);

    print_out0(VERBOSE_MINIMAL, mpi_rank, "Done.\n");

    return 0;
}

int read_options(int argc, char** argv, sim_offload_data* sim) {
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
            printf("\nUnrecognized option. The valid parameters are:\n");
            printf("--in input file without .h5 (default: ascot)\n");
            printf("--out output file without .h5 (default: same as input)\n");
            printf("--mpi_size number of independent processes\n");
            printf("--mpi_rank rank of independent process\n");
            abort();
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

/** @brief Generate an identification muber 
 *  
 *  qid is a 32 bit unsigned integer, which is represented
 *  in string format. The string is formed by 10 numbers and it
 *  is padded with leading zeroes.
 */
void generate_qid(char* qid) {
    /* Generate 32 bit random integer. */
    srand48( time(NULL) );
    long int qint = -1;
    while(qint < 0) {
        qint = mrand48();
    }
    
    /* Turn it into a string */
    sprintf(qid, "%010lu", (long unsigned int)qint);
}
