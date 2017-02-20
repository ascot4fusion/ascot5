/**
 * @file ascot5_gc.c
 * @brief ASCOT5 with guiding center orbit following
 */
#include <getopt.h>
#include <math.h>
#ifdef MPI
  #include <mpi.h>
#endif
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "ascot5.h"
#include "wall.h"
#include "distributions.h"
#include "B_field.h"
#include "ascot4_interface.h"
#include "plasma_1d.h"
#include "interact.h"
#include "simulate.h"
#include "simulate_gc_rk4.h"
#include "particle.h"
#include "endcond.h"
#include "hdf5_histogram.h"

int read_options(int argc, char** argv, sim_offload_data* sim);

int main(int argc, char** argv) {
    /* Prepare simulation parameters and data for offload */
    sim_offload_data sim;
    real* B_offload_array;
    real* plasma_offload_array;
    real* wall_offload_array;
    real* dist_offload_array_mic0;
    real* dist_offload_array_mic1;
    real* dist_offload_array_host;

    /* Default values */
    sim.tstep = 1e-8;
    sim.tcollstep = 1e-6;
    sim.trstep = 1e-5;
    sim.tmax = 1.0;
    sim.active_endcond = endcond_tmax | endcond_emin | endcond_therm | endcond_wall | endcond_rhomax;
    sim.emin = 1e4*CONST_E;
    sim.dist_offload_data.n_r = 20;
    sim.dist_offload_data.min_r = 3;
    sim.dist_offload_data.max_r = 8.5;
    sim.dist_offload_data.n_z = 40;
    sim.dist_offload_data.min_z = -4.25;
    sim.dist_offload_data.max_z = 3.6;
    sim.dist_offload_data.n_vpara = 25;
    sim.dist_offload_data.min_vpara = -1.5e7;
    sim.dist_offload_data.max_vpara = 1.5e7;
    sim.dist_offload_data.n_vperp = 35;
    sim.dist_offload_data.min_vperp = 0;
    sim.dist_offload_data.max_vperp = 0;

    read_options(argc, argv, &sim);

    B_field_init_offload(&sim.B_offload_data, &B_offload_array);
    #if VERBOSE >= 1
    printf("Initialized magnetic field, %.1f MB.\n", sim.B_offload_data.offload_array_length * sizeof(real) / (1024.0*1024.0));
    #endif

    plasma_1d_init_offload(&sim.plasma_offload_data, &plasma_offload_array);
    #if VERBOSE >= 1
    printf("Initialized background plasma, %.1f MB.\n", sim.plasma_offload_data.offload_array_length * sizeof(real) / (1024.0*1024.0));
    #endif

    wall_init_offload(&sim.wall_offload_data, &wall_offload_array);
    #if VERBOSE >= 1
    printf("Initialized wall, %.1f MB.\n", sim.wall_offload_data.offload_array_length * sizeof(real) / (1024.0*1024.0));
    #endif

    #ifndef NOTARGET
    dist_rzvv_init_offload(&sim.dist_offload_data, &dist_offload_array_mic0);
    dist_rzvv_init_offload(&sim.dist_offload_data, &dist_offload_array_mic1);
    #else
    dist_rzvv_init_offload(&sim.dist_offload_data, &dist_offload_array_host);
    #endif
    #if VERBOSE >= 1
    printf("Initialized RZVV dist, %.1f MB.\n", sim.dist_offload_data.offload_array_length * sizeof(real) / (1024.0*1024.0));
    #endif

    int mpi_rank, mpi_size;
    #ifdef MPI
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    #else
    mpi_rank = 0;
    mpi_size = 1;
    #endif

    #if VERBOSE >= 1
    printf("Initialized MPI, rank %d, size %d.\n", mpi_rank, mpi_size);
    #endif

    int n;
    particle* p;
    ascot4_read_particles(&p, &n, "input.particles");

    char filename[256];
    if(mpi_size == 1) {
        sprintf(filename, "ascot.h5");
    }
    else {
        sprintf(filename, "ascot_%06d.h5", mpi_rank);
    }
    ascot4_write_inistate(n, p, filename);

    #if VERBOSE >= 1
    printf("Read %d particles.\n", n);
    #endif

    int start_index = mpi_rank * (n / mpi_size);
    p += start_index;

    if(mpi_rank == mpi_size-1) {
        n = n - mpi_rank * (n / mpi_size);
    }
    else {
        n = n / mpi_size;
    }

    #ifndef NOTARGET
    int n_mic = n / 2;
    int n_host = 0;
    #else
    int n_mic = 0;
    int n_host = n;
    #endif

    double mic0_start, mic0_end, mic1_start, mic1_end, host_start, host_end;

    fflush(stdout);
#ifdef _OMP
    omp_set_nested(1);
#endif
    #pragma omp parallel sections num_threads(3)
    {
    #ifndef NOTARGET
    #pragma omp section
    {
#ifdef _OMP
        mic0_start = omp_get_wtime();
#endif
        #pragma omp target device(0) map( \
            p[0:n_mic], \
            B_offload_array[0:sim.B_offload_data.offload_array_length], \
         dist_offload_array_mic0[0:sim.dist_offload_data.offload_array_length],\
         plasma_offload_array[0:sim.plasma_offload_data.offload_array_length], \
            wall_offload_array[0:sim.wall_offload_data.offload_array_length] \
         )
        simulate_gc_rk4(1, n_mic, p, sim, B_offload_array, plasma_offload_array,
                 wall_offload_array, dist_offload_array_mic0);
#ifdef _OMP
        mic0_end = omp_get_wtime();
#endif
    }

    #pragma omp section
    {
#ifdef _OMP
        mic1_start = omp_get_wtime();
#endif
        #pragma omp target device(1) map( \
            p[n_mic:n_mic], \
            B_offload_array[0:sim.B_offload_data.offload_array_length], \
         dist_offload_array_mic1[0:sim.dist_offload_data.offload_array_length],\
         plasma_offload_array[0:sim.plasma_offload_data.offload_array_length], \
            wall_offload_array[0:sim.wall_offload_data.offload_array_length] \
         )
        simulate_gc_rk4(2, n_mic, p+n_mic, sim, B_offload_array, plasma_offload_array,
                 wall_offload_array, dist_offload_array_mic1);
#ifdef _OMP
        mic1_end = omp_get_wtime();
#endif
    } 
    #else
    #pragma omp section
    {
#ifdef _OMP
        host_start = omp_get_wtime();
#endif
        simulate_gc_rk4(0, n_host, p+2*n_mic, sim, B_offload_array,
                 plasma_offload_array,
                 wall_offload_array,dist_offload_array_host);
#ifdef _OMP
	host_end = omp_get_wtime();
#endif
    }
    #endif
    }
    /* Code excution returns to host. */

    #if VERBOSE >= 1
        printf("mic0 %lf s, mic1 %lf s, host %lf s\n", mic0_end-mic0_start,
               mic1_end-mic1_start, host_end-host_start);
    #endif

    /* Combine histograms */
    #ifndef NOTARGET
    dist_rzvv_sum(&sim.dist_offload_data, dist_offload_array_mic0,
                  dist_offload_array_mic1);
    ascot4_write_dist_rzvv(&sim.dist_offload_data, dist_offload_array_mic0,
                           filename);
    #else
    ascot4_write_dist_rzvv(&sim.dist_offload_data, dist_offload_array_host,
                           filename);
    #endif
    ascot4_write_endstate(n, p, filename);

    B_field_free_offload(&sim.B_offload_data, &B_offload_array);
    plasma_1d_free_offload(&sim.plasma_offload_data, &plasma_offload_array);
    wall_free_offload(&sim.wall_offload_data, &wall_offload_array);
    #ifndef NOTARGET
    dist_rzvv_free_offload(&sim.dist_offload_data, &dist_offload_array_mic0);
    dist_rzvv_free_offload(&sim.dist_offload_data, &dist_offload_array_mic1);
    #else
    dist_rzvv_free_offload(&sim.dist_offload_data, &dist_offload_array_host);
    #endif
    free(p);

    printf("Done\n");
    return 0;
}

int read_options(int argc, char** argv, sim_offload_data* sim) {
    struct option longopts[] = {
        {"dt", required_argument, 0, 1},
        {"dtcoll", required_argument, 0, 2},
        {"dtoutput", required_argument, 0, 3},
        {"tmax", required_argument, 0, 4},
        {0, 0, 0, 0}
    };

    int c;
    while((c = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
        switch(c) {
            case 1:
                sim->tstep = atof(optarg);
                break;
            case 2:
                sim->tcollstep = atof(optarg);
                break;
            case 3:
                sim->trstep = atof(optarg);
                break;
            case 4:
                sim->tmax = atof(optarg);
                break;
        }
    }
    return 0;
}
