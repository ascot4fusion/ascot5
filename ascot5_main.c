/**
 * @file ascot5_main.c
 * @brief ASCOT5
 */
#include <getopt.h>
#include <math.h>
#ifdef MPI
  #include <mpi.h>
#endif
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ascot5.h"
#include "consts.h"
#include "wall.h"
#include "distributions.h"
#include "diag.h"
#include "B_field.h"
#include "ascot4_interface.h"
#include "plasma_1d.h"
#include "interact.h"
#include "simulate.h"
#include "particle.h"
#include "endcond.h"
#include "hdf5io/hdf5_diag.h"
#include "hdf5io/hdf5_input.h"
#include "hdf5io/hdf5_orbits.h"
#include "hdf5io/hdf5_particlestate.h"

int read_options(int argc, char** argv, sim_offload_data* sim);

int main(int argc, char** argv) {
    /* Prepare simulation parameters and data for offload */
    sim_offload_data sim;
    real* B_offload_array;
    real* E_offload_array;
    real* plasma_offload_array;
    real* wall_offload_array;
    real* diag_offload_array_mic0;
    real* diag_offload_array_mic1;
    real* diag_offload_array_host;
    int n;
    input_particle* p;

    int err = 0;

    read_options(argc, argv, &sim);
    
    err = hdf5_input(&sim, &B_offload_array, &E_offload_array, &plasma_offload_array, 
		     &wall_offload_array, &p, &n);
    if(err) {return 0;};
    
    #ifndef NOTARGET
        diag_init_offload(&sim.diag_offload_data, &diag_offload_array_mic0);
	diag_init_offload(&sim.diag_offload_data, &diag_offload_array_mic1);
    #else
	diag_init_offload(&sim.diag_offload_data, &diag_offload_array_host);
    #endif
	
    #if VERBOSE >= 1
    printf("Initialized diagnostics, %.1f MB.\n", sim.diag_offload_data.offload_array_length * sizeof(real) / (1024.0*1024.0));
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

    if(mpi_size == 1) {
	char temp[256];
	strcat(sim.hdf5_out, ".h5");
    }
    else {
	char temp[256];
        sprintf(temp, "_%06d.h5", mpi_rank);
	strcat(sim.hdf5_out, temp);
    }
    err = hdf5_checkoutput(&sim);
    if(err) {return 0;};

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
        E_offload_array[0:sim.E_offload_data.offload_array_length], \
        diag_offload_array_mic0[0:sim.diag_offload_data.offload_array_length], \
        plasma_offload_array[0:sim.plasma_offload_data.offload_array_length], \
        wall_offload_array[0:sim.wall_offload_data.offload_array_length] \
                )
                simulate(1, n_mic, p, &sim, B_offload_array, E_offload_array,
                         plasma_offload_array, wall_offload_array,
                         diag_offload_array_mic0);

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
        p[n_mic:2*n_mic], \
        B_offload_array[0:sim.B_offload_data.offload_array_length], \
        E_offload_array[0:sim.E_offload_data.offload_array_length], \
        diag_offload_array_mic1[0:sim.diag_offload_data.offload_array_length], \
        plasma_offload_array[0:sim.plasma_offload_data.offload_array_length], \
        wall_offload_array[0:sim.wall_offload_data.offload_array_length] \
                )
                simulate(2, n_mic, p+n_mic, &sim, B_offload_array,
                         E_offload_array, plasma_offload_array,
                         wall_offload_array, diag_offload_array_mic1);

                #ifdef _OMP
                mic1_end = omp_get_wtime();
                #endif
            }

        #endif
            #pragma omp section
            {
                #ifdef _OMP
                host_start = omp_get_wtime();
                #endif
		
		sim_data sim_host;
        
                simulate_begin(0, n_host, p+2*n_mic, &sim, &sim_host,
		    B_offload_array,
		    E_offload_array,
		    plasma_offload_array, wall_offload_array,
		    diag_offload_array_host);

		hdf5_particlestate_write(sim.hdf5_out, "inistate", n, p);
		
		simulate_continue(0, n_host, p+2*n_mic, &sim_host);
		
		hdf5_particlestate_write(sim.hdf5_out, "endstate", n, p);

		hdf5_orbits_write(&sim_host, sim.hdf5_out);

		simulate_end(&sim_host);

		#ifdef _OMP
                host_end = omp_get_wtime();
                #endif
            }
    }
    /* Code excution returns to host. */
    
    #if VERBOSE >= 1
    printf("mic0 %lf s, mic1 %lf s, host %lf s\n", mic0_end-mic0_start,
           mic1_end-mic1_start, host_end-host_start);
    #endif
    
    /* Combine histograms */
    #ifndef NOTARGET
        diag_sum(&sim.diag_offload_data, diag_offload_array_mic0,diag_offload_array_mic1);
        hdf5_diag_write(&sim, diag_offload_array_mic0, sim.hdf5_out);
    #else
	hdf5_diag_write(&sim, diag_offload_array_host, sim.hdf5_out);
    #endif
    

    B_field_free_offload(&sim.B_offload_data, &B_offload_array);
    plasma_1d_free_offload(&sim.plasma_offload_data, &plasma_offload_array);
    wall_free_offload(&sim.wall_offload_data, &wall_offload_array);
    #ifndef NOTARGET
        diag_free_offload(&sim.diag_offload_data, &diag_offload_array_mic0);
        diag_free_offload(&sim.diag_offload_data, &diag_offload_array_mic1);
    #else
        diag_free_offload(&sim.diag_offload_data, &diag_offload_array_host);
    #endif
    free(p);

    printf("Done\n");
    return 0;
}

int read_options(int argc, char** argv, sim_offload_data* sim) {
    struct option longopts[] = {
        {"in", required_argument, 0, 1},
        {"out", required_argument, 0, 2},
        {0, 0, 0, 0}
    };

    sim->hdf5_in[0]  = '\0';
    sim->hdf5_out[0] = '\0';

    int c;
    while((c = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
        switch(c) {
	    case 1:
	        strcpy(sim->hdf5_in, optarg);
	        break;
	    case 2:
                strcpy(sim->hdf5_out, optarg);
                break;
	    default:
		printf("\nUnrecognized option. The valid parameters are:\n");
		printf("-in hdf5 input file (default ascot)\n");
		printf("-out hdf5 output file (default same as input)\n");
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
    
    return 0;
}
