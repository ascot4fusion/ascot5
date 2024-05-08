/**
 * @file mpi_interface.c
 * @brief MPI interface
 *
 * This module provides interfaces for communication with MPI.
 */
#ifdef MPI
#include <mpi.h>
#endif
#include <stddef.h>
#include <stdlib.h>
#include "ascot5.h"
#include "diag.h"
#include "mpi_interface.h"
#include "particle.h"
#include "simulate.h"

/**
 * @brief Initialize MPI
 *
 * This function initializes MPI and sets the mpi_rank and mpi_size accordingly.
 * If compiled without MPI, mpirank is set to zero and mpi_size is set to one.
 *
 * @param argc count of the command line arguments
 * @param argv pointers to the command line arguments
 * @param sim pointer to simulation offload struct
 * @param mpi_rank pointer to the determined rank of this MPI process
 * @param mpi_size pointer to the determined number of MPI processes
 * @param mpi_root pointer to the determined rank of the root MPI process
 */
void mpi_interface_init(int argc, char** argv, int* mpi_rank, int* mpi_size,
                        int* mpi_root) {
#ifdef MPI
    /* MPI is used. Init MPI and set rank and size as determined by MPI.
     * Rank of the root process is zero */
    int provided;
    int initialized;
    MPI_Initialized(&initialized);
    if(initialized == 0)
        MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, mpi_size);
    *mpi_root = 0;
#else
    /* MPI is not used and no command line arguments are given. Assume that
     * the entire simulation is a single process. */
    *mpi_rank = 0;
    *mpi_size = 1;
    *mpi_root = 0;
#endif
}

/**
 * @brief Finalize MPI
 *
 * This function finalizes the MPI environment, to be called at the end of
 * execution.
 */
void mpi_interface_finalize() {
#ifdef MPI
    MPI_Finalize();
#endif
}

/**
 * @brief Divide markers to mpi processes
 *
 * This function divides ntotal markers evenly and returns the starting index
 * and number of markers to be simulated for each process.
 *
 * @param start_index pointer to variable for starting index in input markers
 * @param n pointer to variable for number of markers for this process
 * @param n_tot number of markers in the simulation
 * @param mpi_rank rank of this MPI process
 * @param mpi_size total number of MPI processes
 */
void mpi_my_particles(int* start_index, int* n, int n_tot, int mpi_rank,
                      int mpi_size) {
    if(mpi_rank == mpi_size-1) {
        *n = n_tot - mpi_rank * (n_tot / mpi_size);
    }
    else {
        *n = n_tot / mpi_size;
    }

    *start_index = mpi_rank * (n_tot / mpi_size);
}

/**
 * @brief Gather all particle states to the root process
 *
 * This function gathers the particle states from each process to an array
 * in the root process. In condor-style execution only particle states for
 * current process are stored. An allocated array for gathered particle states
 * is stored in psgathered and number of states in ngathered.
 *
 * @todo Could be done more cleanly with custom datatypes
 *
 * @param ps pointer to array particle states for this process
 * @param ps_gather pointer to pointer to array where markers are gathered
 * @param n_gather pointer to variable for number of gathered markers
 * @param n_tot total number of markers in the simulation
 * @param mpi_rank rank of this MPI process
 * @param mpi_size total number of MPI processes
 * @param mpi_root rank of the root process
 */
void mpi_gather_particlestate(
    particle_state* ps, particle_state** ps_gather, int* n_gather, int n_tot,
    int mpi_rank, int mpi_size, int mpi_root) {
#ifdef MPI
    const int n_real = 32;
    const int n_int = 5;
    const int n_err = 1;

    particle_state* ps_all = malloc(n_tot * sizeof(particle_state));

    if(mpi_rank == mpi_root) {
        int start_index, n;
        mpi_my_particles(&start_index, &n, n_tot, mpi_rank, mpi_size);

        for(int j = 0; j < n; j++) {
            ps_all[j] = ps[j];
        }

        for(int i = 1; i < mpi_size; i++) {
            mpi_my_particles(&start_index, &n, n_tot, i, mpi_size);

            real* realdata;
            realdata = malloc(n_real * n * sizeof(realdata));
            integer* intdata;
            intdata = malloc(n_int * n * sizeof(intdata));
            a5err* errdata;
            errdata = malloc(n_err * n * sizeof(intdata));

            MPI_Recv(realdata, n_real*n, mpi_type_real, i, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(intdata, n_int*n, mpi_type_integer, i, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(errdata, n_err*n, mpi_type_a5err, i, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for(int j = 0; j < n; j++) {
                ps_all[start_index+j].r          = realdata[0*n+j];
                ps_all[start_index+j].phi        = realdata[1*n+j];
                ps_all[start_index+j].z          = realdata[2*n+j];
                ps_all[start_index+j].ppar       = realdata[3*n+j];
                ps_all[start_index+j].mu         = realdata[4*n+j];
                ps_all[start_index+j].zeta       = realdata[5*n+j];
                ps_all[start_index+j].rprt       = realdata[6*n+j];
                ps_all[start_index+j].phiprt     = realdata[7*n+j];
                ps_all[start_index+j].zprt       = realdata[8*n+j];
                ps_all[start_index+j].p_r        = realdata[9*n+j];
                ps_all[start_index+j].p_phi      = realdata[10*n+j];
                ps_all[start_index+j].p_z        = realdata[11*n+j];
                ps_all[start_index+j].mass       = realdata[12*n+j];
                ps_all[start_index+j].charge     = realdata[13*n+j];
                ps_all[start_index+j].anum       = intdata[0*n+j];
                ps_all[start_index+j].znum       = intdata[1*n+j];
                ps_all[start_index+j].weight     = realdata[14*n+j];
                ps_all[start_index+j].time       = realdata[15*n+j];
                ps_all[start_index+j].cputime    = realdata[16*n+j];
                ps_all[start_index+j].rho        = realdata[17*n+j];
                ps_all[start_index+j].theta      = realdata[18*n+j];
                ps_all[start_index+j].id         = intdata[2*n+j];
                ps_all[start_index+j].endcond    = intdata[3*n+j];
                ps_all[start_index+j].walltile   = intdata[4*n+j];
                ps_all[start_index+j].B_r        = realdata[19*n+j];
                ps_all[start_index+j].B_phi      = realdata[20*n+j];
                ps_all[start_index+j].B_z        = realdata[21*n+j];
                ps_all[start_index+j].B_r_dr     = realdata[22*n+j];
                ps_all[start_index+j].B_phi_dr   = realdata[23*n+j];
                ps_all[start_index+j].B_z_dr     = realdata[24*n+j];
                ps_all[start_index+j].B_r_dphi   = realdata[25*n+j];
                ps_all[start_index+j].B_phi_dphi = realdata[26*n+j];
                ps_all[start_index+j].B_z_dphi   = realdata[27*n+j];
                ps_all[start_index+j].B_r_dz     = realdata[28*n+j];
                ps_all[start_index+j].B_phi_dz   = realdata[29*n+j];
                ps_all[start_index+j].B_z_dz     = realdata[30*n+j];
                ps_all[start_index+j].mileage    = realdata[31*n+j];
                ps_all[start_index+j].err        = errdata[j];
            }

            free(realdata);
            free(intdata);
            free(errdata);
        }
    }
    else {

        int start_index, n;
        mpi_my_particles(&start_index, &n, n_tot, mpi_rank, mpi_size);

        real* realdata;
        realdata = malloc(n_real * n * sizeof(realdata));
        integer* intdata;
        intdata = malloc(n_int * n * sizeof(intdata));
        a5err* errdata;
        errdata = malloc(n_err * n * sizeof(intdata));

        for(int j = 0; j < n; j++) {
            realdata[0*n+j]  = ps[j].r;
            realdata[1*n+j]  = ps[j].phi;
            realdata[2*n+j]  = ps[j].z;
            realdata[3*n+j]  = ps[j].ppar;
            realdata[4*n+j]  = ps[j].mu;
            realdata[5*n+j]  = ps[j].zeta;
            realdata[6*n+j]  = ps[j].rprt;
            realdata[7*n+j]  = ps[j].phiprt;
            realdata[8*n+j]  = ps[j].zprt;
            realdata[9*n+j]  = ps[j].p_r;
            realdata[10*n+j] = ps[j].p_phi;
            realdata[11*n+j] = ps[j].p_z;
            realdata[12*n+j] = ps[j].mass;
            realdata[13*n+j] = ps[j].charge;
            intdata[0*n+j]   = ps[j].anum;
            intdata[1*n+j]   = ps[j].znum;
            realdata[14*n+j] = ps[j].weight;
            realdata[15*n+j] = ps[j].time;
            realdata[16*n+j] = ps[j].cputime;
            realdata[17*n+j] = ps[j].rho;
            realdata[18*n+j] = ps[j].theta;
            intdata[2*n+j]   = ps[j].id;
            intdata[3*n+j]   = ps[j].endcond;
            intdata[4*n+j]   = ps[j].walltile;
            realdata[19*n+j] = ps[j].B_r;
            realdata[20*n+j] = ps[j].B_phi;
            realdata[21*n+j] = ps[j].B_z;
            realdata[22*n+j] = ps[j].B_r_dr;
            realdata[23*n+j] = ps[j].B_phi_dr;
            realdata[24*n+j] = ps[j].B_z_dr;
            realdata[25*n+j] = ps[j].B_r_dphi;
            realdata[26*n+j] = ps[j].B_phi_dphi;
            realdata[27*n+j] = ps[j].B_z_dphi;
            realdata[28*n+j] = ps[j].B_r_dz;
            realdata[29*n+j] = ps[j].B_phi_dz;
            realdata[30*n+j] = ps[j].B_z_dz;
            realdata[31*n+j] = ps[j].mileage;
            errdata[j] = ps[j].err;
        }

        MPI_Send(realdata, n_real*n, mpi_type_real, 0, 0, MPI_COMM_WORLD);
        MPI_Send(intdata, n_int*n, mpi_type_integer, 0, 0, MPI_COMM_WORLD);
        MPI_Send(errdata, n_err*n, mpi_type_a5err, 0, 0, MPI_COMM_WORLD);

        free(realdata);
        free(intdata);
        free(errdata);
    }

    *ps_gather = ps_all;
    *n_gather = n_tot;

#else

    int start_index, n;
    mpi_my_particles(&start_index, &n, n_tot, mpi_rank, mpi_size);

    /* Only store particles for this process in Condor-style execution */
    particle_state* ps_all = malloc(n * sizeof(particle_state));

    for(int j = 0; j < n; j++) {
        ps_all[j] = ps[j];
    }

    *ps_gather = ps_all;
    *n_gather = n;

#endif
}

/**
 * @brief Gather all diagnostics to the root process
 *
 * This function gathers the distributions and orbits to the root process.
 * Distributions are summed and orbit data is appended to the root process
 * diagnostics array.
 *
 * @param data diagnostics offload data
 * @param offload_array pointer to diagnostics offload array
 * @param ntotal total number of markers in the simulation
 * @param mpi_rank rank of this MPI process
 * @param mpi_size total number of MPI processes
 * @param mpi_root rank of the root process
 */
void mpi_gather_diag(diag_offload_data* data, real* offload_array, int ntotal,
                     int mpi_rank, int mpi_size, int mpi_root) {
#ifdef MPI

    if(data->dist5D_collect || data->distrho5D_collect
       || data->dist6D_collect || data->distrho6D_collect) {
        if(mpi_rank == mpi_root) {
            MPI_Reduce(MPI_IN_PLACE, offload_array,
                data->offload_dist_length, mpi_type_real, MPI_SUM,
                0, MPI_COMM_WORLD);
        } else {
            MPI_Reduce(offload_array, offload_array,
                data->offload_dist_length, mpi_type_real, MPI_SUM,
                0, MPI_COMM_WORLD);
        }
    }

    if(data->diagorb_collect) {
        if(mpi_rank == mpi_root) {
            for(int i = 1; i < mpi_size; i++) {
                int start_index, n;
                mpi_my_particles(&start_index, &n, ntotal, i, mpi_size);

                for(int j = 0; j < data->diagorb.Nfld; j++) {
                    MPI_Recv(&offload_array[data->offload_diagorb_index
                                        +j*data->diagorb.Nmrk*data->diagorb.Npnt
                                        +start_index*data->diagorb.Npnt],
                        n*data->diagorb.Npnt, mpi_type_real, i, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        else {
            int start_index, n;
            mpi_my_particles(&start_index, &n, ntotal, mpi_rank, mpi_size);

            for(int j = 0; j < data->diagorb.Nfld; j++) {
                MPI_Send(&offload_array[data->offload_diagorb_index
                                      +j*data->diagorb.Nmrk*data->diagorb.Npnt],
                n*data->diagorb.Npnt, mpi_type_real, 0, 0, MPI_COMM_WORLD);
            }
        }
    }

    if(data->diagtrcof_collect) {
        /* 3 fields for transport coefficients, id, D, K */
        int nfield = 3;

        if(mpi_rank == mpi_root) {
            for(int i = 1; i < mpi_size; i++) {
                int start_index, n;
                mpi_my_particles(&start_index, &n, ntotal, i, mpi_size);

                for(int j = 0; j < nfield; j++) {
                    MPI_Recv(&offload_array[data->offload_diagtrcof_index
                                        +j*data->diagtrcof.Nmrk
                                        +start_index],
                        n, mpi_type_real, i, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        else {
            int start_index, n;
            mpi_my_particles(&start_index, &n, ntotal, mpi_rank, mpi_size);

            for(int j = 0; j < nfield; j++) {
                MPI_Send(&offload_array[data->offload_diagtrcof_index
                                      +j*data->diagtrcof.Nmrk],
                n, mpi_type_real, 0, 0, MPI_COMM_WORLD);
            }
        }
    }

#endif
}
