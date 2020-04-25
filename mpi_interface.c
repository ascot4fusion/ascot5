/**
 * @file mpi_interface.c
 * @brief MPI interface
 *
 * This module provides interfaces for communication with MPI.
 */
#include <mpi.h>
#include <stddef.h>
#include <stdlib.h>
#include "ascot5.h"
#include "mpi_interface.h"

void mpi_interface_init(int argc, char** argv, int* mpi_rank, int* mpi_size) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, mpi_size);
}

void mpi_my_particles(int* start_index, int* n, int ntotal, int mpi_rank,
                      int mpi_size) {
    if(mpi_rank == mpi_size-1) {
        *n = ntotal - mpi_rank * (ntotal / mpi_size);
    }
    else {
        *n = ntotal / mpi_size;
    }

    *start_index = mpi_rank * (ntotal / mpi_size);
}

particle_state* mpi_gather_particlestates(particle_state* ps, int ntotal, int mpi_rank,
                               int mpi_size) {
/*
    MPI_Datatype mpi_type_tmp, mpi_type_particlestate;
    MPI_Aint lb, extent;

    MPI_Type_create_struct(mpi_particlestate_count,
        mpi_particlestate_blocklengths, mpi_particlestate_displacements,
        mpi_particlestate_types, &mpi_type_tmp);
    MPI_Type_get_extent(mpi_type_tmp, &lb, &extent);
    MPI_Type_create_resized(mpi_type_tmp, lb, extent, &mpi_type_particlestate);
    MPI_Type_commit(&mpi_type_particlestate);

    int start_index, n;
    mpi_my_particles(&start_index, &n, mpi_rank, mpi_size);

    MPI_Gather(ps+start_index, n, mpi_type_particlestate, ps, n,
        mpi_type_particlestate, 0, MPI_COMM_WORLD);
*/
    const int n_real = 31;
    const int n_int = 3;
    const int n_err = 1;

    if(mpi_rank == 0) {
        int start_index, n;
        mpi_my_particles(&start_index, &n, ntotal, mpi_rank, mpi_size);

        particle_state* ps_all = (particle_state*) malloc(ntotal*sizeof(particle_state));

        for(int j = 0; j < n; j++) {
            ps_all[j] = ps[j];
        }

        for(int i = 1; i < mpi_size; i++) {
            mpi_my_particles(&start_index, &n, ntotal, i, mpi_size);
            printf("%d %d %d %d %d\n", start_index, n, ntotal, i, mpi_size);

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
                ps_all[start_index+j].r      = realdata[0*n+j];
                ps_all[start_index+j].phi    = realdata[1*n+j];
                ps_all[start_index+j].z      = realdata[2*n+j];
                ps_all[start_index+j].vpar   = realdata[3*n+j];
                ps_all[start_index+j].mu     = realdata[4*n+j];
                ps_all[start_index+j].zeta   = realdata[5*n+j];
                ps_all[start_index+j].rprt   = realdata[6*n+j];
                ps_all[start_index+j].phiprt = realdata[7*n+j];
                ps_all[start_index+j].zprt   = realdata[8*n+j];
                ps_all[start_index+j].rdot   = realdata[9*n+j];
                ps_all[start_index+j].phidot = realdata[10*n+j];
                ps_all[start_index+j].zdot   = realdata[11*n+j];
                ps_all[start_index+j].mass   = realdata[12*n+j];
                ps_all[start_index+j].charge = realdata[13*n+j];
                ps_all[start_index+j].weight = realdata[14*n+j];
                ps_all[start_index+j].time   = realdata[15*n+j];
                ps_all[start_index+j].cputime = realdata[16*n+j];
                ps_all[start_index+j].rho    = realdata[17*n+j];
                ps_all[start_index+j].theta  = realdata[18*n+j];
                ps_all[start_index+j].id       = intdata[0*n+j];
                ps_all[start_index+j].endcond  = intdata[1*n+j];
                ps_all[start_index+j].walltile = intdata[2*n+j];
                ps_all[start_index+j].B_r    = realdata[19*n+j];
                ps_all[start_index+j].B_phi  = realdata[20*n+j];
                ps_all[start_index+j].B_z    = realdata[21*n+j];
                ps_all[start_index+j].B_r_dr = realdata[22*n+j];
                ps_all[start_index+j].B_phi_dr = realdata[23*n+j];
                ps_all[start_index+j].B_z_dr   = realdata[24*n+j];
                ps_all[start_index+j].B_r_dphi = realdata[25*n+j];
                ps_all[start_index+j].B_phi_dphi = realdata[26*n+j];
                ps_all[start_index+j].B_z_dphi = realdata[27*n+j];
                ps_all[start_index+j].B_r_dz   = realdata[28*n+j];
                ps_all[start_index+j].B_phi_dz = realdata[29*n+j];
                ps_all[start_index+j].B_z_dz   = realdata[30*n+j];
                ps_all[start_index+j].err      = errdata[j];
            }

            free(realdata);
            free(intdata);
            free(errdata);
        }

        return ps_all;
    }
    else {

        int start_index, n;
        mpi_my_particles(&start_index, &n, ntotal, mpi_rank, mpi_size);

        real* realdata;
        realdata = malloc(n_real * n * sizeof(realdata));
        integer* intdata;
        intdata = malloc(n_int * n * sizeof(intdata));
        a5err* errdata;
        errdata = malloc(n_err * n * sizeof(intdata));

        for(int j = 0; j < n; j++) {
            realdata[0*n+j] = ps[j].r;
            realdata[1*n+j] = ps[j].phi;
            realdata[2*n+j] = ps[j].z;
            realdata[3*n+j] = ps[j].vpar;
            realdata[4*n+j] = ps[j].mu;
            realdata[5*n+j] = ps[j].zeta;
            realdata[6*n+j] = ps[j].rprt;
            realdata[7*n+j] = ps[j].phiprt;
            realdata[8*n+j] = ps[j].zprt;
            realdata[9*n+j] = ps[j].rdot;
            realdata[10*n+j] = ps[j].phidot;
            realdata[11*n+j] = ps[j].zdot;
            realdata[12*n+j] = ps[j].mass;
            realdata[13*n+j] = ps[j].charge;
            realdata[14*n+j] = ps[j].weight;
            realdata[15*n+j] = ps[j].time;
            realdata[16*n+j] = ps[j].cputime;
            realdata[17*n+j] = ps[j].rho;
            realdata[18*n+j] = ps[j].theta;
            intdata[0*n+j] = ps[j].id;
            intdata[1*n+j] = ps[j].endcond;
            intdata[2*n+j] = ps[j].walltile;
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
            errdata[j] = ps[j].err;
        }

        MPI_Send(realdata, n_real*n, mpi_type_real, 0, 0, MPI_COMM_WORLD);
        MPI_Send(intdata, n_int*n, mpi_type_integer, 0, 0, MPI_COMM_WORLD);
        MPI_Send(errdata, n_err*n, mpi_type_a5err, 0, 0, MPI_COMM_WORLD);

        free(realdata);
        free(intdata);
        free(errdata);

        return (particle_state*) malloc(1);
    }

}
