/**
 * @file mpi_interface.c
 * @brief MPI interface
 *
 * This module provides interfaces for communication with MPI.
 */
#include <mpi.h>

void mpi_interface_init(int argc, char** argv, int* mpi_rank, int* mpi_size) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, mpi_size);
}

void mpi_my_particles(int* start_index, int* n, int mpi_rank, int mpi_size) {
    *start_index = mpi_rank * (*n / mpi_size);

    if(mpi_rank == mpi_size-1) {
        *n = *n - mpi_rank * (*n / mpi_size);
    }
    else {
        *n = *n / mpi_size;
    }
}
