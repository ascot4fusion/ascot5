/**
 * @file mpi_interface.h
 * @brief Header file for mpi_interface.c
 */
#ifndef MPI_INTERFACE_H
#define MPI_INTERFACE_H

void mpi_interface_init(int argc, char** argv, int* mpi_rank, int* mpi_size);
void mpi_my_particles(int* start_index, int* n, int mpi_rank, int mpi_size);

#endif
