/**
 * @file mpi_interface.h
 * @brief Header file for mpi_interface.c
 */
#ifndef MPI_INTERFACE_H
#define MPI_INTERFACE_H

#ifdef MPI
#include <mpi.h>
#endif
#include "diag.h"
#include "particle.h"
#include "simulate.h"

#define mpi_type_integer MPI_LONG
#define mpi_type_real MPI_DOUBLE
#define mpi_type_a5err MPI_UNSIGNED_LONG_LONG

void mpi_interface_init(int argc, char** argv, sim_offload_data* sim,
                        int* mpi_rank, int* mpi_size);
void mpi_interface_finalize();

void mpi_my_particles(int* start_index, int* n, int ntotal, int mpi_rank,
                      int mpi_size);

void mpi_gather_particlestate(particle_state* ps, particle_state* ps_all,
                              int ntotal, int mpi_rank, int mpi_size);
void mpi_gather_diag(diag_offload_data* data, real* offload_array, int ntotal,
                     int mpi_rank, int mpi_size);

#endif
