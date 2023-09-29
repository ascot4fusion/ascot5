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

void mpi_interface_init(int argc, char** argv, sim_offload_data* sim,
                        int* mpi_rank, int* mpi_size, int* mpi_root);
void mpi_interface_finalize();

void mpi_my_particles(int* start_index, int* n, int ntotal, int mpi_rank,
                      int mpi_size);

void mpi_gather_particlestate(particle_state* ps, particle_state** psgathered,
                              int* ngathered, int ntotal, int mpi_rank,
                              int mpi_size, int mpi_root);
void mpi_gather_diag(diag_offload_data* data, real* offload_array, int ntotal,
                     int mpi_rank, int mpi_size, int mpi_root);

#endif
