/**
 * @file mpi_interface.h
 * @brief Header file for mpi_interface.c
 */
#ifndef MPI_INTERFACE_H
#define MPI_INTERFACE_H

#include <mpi.h>
#include "particle.h"

/*
#if defined SINGLEPRECISION
typedef MPI_INT mpi_type_integer;
typedef MPI_FLOAT mpi_type_real;
#else
typedef MPI_LONG mpi_type_integer;
typedef MPI_DOUBLE mpi_type_real;
#endif

typedef MPI_UNSIGNED_LONG_LONG mpi_type_a5err;
*/

#define mpi_type_integer MPI_LONG
#define mpi_type_real MPI_DOUBLE
#define mpi_type_a5err MPI_UNSIGNED_LONG_LONG
/*
const int mpi_particlestate_count = 35;
const int mpi_particlestate_blocklengths[35] = {1};
const MPI_Aint mpi_particlestate_displacements[35]={offsetof(particle_state,r),
                                        offsetof(particle_state,phi),
                                        offsetof(particle_state,z),
                                        offsetof(particle_state,vpar),
                                        offsetof(particle_state,mu),
                                        offsetof(particle_state,theta),
                                        offsetof(particle_state,rprt),
                                        offsetof(particle_state,phiprt),
                                        offsetof(particle_state,zprt),
                                        offsetof(particle_state,rdot),
                                        offsetof(particle_state,phidot),
                                        offsetof(particle_state,zdot),
                                        offsetof(particle_state,mass),
                                        offsetof(particle_state,charge),
                                        offsetof(particle_state,weight),
                                        offsetof(particle_state,time),
                                        offsetof(particle_state,cputime),
                                        offsetof(particle_state,rho),
                                        offsetof(particle_state,pol),
                                        offsetof(particle_state,id),
                                        offsetof(particle_state,endcond),
                                        offsetof(particle_state,walltile),
                                        offsetof(particle_state,B_r),
                                        offsetof(particle_state,B_phi),
                                        offsetof(particle_state,B_z),
                                        offsetof(particle_state,B_r_dr),
                                        offsetof(particle_state,B_phi_dr),
                                        offsetof(particle_state,B_z_dr),
                                        offsetof(particle_state,B_r_dphi),
                                        offsetof(particle_state,B_phi_dphi),
                                        offsetof(particle_state,B_z_dphi),
                                        offsetof(particle_state,B_r_dz),
                                        offsetof(particle_state,B_phi_dz),
                                        offsetof(particle_state,B_z_dz),
                                        offsetof(particle_state,err)};
const MPI_Datatype mpi_particlestate_types[35] = {mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_integer,
                                         mpi_type_integer,
                                         mpi_type_integer,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_real,
                                         mpi_type_a5err};
*/
void mpi_interface_init(int argc, char** argv, int* mpi_rank, int* mpi_size);
void mpi_my_particles(int* start_index, int* n, int ntotal, int mpi_rank, int mpi_size);
particle_state* mpi_gather_particlestates(particle_state* ps, int ntotal, int mpi_rank, int mpi_size);

#endif
