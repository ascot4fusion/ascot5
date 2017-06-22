/**
 * @file ascot4_interface.h
 * @brief Header file for ascot4_interface.c
 */
#ifndef ASCOT4_INTERFACE_H
#define ASCOT4_INTERFACE_H
#include "ascot5.h"
#include "B_field.h"
#include "particle.h"
#include "distributions.h"

void ascot4_write_B(B_field_data* Bdata);
void ascot4_write_particles(particle* p, int n, char* filename,
                            char* comment);
void ascot4_read_particles(input_particle** p, int* n, char* filename);
void ascot4_write_dist_rzvv(dist_rzvv_offload_data* dist, real* hist,
                            char* filename);
void ascot4_write_inistate(int n, input_particle* p, char* filename);
void ascot4_write_endstate(int n, input_particle* p, char* filename);

#endif
