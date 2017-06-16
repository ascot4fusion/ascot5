/**
 * @file orbit_write.h
 * @brief Functions to write particle and guiding center information.
 */
#ifndef ORBIT_WRITE_H
#define ORBIT_WRITE_H

#include <stdio.h>
#include "particle.h"
#include "B_field.h"

#define WRITEORB_SLOTS 100000

typedef struct{ 

    int id[WRITEORB_SLOTS];
    real time[WRITEORB_SLOTS];
    real r[WRITEORB_SLOTS];
    real phi[WRITEORB_SLOTS];
    real z[WRITEORB_SLOTS];
    real vpar[WRITEORB_SLOTS];
    real mu[WRITEORB_SLOTS];
} writeorbit_guidingcenter;

void write_particle(FILE* out, particle* p);

void write_guidingcenter(FILE* out, particle* p);

void write_fo_as_particle(FILE* out, particle_simd_fo* p);

void write_fo_as_guidingcenter(FILE* out, particle_simd_fo* p, B_field_data* Bdata);

void write_gc_as_particle(FILE* out, particle* p);

void write_gc_as_guidingcenter(FILE* out, particle_simd_gc* p);

void writeorbit_store_gc2guidingcenter(particle_simd_gc p, writeorbit_guidingcenter* w, int* write, int Nwrite, int slot);

#endif
