/**
 * @file orbit_write.h
 * @brief Functions to write particle and guiding center information.
 */
#include <stdio.h>
#include "particle.h"
#include "B_field.h"

void write_particle(FILE* out, particle* p);


void write_guidingcenter(FILE* out, particle* p);

void write_fo_as_particle(FILE* out, particle_simd_fo* p);

void write_fo_as_guidingcenter(FILE* out, particle_simd_fo* p, B_field_data* Bdata);

void write_gc_as_particle(FILE* out, particle* p);

void write_gc_as_guidingcenter(FILE* out, particle_simd_gc* p);
