/**
 * @file orbit_following.h
 */
#ifndef ORBIT_FOLLOWING_H
#define ORBIT_FOLLOWING_H

#include "bfield.h"
#include "boozer.h"
#include "mhd.h"
#include "particle.h"

void step_ml_cashkarp(particle_simd_ml* p, real* h, real* hnext,
                      real tol, Bfield* bfield);
void step_ml_cashkarp_mhd(particle_simd_ml* p, real* h, real* hnext,
                          real tol, Bfield* bfield,
                          Boozer* boozer,
                          Mhd* mhd);

void step_fo_vpa(particle_simd_fo* p, real* h, Bfield* bfield,
                 Efield* efield, int aldforce);
void step_fo_vpa_mhd(
    particle_simd_fo* p, real* h, Bfield* bfield, Efield* efield,
    Boozer* boozer, Mhd* mhd, int aldforce);

void step_gc_cashkarp(particle_simd_gc* p, real* h, real* hnext, real tol,
                      Bfield* bfield, Efield* efield, int aldforce);
void step_gc_cashkarp_mhd(
    particle_simd_gc* p, real* h, real* hnext, real tol, Bfield* bfield,
    Efield* efield, Boozer* boozer, Mhd* mhd, int aldforce);

void step_gc_rk4(particle_simd_gc* p, real* h, Bfield* bfield,
                 Efield* efield, int aldforce);
void step_gc_rk4_mhd(particle_simd_gc* p, real* h, Bfield* bfield,
                     Efield* efield, Boozer* boozer,
                     Mhd* mhd, int aldforce);

#endif
