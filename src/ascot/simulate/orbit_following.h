/**
 * @file orbit_following.h
 */
#ifndef ORBIT_FOLLOWING_H
#define ORBIT_FOLLOWING_H

#include "data/bfield.h"
#include "data/boozer.h"
#include "data/mhd.h"
#include "data/marker.h"

/**
 * Trace markers representing magnetic field lines for a single step.
 *
 * This function calculates a magnetic field line step for a struct of NSIMD
 * markers with Cash-Karp (adaptive RK5) simultaneously using SIMD instructions.
 * All arrays in the function are of NSIMD length so vectorization can be
 * performed directly without gather and scatter operations. Informs whether
 * time step was accepted or rejected and provides a suggestion for the next
 * time step.
 *
 * @param p Markers that are advanced.
 * @param h NSIMD length array containing time step lengths
 * @param hnext Suggestion for the next step size.
 *        Negative sign indicates a failed step and the suggested value for the
 *        next step is the absoulte value.
 * @param tol Error tolerance for acceptance.
 * @param bfield Magnetic field data.
 */
void step_fl_cashkarp(
    MarkerFieldLine *p, real *h, real *hnext, real tol, Bfield *bfield);

void step_fl_cashkarp_mhd(
    MarkerFieldLine *p, real *h, real *hnext, real tol, Bfield *bfield,
    Boozer *boozer, Mhd *mhd);

void step_go_vpa(
    MarkerGyroOrbit *p, real *h, Bfield *bfield, Efield *efield, int aldforce);
void step_go_vpa_mhd(
    MarkerGyroOrbit *p, real *h, Bfield *bfield, Efield *efield,
    Boozer *boozer, Mhd *mhd, int aldforce);

void step_gc_cashkarp(
    MarkerGuidingCenter *p, real *h, real *hnext, real tol, Bfield *bfield,
    Efield *efield, int aldforce);
void step_gc_cashkarp_mhd(
    MarkerGuidingCenter *p, real *h, real *hnext, real tol, Bfield *bfield,
    Efield *efield, Boozer *boozer, Mhd *mhd, int aldforce);

void step_gc_rk4(
    MarkerGuidingCenter *p, real *h, Bfield *bfield, Efield *efield, int aldforce);
void step_gc_rk4_mhd(
    MarkerGuidingCenter *p, real *h, Bfield *bfield, Efield *efield,
    Boozer *boozer, Mhd *mhd, int aldforce);

#endif
