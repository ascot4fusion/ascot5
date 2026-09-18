/**
 * @file orbit_following.h
 * Orbit integrators.
 *
 * These functions solve the deterministic motion of the marker due to the
 * electric and magnetic fields. Also ballistic motion is solved when
 * applicable.
 */
#ifndef ORBIT_FOLLOWING_H
#define ORBIT_FOLLOWING_H

#include "data/bfield.h"
#include "data/boozer.h"
#include "data/marker.h"
#include "data/mhd.h"

/**
 * Trace markers representing magnetic field lines for a single step.
 *
 * This function calculates a magnetic field line step simultaneously for a
 * vector of markers with the Cash-Karp (adaptive RK5) method. All arrays in the
 * function are of same length so vectorization can be performed directly.
 * Informs whether time step was accepted or rejected and provides a suggestion
 * for the next time step.
 *
 * @param mrk Markers that are advanced.
 * @param h Integration step for each marker [m].
 * @param hnext Suggestion for the next step size [m].
 *        Negative sign indicates a failed step. In this case, the suggested
 *        value for the next step is the absolute value.
 * @param tol Error tolerance for acceptance.
 * @param bfield Magnetic field data.
 */
void step_fl_cashkarp(
    MarkerFieldLine *mrk, const real *h, real *hnext, real tol, Bfield *bfield);

/**
 * Trace markers representing magnetic field lines for a single step when MHD
 * perturbations are present.
 *
 * This function is identical to step_fl_cashkarp but with MHD perturbations.
 *
 * @param mrk Markers that are advanced.
 * @param h Time step for each marker [m].
 * @param hnext Suggestion for the next step size [m].
 *        Negative sign indicates a failed step. In this case, the suggested
 *        value for the next step is the absolute value.
 * @param tol Error tolerance for acceptance.
 * @param bfield Magnetic field data.
 * @param boozer Boozer data.
 * @param mhd MHD data.
 */
void step_fl_cashkarp_mhd(
    MarkerFieldLine *mrk, const real *h, real *hnext, real tol, Bfield *bfield,
    Boozer *boozer, Mhd *mhd);

/**
 * Integrate full orbit step with VPA.
 *
 * The integration is performed for a vector of markers simultaneously using the
 * volume preserving algorithm (Boris method for relativistic particles) see
 * Zhang 2015 https://doi.org/10.1063/1.4916570. All arrays in the function are
 * of same length so vectorization can be performed directly.
 *
 * This algorithm is valid for neutral particles as well, in which case the
 * motion reduces to ballistic motion where momentum remains constant.
 *
 * @param mrk Markers that are advanced.
 * @param h Time step for each marker [s].
 * @param bfield Magnetic field data.
 * @param efield Electric field data.
 * @param aldforce Toggle for Abraham-Lorentz-Dirac force.
 */
void step_go_vpa(
    MarkerGyroOrbit *mrk, const real *h, Bfield *bfield, Efield *efield,
    int aldforce);

/**
 * Integrate full orbit step with VPA when MHD perturbation is present.
 *
 * Identical to step_go_vpa but with MHD present.
 *
 * @param mrk Markers that are advanced.
 * @param h Time step for each marker [s].
 * @param bfield Magnetic field data.
 * @param efield Electric field data.
 * @param boozer Boozer data.
 * @param mhd MHD data.
 * @param aldforce Toggle for Abraham-Lorentz-Dirac force.
 */
void step_go_vpa_mhd(
    MarkerGyroOrbit *mrk, const real *h, Bfield *bfield, Efield *efield,
    Boozer *boozer, Mhd *mhd, int aldforce);

/**
 * Integrate guiding center step with fixed time-step.
 *
 * The integration is performed for a vector of markers simultaneously with RK4.
 * All arrays in the function are of same length so vectorization can be
 * performed directly.
 *
 * @param mrk Markers that are advanced.
 * @param h Time step for each marker [s].
 * @param bfield Magnetic field data.
 * @param efield Electric field data.
 * @param aldforce Toggle for Abraham-Lorentz-Dirac force.
 */
void step_gc_rk4(
    MarkerGuidingCenter *mrk, const real *h, Bfield *bfield, Efield *efield,
    int aldforce);

/**
 * Integrate guiding center step with fixed time-step and with MHD present.
 *
 * Identical to step_gc_rk4 but with MHD present.
 *
 * @param mrk Markers that are advanced.
 * @param h Time step for each marker [s].
 * @param bfield Magnetic field data.
 * @param efield Electric field data.
 * @param boozer Boozer data.
 * @param mhd MHD data.
 * @param aldforce Toggle for Abraham-Lorentz-Dirac force.
 */
void step_gc_rk4_mhd(
    MarkerGuidingCenter *mrk, const real *h, Bfield *bfield, Efield *efield,
    Boozer *boozer, Mhd *mhd, int aldforce);

/**
 * Integrate guiding center step with adaptive time-step.
 *
 * The integration is performed for a vector of markers simultaneously with
 * Cash-Karp method. All arrays in the function are of same length so
 * vectorization can be performed directly. Informs whether time step was
 * accepted or rejected and provides a suggestion for the next time step.
 *
 * @param mrk Markers that are advanced.
 * @param h Time step for each marker [s].
 * @param hnext Suggestion for the next step size [s].
 *        Negative sign indicates a failed step. In this case, the suggested
 *        value for the next step is the absolute value.
 * @param tol Error tolerance for acceptance.
 * @param bfield Magnetic field data.
 * @param efield Electric field data.
 * @param aldforce Toggle for Abraham-Lorentz-Dirac force.
 */
void step_gc_cashkarp(
    MarkerGuidingCenter *mrk, const real *h, real *hnext, real tol,
    Bfield *bfield, Efield *efield, int aldforce);

/**
 * Integrate guiding center step with adaptive time-step and with MHD present.
 *
 * Identical to step_gc_cashkarp but with MHD present.
 *
 * @param mrk Markers that are advanced.
 * @param h Time step for each marker [s].
 * @param hnext Suggestion for the next step size [s].
 *        Negative sign indicates a failed step. In this case, the suggested
 *        value for the next step is the absolute value.
 * @param tol Error tolerance for acceptance.
 * @param bfield Magnetic field data.
 * @param efield Electric field data.
 * @param boozer Boozer data.
 * @param mhd MHD data.
 * @param aldforce Toggle for Abraham-Lorentz-Dirac force.
 */
void step_gc_cashkarp_mhd(
    MarkerGuidingCenter *mrk, const real *h, real *hnext, real tol,
    Bfield *bfield, Efield *efield, Boozer *boozer, Mhd *mhd, int aldforce);

#endif
