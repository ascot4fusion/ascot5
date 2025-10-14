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

/**
 * @brief Integrate a magnetic field line step for a struct of markers
 *
 * This function calculates a magnetic field line step for a struct of NSIMD
 * markers with Cash-Karp (adaptive RK5) simultaneously using SIMD instructions.
 * All arrays in the function are of NSIMD length so vectorization can be
 * performed directly without gather and scatter operations. Informs whether
 * time step was accepted or rejected and provides a suggestion for the next
 * time step.
 *
 * @param p marker struct that will be integrated
 * @param h NSIMD length array containing time step lengths
 * @param hnext suggestion for the next time step. Negative if rejected.
 * @param tol error tolerance
 * @param bfield pointer to magnetic field data
 * @param boozerdata pointer to Boozer data
 * @param mhddata pointer to MHD data
 */
void step_fl_cashkarp_mhd(
    MarkerFieldLine *p, real *h, real *hnext, real tol, Bfield *bfield,
    Boozer *boozer, Mhd *mhd);

/**
 * @brief Integrate a full orbit step for a struct of particles with VPA
 *
 * The integration is performed for a struct of NSIMD particles using the
 * volume preserving algorithm (Boris method for relativistic particles) see
 * Zhang 2015.
 *
 * This algorithm is valid for neutral particles as well, in which case the
 * motion reduces to ballistic motion where momentum remains constant.
 *
 * @param p MarkerGyroOrbit struct that will be updated
 * @param h pointer to array containing time steps
 * @param bfield pointer to magnetic field data
 * @param efield pointer to electric field data
 * @param aldforce indicates whether Abraham-Lorentz-Dirac force is enabled
 */
void step_go_vpa(
    MarkerGyroOrbit *p, real *h, Bfield *bfield, Efield *efield, int aldforce);

/**
 * @brief Integrate a full orbit step with VPA and MHd modes present.
 *
 * Same as previous method but with MHD present.
 *
 * @param p MarkerGyroOrbit struct that will be updated
 * @param h pointer to array containing time steps
 * @param bfield pointer to magnetic field data
 * @param efield pointer to electric field data
 * @param boozer pointer to boozer data
 * @param mhd pointer to MHD data
 * @param aldforce indicates whether Abraham-Lorentz-Dirac force is enabled
 */
void step_go_vpa_mhd(
    MarkerGyroOrbit *p, real *h, Bfield *bfield, Efield *efield, Boozer *boozer,
    Mhd *mhd, int aldforce);

/**
 * @brief Integrate a guiding center step for a struct of markers
 *
 * This function calculates a guiding center step for a struct of NSIMD
 * markers with Cash-Karp (adaptive RK5) simultaneously using SIMD instructions.
 * All arrays in the function are of NSIMD length so vectorization can be
 * performed directly without gather and scatter operations. Informs whther time
 * step was accepted or rejected and provides a suggestion for the next time
 * step.
 *
 * @param p marker struct that will be updated
 * @param h array containing time step lengths
 * @param hnext suggestion for the next time step. Negative sign indicates
 * current step was rejected
 * @param tol error tolerance
 * @param bfield pointer to magnetic field data
 * @param efield pointer to electric field data
 * @param aldforce indicates whether Abraham-Lorentz-Dirac force is enabled
 */
void step_gc_cashkarp(
    MarkerGuidingCenter *p, real *h, real *hnext, real tol, Bfield *bfield,
    Efield *efield, int aldforce);

/**
 * @brief Integrate a guiding center step for a struct of markers with MHD
 *
 * Rejected step has a negative suggestion for the next time-step. The negative
 * sign is only used to indicate a rejected step and absolute value should be
 * used for the next time-step.
 *
 * @param p marker struct that will be updated
 * @param h array containing time step lengths
 * @param hnext suggestion for the next time step.
 * @param tol error tolerance
 * @param bfield pointer to magnetic field data
 * @param efield pointer to electric field data
 * @param boozer pointer to Boozer data
 * @param mhd pointer to MHD data
 * @param aldforce indicates whether Abraham-Lorentz-Dirac force is enabled
 */
void step_gc_cashkarp_mhd(
    MarkerGuidingCenter *p, real *h, real *hnext, real tol, Bfield *bfield,
    Efield *efield, Boozer *boozer, Mhd *mhd, int aldforce);

/**
 * @brief Integrate a guiding center step for a struct of markers with RK4
 *
 * This function calculates a guiding center step for a struct of NSIMD
 * markers with RK4 simultaneously using SIMD instructions. All arrays in the
 * function are of NSIMD length so vectorization can be performed directly
 * without gather and scatter operations.
 *
 * @param p simd_gc struct that will be updated
 * @param h pointer to array containing time steps
 * @param bfield pointer to magnetic field data
 * @param efield pointer to electric field data
 * @param aldforce indicates whether Abraham-Lorentz-Dirac force is enabled
 */
void step_gc_rk4(
    MarkerGuidingCenter *p, real *h, Bfield *bfield, Efield *efield,
    int aldforce);

/**
 * @brief Integrate a guiding center step with RK4 with MHD modes present.
 *
 * Same as previous function but with MHD present
 *
 * @param p simd_gc struct that will be updated
 * @param h pointer to array containing time steps
 * @param bfield pointer to magnetic field data
 * @param efield pointer to electric field data
 * @param boozer pointer to boozer data
 * @param mhd pointer to MHD data
 * @param aldforce indicates whether Abraham-Lorentz-Dirac force is enabled
 */
void step_gc_rk4_mhd(
    MarkerGuidingCenter *p, real *h, Bfield *bfield, Efield *efield,
    Boozer *boozer, Mhd *mhd, int aldforce);

#endif
