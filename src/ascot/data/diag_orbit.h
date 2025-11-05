/**
 * @file diag_orb.h
 * @brief Header file for diag_orb.c.
 *
 * This file also contains definitions for orbit diagnostics data structures.
 */
#ifndef DIAG_ORB_H
#define DIAG_ORB_H

#include "diag.h"
#include "marker.h"
#include "bfield.h"
#include <stdio.h>

#define DIAG_ORB_FO 1 /**< Data stored in FO mode */
#define DIAG_ORB_GC 2 /**< Data stored in GC mode */
#define DIAG_ORB_ML 3 /**< Data stored in ML mode */

DECLARE_TARGET_SIMD_UNIFORM(ang0)
/**
 * @brief Check if marker has crossed a plane.
 *
 * This helper function checks whether the angle, either toroidal or poloidal,
 * that defines a Poincare plane is between marker's initial and final angles
 * (of single timestep).
 *
 * @param fang marker final angle in radians.
 * @param iang marker initial angle in radians.
 * @param ang0 Poincare plane angle.
 *
 * @return zero if no-crossing, number k, ang0 = k + (fang - iang), otherwise.
 */
real diag_orb_check_plane_crossing(real fang, real iang, real ang0);

DECLARE_TARGET_SIMD_UNIFORM(r0)
/**
 * @brief Check if marker has crossed given rho
 *
 * This helper function checks whether given rho that defines a Poincare plane
 * is between marker's initial and final rhos (of single timestep).
 *
 * @param frho marker final rho in metres.
 * @param irho marker initial rho in metres.
 * @param rho0 Poincare plane rho.
 *
 * @return zero if no-crossing, number k, rho0 = k * (frho - irho), otherwise.
 */
real diag_orb_check_radial_crossing(real fr, real ir, real r0);

void DiagOrbit_offload(DiagOrbit *orbit);

/**
 * Record orbit for gyro-orbit markers.
 *
 * @param orbit Orbit diagnostics.
 * @param bfield Magnetic field data.
 * @param mrk_f Marker at the end of the time-step.
 * @param mrk_i Marker at the beginning of the time-step.
 */
void DiagOrbit_update_go(
    DiagOrbit *orbit, Bfield *bfield, MarkerGyroOrbit *mrk_f,
    MarkerGyroOrbit *mrk_i);

/**
 * Record orbit for guiding-center markers.
 *
 * @param orbit Orbit diagnostics.
 * @param bfield Magnetic field data.
 * @param mrk_f Marker at the end of the time-step.
 * @param mrk_i Marker at the beginning of the time-step.
 */
void DiagOrbit_update_gc(
    DiagOrbit *orbit, Bfield *bfield, MarkerGuidingCenter *mrk_f,
    MarkerGuidingCenter *mrk_i);

/**
 * Record orbit for field-line markers.
 *
 * @param orbit Orbit diagnostics.
 * @param bfield Magnetic field data.
 * @param mrk_f Marker at the end of the time-step.
 * @param mrk_i Marker at the beginning of the time-step.
 */
void DiagOrbit_update_fl(
    DiagOrbit *orbit, Bfield *bfield, MarkerFieldLine *mrk_f,
    MarkerFieldLine *mrk_i);

#endif
