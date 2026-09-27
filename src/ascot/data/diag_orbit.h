/**
 * @file diag_orbit.h
 * Orbit diagnostic.
 */
#ifndef DIAG_ORB_H
#define DIAG_ORB_H

#include "diag.h"
#include "marker.h"
#include "bfield.h"
#include <stdio.h>

#define DIAG_ORB_GO 1 /**< Data stored in FO mode */
#define DIAG_ORB_GC 2 /**< Data stored in GC mode */
#define DIAG_ORB_FL 3 /**< Data stored in ML mode */


/**
 *
 * @param orbit Orbit diagnostics.
 */
void DiagOrbit_offload(DiagOrbit *orbit);

/**
 *
 * @param orbit Orbit diagnostics.
 */
void DiagOrbit_onload(DiagOrbit *orbit);

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
