/**
 * @file diag_moment.h
 * Diagnostic that collects particle moments.
 *
 * The moments are collected into a struct indexed by rho where a single rho
 * can correspond to multiple moments.
 */
#ifndef DIAG_MOMENT_H
#define DIAG_MOMENT_H

#include "defines.h"
#include "bfield_data.h"
#include "diag.h"
#include "marker.h"
#include <stdlib.h>

/**
 * Offload data to the accelerator.
 *
 * @param moment The struct to offload.
 */
void DiagMoment_offload(DiagMoment *moment);

/**
 * Onload data to the accelerator.
 *
 * @param moment The struct to offload.
 */
void DiagMoment_onload(DiagMoment *moment);

/**
 * Update moments for gyro-orbit markers.
 *
 * @param moment The moment data.
 * @param bfield Magnetic field data.
 * @param mrk_f Marker at the end of the time-step.
 * @param mrk_i Marker at the beginning of the time-step.
 */
void DiagMoment_update_go(
    DiagMoment *moment, Bfield *bfield, MarkerGyroOrbit *mrk_f, MarkerGyroOrbit *mrk_i);

/**
 * Update moment for guiding-center markers.
 *
 * @param moment The moment data.
 * @param bfield Magnetic field data.
 * @param mrk_f Marker at the end of the time-step.
 * @param mrk_i Marker at the beginning of the time-step.
 */
void DiagMoment_update_gc(
    DiagMoment *moment, Bfield *bfield, MarkerGuidingCenter *mrk_f, MarkerGuidingCenter *mrk_i);
#endif
