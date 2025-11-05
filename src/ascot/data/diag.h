/**
 * @file diag.h
 * Interface for simulation diagnostics.
 *
 * Any output gathered during the simulation loop is collected here.
 */
#ifndef DIAG_H
#define DIAG_H
#include "bfield.h"
#include "defines.h"
#include "diag_data.h"
#include "marker.h"
#include "options.h"

/**
 * Offload data to the accelerator.
 *
 * @param diag The struct to offload.
 */
void Diag_offload(Diagnostics *diag);

/**
 * Collect diagnostics for gyro-orbit markers.
 *
 * @param diag Diagnostics data.
 * @param bfield Magnetic field data.
 * @param mrk_f Marker at the end of the time-step.
 * @param mrk_i Marker at the beginning of the time-step.
 */
void Diag_update_go(
    Diagnostics *diag, Bfield *bfield, MarkerGyroOrbit *mrk_f,
    MarkerGyroOrbit *mrk_i);

/**
 * Collect diagnostics for guiding-center markers.
 *
 * @param diag Diagnostics data.
 * @param bfield Magnetic field data.
 * @param mrk_f Marker at the end of the time-step.
 * @param mrk_i Marker at the beginning of the time-step.
 */
void Diag_update_gc(
    Diagnostics *diag, Bfield *bfield, MarkerGuidingCenter *mrk_f,
    MarkerGuidingCenter *mrk_i);

/**
 * Collect diagnostics for field-line markers.
 *
 * Histograms are not collected for magnetic field lines.
 *
 * @param diag Diagnostics data.
 * @param bfield Magnetic field data.
 * @param mrk_f Marker at the end of the time-step.
 * @param mrk_i Marker at the beginning of the time-step.
 */
void Diag_update_fl(
    Diagnostics *diag, Bfield *bfield, MarkerFieldLine *mrk_f,
    MarkerFieldLine *mrk_i);

#endif
