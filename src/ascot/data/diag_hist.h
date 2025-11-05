/**
 * @file hist.h
 * Diagnostic that collects particle histograms.
 *
 * A histogram can have an arbitrary number of dimensions consisting of any
 * of the predefined coordinate axes.
 */
#ifndef DIAG_HIST_H
#define DIAG_HIST_H

#include "defines.h"
#include "diag.h"
#include "marker.h"
#include <stdlib.h>

/**
 * Offload data to the accelerator.
 *
 * @param hist The struct to offload.
 */
void DiagHist_offload(DiagHist *hist);

/**
 * Update histogram for gyro-orbit markers.
 *
 * @param hist The histogram data.
 * @param mrk_f Marker at the end of the time-step.
 * @param mrk_i Marker at the beginning of the time-step.
 *
 */
void DiagHist_update_go(
    DiagHist *hist, MarkerGyroOrbit *mrk_f, MarkerGyroOrbit *mrk_i);

void DiagHist_update_gc(
    DiagHist *hist, MarkerGuidingCenter *mrk_f, MarkerGuidingCenter *mrk_i);
#endif
