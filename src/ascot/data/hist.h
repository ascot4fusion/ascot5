/**
 * @file hist.h
 * @brief Header file for hist.c
 */
#ifndef HIST_H
#define HIST_H

#include <stdlib.h>
#include "defines.h"
#include "marker.h"
#include "diag.h"



int hist_init(histogram* data, int dimensions, hist_coordinate* coordinates,
              real* binmin, real* binmax, size_t* nbin);
void hist_free(histogram* data);
void hist_offload(histogram* data);
void hist_update_fo(histogram* hist, MarkerGyroOrbit* p_f,
                    MarkerGyroOrbit* p_i);
void hist_update_gc(histogram* hist, MarkerGuidingCenter* p_f,
                    MarkerGuidingCenter* p_i);
#endif
