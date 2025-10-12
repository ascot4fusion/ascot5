/**
 * @file diag_transcoef.h
 * @brief Header file for diag_transcoef.c.
 *
 * Contains definitions for transport coefficient data structures.
 */
#ifndef DIAG_TRANSCOEF_H
#define DIAG_TRANSCOEF_H

#include "defines.h"
#include "marker.h"
#include "options.h"
#include "diag.h"


void diag_transcoef_init(diag_transcoef_data* data, Options* options,
                         size_t nmarkers);
void diag_transcoef_free(diag_transcoef_data* data);
void diag_transcoef_update_fo(diag_transcoef_data* data, Options* options,
                              MarkerGyroOrbit* p_f, MarkerGyroOrbit* p_i);
void diag_transcoef_update_gc(diag_transcoef_data* data, Options* options,
                              MarkerGuidingCenter* p_f, MarkerGuidingCenter* p_i);
void diag_transcoef_update_ml(diag_transcoef_data* data, Options* options,
                              MarkerFieldLine* p_f, MarkerFieldLine* p_i);
#endif
