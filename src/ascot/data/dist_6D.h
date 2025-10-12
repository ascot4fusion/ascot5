/**
 * @file dist_6D.h
 * @brief Header file for dist_6D.c
 */
#ifndef DIST_6D_H
#define DIST_6D_H

#include <stdlib.h>
#include "defines.h"
#include "marker.h"
#include "diag.h"



int dist_6D_init(dist_6D_data* data);
void dist_6D_free(dist_6D_data* data);
void dist_6D_offload(dist_6D_data* data);
void dist_6D_update_fo(dist_6D_data* dist,
                       MarkerGyroOrbit* p_f, MarkerGyroOrbit* p_i);
void dist_6D_update_gc(dist_6D_data* dist,
                       MarkerGuidingCenter* p_f, MarkerGuidingCenter* p_i);

#endif
