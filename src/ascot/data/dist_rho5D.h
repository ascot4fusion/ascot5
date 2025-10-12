/**
 * @file dist_rho5D.h
 * @brief Header file for dist_rho5D.c
 */
#ifndef DIST_RHO5D_H
#define DIST_RHO5D_H

#include <stdlib.h>
#include "defines.h"
#include "marker.h"
#include "diag.h"



int dist_rho5D_init(dist_rho5D_data* data);
void dist_rho5D_free(dist_rho5D_data* data);
void dist_rho5D_offload(dist_rho5D_data* data);
void dist_rho5D_update_fo(dist_rho5D_data* dist, MarkerGyroOrbit* p_f,
                          MarkerGyroOrbit* p_i);
void dist_rho5D_update_gc(dist_rho5D_data* dist, MarkerGuidingCenter* p_f,
                          MarkerGuidingCenter* p_i);

#endif
