/**
 * @file dist_rho6D.h
 * @brief Header file for dist_rho6D.c
 */
#ifndef DIST_RHO6D_H
#define DIST_RHO6D_H

#include <stdlib.h>
#include "defines.h"
#include "marker.h"
#include "diag.h"



int dist_rho6D_init(dist_rho6D_data* dist_data);
void dist_rho6D_free(dist_rho6D_data* dist_data);
void dist_rho6D_offload(dist_rho6D_data* dist_data);
void dist_rho6D_update_fo(dist_rho6D_data* dist, MarkerGyroOrbit* p_f,
                          MarkerGyroOrbit* p_i);
void dist_rho6D_update_gc(dist_rho6D_data* dist, MarkerGuidingCenter* p_f,
                          MarkerGuidingCenter* p_i);

#endif
