/**
 * @file dist_com.h
 * @brief Header file for dist_com.c
 */
#ifndef DIST_COM_H
#define DIST_COM_H

#include <stdlib.h>
#include "defines.h"
#include "marker.h"
#include "bfield.h"
#include "diag.h"



int dist_COM_init(dist_COM_data* data);
void dist_COM_free(dist_COM_data* data);
void dist_COM_offload(dist_COM_data* data);
void dist_COM_update_fo(dist_COM_data* dist, Bfield *bfield,
                        MarkerGyroOrbit* p_f, MarkerGyroOrbit* p_i);
void dist_COM_update_gc(dist_COM_data* dist, Bfield *bfield,
                        MarkerGuidingCenter* p_f, MarkerGuidingCenter* p_i);

#endif
