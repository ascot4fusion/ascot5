/**
 *  * @file nbi_diag.h
 *   * @brief Header file for nbi_diag.c
 *    */
#ifndef NBI_DIAG_H
#define NBI_DIAG_H

#include "ascot5.h"

typedef struct {
    int n_x;
    real min_x;
    real max_x;

    int n_y;
    real min_y;
    real max_y;

    int n_z;
    real min_z;
    real max_z;

    real full_energy;

    real* hist_neutral;
} nbi_diag_data;

void nbi_diag_init(nbi_diag_data* data, int n_x, real min_x, real max_x,
                   int n_y, real min_y, real max_y, int n_z, real min_z,
                   real max_z, real full_energy);
void nbi_diag_update_dist(nbi_diag_data* data, real x, real y, real z, real energy, real weight);
void nbi_diag_normalize_dist(nbi_diag_data* data, real weight);

#endif
