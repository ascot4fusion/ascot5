/**
 * @file wall_2d.h
 * @brief Header file for wall_2d.c
 */
#ifndef WALL_2D_H
#define WALL_2D_H
#include "ascot5.h"

typedef struct {
    int n;
    int offload_array_length;
} wall_2d_offload_data;

/**
 * @brief 2D wall data parameters
 */
typedef struct {
	int n;          /**< number of points in the wall polygon */
	real* wall_r;   /**< r coordinates for the wall polygon points */
	real* wall_z;   /**< z coordinates for the wall polygon points */
} wall_2d_data;

void wall_2d_init_offload(wall_2d_offload_data* offload_data,
                          real** offload_array);
void wall_2d_free_offload(wall_2d_offload_data* offload_data,
                          real** offload_array);

#pragma omp declare target
void wall_2d_init(wall_2d_data* w, wall_2d_offload_data* offload_data,
                  real* offload_array);
#pragma omp declare simd uniform(w)
int wall_2d_inside(real r, real z, wall_2d_data* w);
#pragma omp declare simd uniform(w)
int wall_2d_hit_wall(real r1, real phi1, real z1, real r2, real phi2, real z2,
                     wall_2d_data* w);
#pragma omp end declare target

#endif
