/**
 * @file wall_2d.h
 * @brief Header file for wall_2d.c
 */
#ifndef WALL_2D_H
#define WALL_2D_H
#include "../ascot5.h"
#include "../offload.h"

/**
 * @brief 2D wall data parameters
 *
 * Note: The start and end point of wall polygon does not have to concide.
 */
typedef struct {
    int n;        /**< Number of points in the wall polygon      */
    real* wall_r; /**< R coordinates for the wall polygon points */
    real* wall_z; /**< z coordinates for the wall polygon points */
    int* flag;    /**< Array of wall element flags               */
} wall_2d_data;

int wall_2d_init(wall_2d_data* data, int nelements, real* r, real* z,
                 int* flag);
void wall_2d_free(wall_2d_data* data);
void wall_2d_offload(wall_2d_data* data);
GPU_DECLARE_TARGET_SIMD_UNIFORM(w)
int wall_2d_inside(real r, real z, wall_2d_data* w);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(w)
int wall_2d_hit_wall(real r1, real phi1, real z1, real r2, real phi2, real z2,
                     wall_2d_data* w, real* w_coll);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(w)
int wall_2d_find_intersection(real r1, real z1, real r2, real z2,
                              wall_2d_data* w, real* w_coll);
DECLARE_TARGET_END
#endif
