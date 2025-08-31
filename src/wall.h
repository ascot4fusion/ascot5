/**
 * @file wall.h
 * @brief Header file for wall.c
 *
 * Contains a list declaring all wall_types, and declaration of
 * wall_offload_data and wall_data structs.
 *
 */
#ifndef WALL_H
#define WALL_H

#include "ascot5.h"
#include "offload.h"
#include "wall/wall_2d.h"
#include "wall/wall_3d.h"

/**
 * @brief Wall model types
 */
typedef enum wall_type {
    wall_type_2D, /**< Axisymmetric wall model consisting of single contour */
    wall_type_3D, /**< 3D wall model consisting of triangles                */
} wall_type;

/**
 * @brief Wall model simulation data
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the `type` field.
 */
typedef struct {
    wall_2d_data* w2d; /**< 2D model or NULL if not active         */
    wall_3d_data* w3d; /**< 3D model or NULL if not active         */
    wall_type type;    /**< Wall model type wrapped by this struct */
} wall_data;

void wall_free(wall_data* data);
void wall_offload(wall_data* data);
GPU_DECLARE_TARGET_SIMD_UNIFORM(w)
int wall_hit_wall(real r1, real phi1, real z1, real r2, real phi2, real z2,
                  wall_data* w, real* w_coll);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(w)
int wall_get_n_elements(wall_data* w);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(w)
int wall_get_flag(wall_data* w, int idx);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(w)
int wall_get_n_elements(wall_data* w);
DECLARE_TARGET_END
#endif
