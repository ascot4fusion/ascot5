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
#include "wall/wall_2d.h"
#include "wall/wall_3d.h"

/**
 * @brief Wall model types
 *
 * Wall model types are used in the magnetic wall interface (wall.c) to direct
 * function calls to correct wall model instances. Each wall model instance must
 * have a corresponding type.
 */
typedef enum wall_type {
    wall_type_2D, /**< Axisymmetric wall model consisting of single contour */
    wall_type_3D, /**< 3D wall model consisting of triangles                */
} wall_type;

/**
 * @brief Wall model offload data
 *
 * This struct holds data necessary for offloading. The struct is initialized in
 * wall_init_offload().
 *
 * The intended usage is that only single offload data is used at the time, and
 * the type of the data is declared with the "type" field.
 */
typedef struct {
    wall_type type;               /**< Wall model type wrapped by this struct */
    wall_2d_offload_data w2d;     /**< 2D model or NULL if not active         */
    wall_3d_offload_data w3d;     /**< 3D model or NULL if not active         */
    int offload_array_length;     /**< Allocated offload array length         */
    int int_offload_array_length; /**< Allocated int offload array length     */
} wall_offload_data;

/**
 * @brief Wall model simulation data
 *
 * This struct holds data necessary for simulation. The struct is initialized
 * from the wall_offload_data in wall_init().
 *
 * The intended usage is that only single wall_data is used at the time, and
 * the type of the data is declared with the "type" field.
 */
typedef struct {
    wall_type type;   /**< Wall model type wrapped by this struct */
    wall_2d_data w2d; /**< 2D model or NULL if not active         */
    wall_3d_data w3d; /**< 3D model or NULL if not active         */
} wall_data;

int wall_init_offload(wall_offload_data* offload_data, real** offload_array,
                      int** int_offload_array);
void wall_free_offload(wall_offload_data* offload_data, real** offload_array,
                       int** int_offload_array);

#pragma omp declare target
int wall_init(wall_data* w, wall_offload_data* offload_data,
              real* offload_array, int* int_offload_array);
#pragma omp declare simd uniform(w)
int wall_hit_wall(real r1, real phi1, real z1, real r2, real phi2, real z2,
                  wall_data* w, real* w_coll);
#pragma omp declare simd uniform(w)
int wall_get_n_elements(wall_data* w);
#pragma omp end declare target

#endif
