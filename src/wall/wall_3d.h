/**
 * @file wall_3d.h
 * @brief Header file for wall_3d.c
 */
#ifndef WALL_3D_H
#define WALL_3D_H
#include "../ascot5.h"
#include "../offload.h"

/** Small value to check if x = 0 (i.e. abs(x) < WALL_EPSILON) */
#define WALL_EPSILON 1e-9

/**
 * @brief 3D wall data parameters
 */
typedef struct {
    int n;               /**< Number of wall triangles                        */
    real xmin;           /**< Minimum extend on x-direction [m]               */
    real xmax;           /**< Maximum extend on x-direction [m]               */
    real xgrid;          /**< Octree cell width in x-direction [m]            */
    real ymin;           /**< Minimum extend on y-direction [m]               */
    real ymax;           /**< Maximum extend on y-direction [m]               */
    real ygrid;          /**<  Octree cell width in y-direction [m]           */
    real zmin;           /**< Minimum extend on z-direction [m]               */
    real zmax;           /**< Maximum extend on z-direction [m]               */
    real zgrid;          /**< Octree cell width in z-direction [m]            */
    int depth;           /**< Depth of the octree                             */
    int ngrid;           /**< Number of cells computational volume is divided
                              to in each direction. ngrid = 2^(depth-1)       */
    real* wall_tris;     /**< Array of wall triangle coordinates */
    int* flag;           /**< Array of wall triangle flags                    */

    /** @brief Array storing information what triangles given octree cell stores
     *
     * First ncell elements store the array position where data for a given cell
     * begins, where the cell index is icell = ix * ngrid**2 + iy * ngrid + iz.
     * The first element in the cell data, i.e. tree_array[tree_array[icell]],
     * is the number of triangles in this cell, ntriangle, and the next
     * ntriangle elements are the triangle indices.
     */
    int* tree_array;
    int tree_array_size; /* Number of elements in tree_array */
} wall_3d_data;

int wall_3d_init(wall_3d_data* data, int nelements, real* vertices, int* flag);
void wall_3d_free(wall_3d_data* data);
void wall_3d_offload(wall_3d_data* data);
GPU_DECLARE_TARGET_SIMD_UNIFORM(w)
int wall_3d_hit_wall(real r1, real phi1, real z1, real r2, real phi2,
                     real z2, wall_3d_data* w, real* w_coll);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(w)
int wall_3d_hit_wall_full(real r1, real phi1, real z1, real r2, real phi2,
                          real z2, wall_3d_data* w, real* w_coll);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD
double wall_3d_tri_collision(real q1[3], real q2[3], real t1[3], real t2[3],
                             real t3[3]);
DECLARE_TARGET_END

DECLARE_TARGET
void wall_3d_init_tree(wall_3d_data* w, real* offload_array);
DECLARE_TARGET_END
DECLARE_TARGET
int wall_3d_tri_in_cube(real t1[3], real t2[3], real t3[3], real bb1[3],
                        real bb2[3]);
DECLARE_TARGET_END
DECLARE_TARGET
int wall_3d_quad_collision(real q1[3], real q2[3], real t1[3], real t2[3],
                           real t3[3], real t4[3]);
DECLARE_TARGET_END

#endif
