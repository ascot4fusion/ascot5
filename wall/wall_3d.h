/**
 * @file wall_3d.h
 * @brief Header file for wall_3d.c
 */
#ifndef WALL_3D_H
#define WALL_3D_H
#include "../ascot5.h"

/** Default depth of octree struct */
#define WALL_OCTREE_DEPTH 7

/**
 * @brief 3D wall offload data
 */
typedef struct {
    int n;                    /**< Number of wall triangles                   */
    real xmin;                /**< Minimum extend on x-direction [m]          */
    real xmax;                /**< Maximum extend on x-direction [m]          */
    real xgrid;               /**< Octree cell width in x-direction [m]       */
    real ymin;                /**< Minimum extend on y-direction [m]          */
    real ymax;                /**< Maximum extend on y-direction [m]          */
    real ygrid;               /**< Octree cell width in y-direction [m]       */
    real zmin;                /**< Minimum extend on z-direction [m]          */
    real zmax;                /**< Maximum extend on z-direction [m]          */
    real zgrid;               /**< Octree cell width in z-direction [m]       */
    int depth;                /**< Depth of the octree                        */
    int ngrid;                /**< Number of cells computational volume is
                                   divided to in each direction.
                                   ngrid = 2^(depth-1)                        */
    int offload_array_length; /**< Length of the offload array                */
} wall_3d_offload_data;

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
    int tree_array_size; /**<  */
    int* tree_array;     /**< Pointer to octree array */
} wall_3d_data;
  
int wall_3d_init_offload(wall_3d_offload_data* offload_data,
                         real** offload_array);
void wall_3d_free_offload(wall_3d_offload_data* offload_data,
                          real** offload_array);

#pragma omp declare target
void wall_3d_init(wall_3d_data* w, wall_3d_offload_data* offload_data,
                  real* offload_array);
#pragma omp declare simd uniform(w)
int wall_3d_hit_wall(real r1, real phi1, real z1, real r2, real phi2,
                     real z2, wall_3d_data* w, real* w_coll);
#pragma omp declare simd uniform(w)
int wall_3d_hit_wall_full(real r1, real phi1, real z1, real r2, real phi2,
                          real z2, wall_3d_data* w, real* w_coll);
#pragma omp declare simd
double wall_3d_tri_collision(real q1[3], real q2[3], real t1[3], real t2[3],
                             real t3[3]);

void wall_3d_init_tree(wall_3d_data* w, real* offload_array);
void wall_3d_init_octree(wall_3d_data* w, real* offload_array);
int wall_3d_tri_in_cube(real t1[3], real t2[3], real t3[3], real bb1[3],
                        real bb2[3]);
int wall_3d_quad_collision(real q1[3], real q2[3], real t1[3], real t2[3],
                           real t3[3], real t4[3]);
#pragma omp end declare target

#endif
