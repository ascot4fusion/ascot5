/**
 * @file wall_3d.h
 * @brief Header file for wall_3d.c
 */
#ifndef WALL_3D_H
#define WALL_3D_H
#include "../ascot5.h"

typedef struct {
    int n;
    real xmin;
    real xmax;
    real xgrid;
    real ymin;
    real ymax;
    real ygrid;
    real zmin;
    real zmax;
    real zgrid;
    int depth;
    int ngrid;
    int offload_array_length;
} wall_3d_offload_data;

/**
 * @brief 3D wall data parameters
 */
typedef struct {
    int n;          /**< number of points in the wall polygon */
    real xmin;
    real xmax;
    real xgrid;
    real ymin;
    real ymax;
    real ygrid;
    real zmin;
    real zmax;
    real zgrid;
    int depth;
    int ngrid;
    real* wall_tris;
    int tree_array_size;
    int* tree_array;
} wall_3d_data;

void wall_3d_init_offload(wall_3d_offload_data* offload_data,
                          real** offload_array);
void wall_3d_free_offload(wall_3d_offload_data* offload_data,
                          real** offload_array);

#pragma omp declare target
void wall_3d_init(wall_3d_data* w, wall_3d_offload_data* offload_data,
                  real* offload_array);
#pragma omp declare simd uniform(w)
int wall_3d_hit_wall(real r1, real phi1, real z1, real r2, real phi2,
                     real z2, wall_3d_data* w);
#pragma omp declare simd uniform(w)
int wall_3d_hit_wall_full(real r1, real phi1, real z1, real r2, real phi2,
                          real z2, wall_3d_data* w);
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
