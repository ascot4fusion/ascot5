/**
 * @file wall.h
 * @brief Header file for wall.c
 */
#ifndef WALL_H
#define WALL_H

#include "wall/wall_2d.h"
#include "wall/wall_3d.h"

typedef struct {
    int type;
    wall_2d_offload_data w2d;
    wall_3d_offload_data w3d;
    int offload_array_length;
} wall_offload_data;

typedef struct {
    int type;
    wall_2d_data w2d;
    wall_3d_data w3d;
} wall_data;

void wall_init_offload(wall_offload_data* offload_data, real** offload_array);
void wall_free_offload(wall_offload_data* offload_data, real** offload_array);

#pragma omp declare target
void wall_init(wall_data* w, wall_offload_data* offload_data,
               real* offload_array);
#pragma omp declare simd uniform(w)
int wall_hit_wall(real r1, real phi1, real z1, real r2, real phi2, real z2,
                  wall_data* w);
#pragma omp end declare target

#endif
