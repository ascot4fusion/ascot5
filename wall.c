/**
 * @file wall.c
 * @brief Wall interface
 */
#include <stdio.h>
#include "wall.h"
#include "wall_2d.h"
#include "wall_3d.h"

void wall_init_offload(wall_offload_data* offload_data, real** offload_array) {
    FILE* f = fopen("input.wall_3d", "r");
    if(f != NULL) {
        /* 3d wall input found */
        offload_data->type = 3;
    } else {
        /* no 3d wall input, assume 2d wall */
        offload_data->type = 1;
    }

    switch(offload_data->type) {
        case 1:
        wall_2d_init_offload(&(offload_data->w2d), offload_array);
        offload_data->offload_array_length = offload_data->w2d.offload_array_length;
        break;

        case 3:
        wall_3d_init_offload(&(offload_data->w3d), offload_array);
        offload_data->offload_array_length = offload_data->w3d.offload_array_length;
        break;
    }
}

void wall_free_offload(wall_offload_data* offload_data, real** offload_array) {
    switch(offload_data->type) {
        case 1:
        wall_2d_free_offload(&(offload_data->w2d), offload_array);
        break;

        case 3:
        wall_3d_free_offload(&(offload_data->w3d), offload_array);
        break;
    }
}

void wall_init(wall_data* w, wall_offload_data* offload_data,
                  real* offload_array) {
    switch(offload_data->type) {
        case 1:
        wall_2d_init(&(w->w2d), &(offload_data->w2d), offload_array);
        break;

        case 3:
        wall_3d_init(&(w->w3d), &(offload_data->w3d), offload_array);
        break;
    }
    w->type = offload_data->type;
}

int wall_hit_wall(real r1, real phi1, real z1, real r2, real phi2, real z2,
                  wall_data* w) {
    switch(w->type) {
        case 1:
        return wall_2d_hit_wall(r1, phi1, z1, r2, phi2, z2, &(w->w2d));
        break;

        case 3:
        return wall_3d_hit_wall(r1, phi1, z1, r2, phi2, z2, &(w->w3d));
        break;
    }
}
