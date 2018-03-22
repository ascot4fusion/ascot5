/**
 * @file wall_2d.c
 * @brief 2D wall collision checks
 */

#include <stdio.h>
#include <stdlib.h>
#include "../ascot5.h"
#include "wall_2d.h"

/**
 * @brief Load 2D wall data and prepare parameters
 *
 * Read a 2D wall from input.wall_2d, store parameters in offload struct and
 * allocate and fill the offload array.
 *
 * @todo Error checking
 * @todo Move reading the file to ascot4_interface
 * 
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void wall_2d_init_offload(wall_2d_offload_data* offload_data,
                          real** offload_array) {
    // Dummy function
}

/**
 * @brief Free offload array and reset parameters
 */
void wall_2d_free_offload(wall_2d_offload_data* offload_data,
                          real** offload_array) {
    free(*offload_array);
    *offload_array = NULL;
}

/**
 * @brief Initialize 2D wall data struct on target
 *
 * Copies the 2D wall parameters from the offload struct to the struct on
 * on target and sets the wall data pointer inside to the offload array.
 */
void wall_2d_init(wall_2d_data* w, wall_2d_offload_data* offload_data,
                  real* offload_array) {
    w->n = offload_data->n;
    w->wall_r = &offload_array[0];
    w->wall_z = &offload_array[offload_data->n];
}

/**
 * @brief Check if coordinates are within 2D polygon wall
 *
 * This function checks if the given coordinates are within the walls defind
 * by a 2D polygon using a modified axis crossing method [1]. Origin is moved
 * to the coordinates and the number of wall segments crossing the positive
 * r-axis are calculated. If this is odd, the point is inside the polygon.
 *
 * [1] D.G. Alciatore, R. Miranda. A Winding Number and Point-in-Polygon
 *     Algorithm. Technical report, Colorado State University, 1995.
 *     http://www.engr.colostate.edu/~dga/dga/papers/point_in_polygon.pdf
 *
 * @param r r coordinate
 * @param z z coordinate
 * @param w 2D wall data structure
 */
int wall_2d_inside(real r, real z, wall_2d_data* w) {
    int hits = 0;

    int i;
    for(i = 0; i < w->n - 1; i++) {
        real wz1 = w->wall_z[i] - z;
        real wz2 = w->wall_z[i+1] - z;
        real wr1 = w->wall_r[i] - r;
        real wr2 = w->wall_r[i+1] - r;
        if(wz1 * wz2 < 0) {
            real ri = wr1 + (wz1*(wr2-wr1)) / (wz1-wz2);
            if(ri > 0) {
                hits++;
            }
        }
    }
    return hits % 2;
}

/**
 * @brief Check if trajectory from (r1, phi1, z1) to (r2, phi2, z2) intersects
 *        the wall
 */
int wall_2d_hit_wall(real r1, real phi1, real z1, real r2, real phi2, real z2,
                     wall_2d_data* w) {
    if(!wall_2d_inside(r2, z2, w))
        return 1;
    else
        return 0;
}
