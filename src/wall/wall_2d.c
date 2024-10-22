/**
 * @file wall_2d.c
 * @brief 2D wall collision checks
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "../print.h"
#include "wall_2d.h"

/**
 * @brief Load 2D wall data and prepare parameters
 *
 * @param data pointer to the data struct
 * @param nelements number of elements in the wall polygon
 * @param r R coordinates for the wall polygon points
 * @param z z coordinates for the wall polygon points
 *
 * @return zero to indicate success
 */
int wall_2d_init(wall_2d_data* data, int nelements, real* r, real* z) {

    data->n = nelements;
    data->wall_r = (real*) malloc( nelements * sizeof(real) );
    data->wall_z = (real*) malloc( nelements * sizeof(real) );
    real rmin = r[0], rmax = r[0];
    real zmin = z[0], zmax = z[0];
    for(int i=0; i < nelements; i++) {
        rmin = fmin(rmin, r[i]);
        rmax = fmax(rmax, r[i]);
        zmin = fmin(zmin, z[i]);
        zmax = fmax(zmax, z[i]);
        data->wall_r[i] = r[i];
        data->wall_z[i] = z[i];
    }

    print_out(VERBOSE_IO, "\n2D wall model (wall_2D)\n");
    print_out(VERBOSE_IO, "Number of wall elements = %d,"
              " R extend = [%2.2f, %2.2f], z extend = [%2.2f, %2.2f]\n",
              nelements, rmin, rmax, zmin, zmax);
    return 0;
}

/**
 * @brief Free allocated resources
 *
 * @param offload_data pointer to the data struct
 */
void wall_2d_free(wall_2d_data* data) {
    free(data->wall_r);
    free(data->wall_z);
}

/**
 * @brief Check if coordinates are within 2D polygon wall
 *
 * This function checks if the given coordinates are within the walls defined
 * by a 2D polygon using a modified axis crossing method [1]. Origin is moved
 * to the coordinates and the number of wall segments crossing the positive
 * r-axis are calculated. If this is odd, the point is inside the polygon.
 *
 * [1] D.G. Alciatore, R. Miranda. A Winding Number and Point-in-Polygon
 *     Algorithm. Technical report, Colorado State University, 1995.
 *     https://www.engr.colostate.edu/~dga/documents/papers/point_in_polygon.pdf
 *
 * @param r r coordinate [m]
 * @param z z coordinate [m]
 * @param w 2D wall data structure
 */
int wall_2d_inside(real r, real z, wall_2d_data* w) {
    int hits = 0;
    for(int i = 0; i < w->n; i++) {
        real wr1, wr2, wz1, wz2;
        if(i == w->n - 1) {
            wz1 = w->wall_z[i] - z;
            wz2 = w->wall_z[0]   - z;
            wr1 = w->wall_r[i] - r;
            wr2 = w->wall_r[0]   - r;
        } else {
            wz1 = w->wall_z[i]   - z;
            wz2 = w->wall_z[i+1] - z;
            wr1 = w->wall_r[i]   - r;
            wr2 = w->wall_r[i+1] - r;
        }
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
 *
 * @param r1 start point R coordinate [m]
 * @param phi1 start point phi coordinate [rad]
 * @param z1 start point z coordinate [m]
 * @param r2 end point R coordinate [m]
 * @param phi2 end point phi coordinate [rad]
 * @param z2 end point z coordinate [m]
 * @param w pointer to data struct on target
 * @param w_coll pointer for storing the parameter in P = P1 + w_coll * (P2-P1),
 *        where P is the point where the collision occurred.
 *
 * @return wall element ID if hit, zero otherwise
 */
int wall_2d_hit_wall(real r1, real phi1, real z1, real r2, real phi2, real z2,
                     wall_2d_data* w, real* w_coll) {
    int tile = 0;
    if(!wall_2d_inside(r2, z2, w)) {
        tile = wall_2d_find_intersection(r1, z1, r2, z2, w, w_coll);
    }
    return tile;
}

/**
 * @brief Find intersection between the wall element and line segment
 *
 * If there are multiple intersections, the one that is closest to P1
 * is returned.
 *
 * @param r1 R1 coordinate of the line segment [P1,P2] [m]
 * @param z1 z1 coordinate of the line segment [P1,P2] [m]
 * @param r2 R2 coordinate of the line segment [P1,P2] [m]
 * @param z2 z2 coordinate of the line segment [P1,P2] [m]
 * @param w pointer to the wall data
 * @param w_coll pointer for storing the parameter in P = P1 + w_coll * (P2-P1),
 *        where P is the point where the collision occurred.
 *
 * @return int wall element id if hit, zero otherwise
 */
int wall_2d_find_intersection(real r1, real z1, real r2, real z2,
                              wall_2d_data* w, real* w_coll) {
    int tile = 0;
    real t0 = 2.0; // Helper variable to pick the closest intersection
    for(int i=0; i<w->n; i++) {
        real r3, z3, r4, z4;
        if(i == w->n-1) {
            r3 = w->wall_r[i];
            z3 = w->wall_z[i];
            r4 = w->wall_r[0];
            z4 = w->wall_z[0];
        } else {
            r3 = w->wall_r[i];
            z3 = w->wall_z[i];
            r4 = w->wall_r[i+1];
            z4 = w->wall_z[i+1];
        }

        real div = (r1 - r2) * (z3 - z4) - (z1 - z2) * (r3 - r4);
        real t   = ( (r1 - r3) * (z3 - z4) - (z1 - z3) * (r3 - r4) ) / div;
        real u   = ( (r1 - r3) * (z1 - z2) - (z1 - z3) * (r1 - r2) ) / div;
        if(0 <= t && t <= 1.0 && 0 <= u && u <= 1.0 && t < t0) {
            t0   = t;
            tile = i + 1;
        }
    }
    *w_coll = t0;
    return tile;
}
