/**
 * @file wall.c
 * @brief Wall interface
 *
 * This is an interface through which wall data is initialized and accessed.
 * Reading e.g. from disk is done elsewhere.
 *
 * To add a new wall instance, make sure these functions are implemented and
 * called from this interface, and that wall.h contains enum type for the new
 * instance.
 *
 * The interface checks which instance given data corresponds to from the
 * "type"-field in wall_offload_data or wall_data that is given as an
 * argument, and calls the relevant function for that instance.
 */
#include <stdio.h>
#include <math.h>
#include "ascot5.h"
#include "print.h"
#include "wall.h"
#include "wall/wall_2d.h"
#include "wall/wall_3d.h"

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void wall_free(wall_data* data) {
    switch(data->type) {
        case wall_type_2D:
            wall_2d_free(&data->w2d);
            break;

        case wall_type_3D:
            wall_3d_free(&data->w3d);
            break;
    }
}

/**
 * @brief Check if a given directed line segment intersects the wall
 *
 * This function is intended to be used to check whether a marker collides with
 * the wall. If there is a collision, this function returns an identification
 * number specific to that wall tile. If the marker hits multiple wall elements,
 * only the first one is returned.
 *
 * This is a SIMD function.
 *
 * @param r1 start point R coordinate [m]
 * @param phi1 start point phi coordinate [rad]
 * @param z1 start point z coordinate [rad]
 * @param r2 end point R coordinate [m]
 * @param phi2 end point phi coordinate [rad]
 * @param z2 end point z coordinate [rad]
 * @param w pointer to data struct on target
 * @param w_coll pointer for storing the parameter in P = P1 + w_coll * (P2-P1),
 *        where P is the point where the collision occurred.
 *
 * @return wall element id if hit, zero otherwise
 */
int wall_hit_wall(real r1, real phi1, real z1, real r2, real phi2, real z2,
                  wall_data* w, real* w_coll) {
    int ret = 0;
    switch(w->type) {
        case wall_type_2D:
            ret = wall_2d_hit_wall(r1, phi1, z1, r2, phi2, z2, &(w->w2d),
                                   w_coll);
            break;

        case wall_type_3D:
            ret = wall_3d_hit_wall(r1, phi1, z1, r2, phi2, z2, &(w->w3d),
                                   w_coll);
            break;
    }
    return ret;
}

/**
 * @brief Return the number of wall elements.
 *
 * @param w pointer to wall data struct on target
 *
 * @return Number of wall elements or zero on failure.
 */
int wall_get_n_elements(wall_data* w) {
    int ret = 0;
    switch(w->type) {
        case wall_type_2D:
            ret = w->w2d.n;
            break;

        case wall_type_3D:
            ret = w->w3d.n;
            break;
    }
    return ret;
}
