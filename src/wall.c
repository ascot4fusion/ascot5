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
<<<<<<< HEAD
void wall_free(wall_data* data) {
    switch(data->type) {
        case wall_type_2D:
            wall_2d_free(&data->w2d);
            break;

        case wall_type_3D:
            wall_3d_free(&data->w3d);
=======
int wall_init_offload(wall_offload_data* offload_data, real** offload_array,
                      int** int_offload_array) {
    int err = 0;
    switch(offload_data->type) {

        case wall_type_2D:
            err = wall_2d_init_offload(&(offload_data->w2d), offload_array,
                                       int_offload_array);
            offload_data->offload_array_length =
                offload_data->w2d.offload_array_length;
            offload_data->int_offload_array_length =
                offload_data->w2d.int_offload_array_length;
            break;

        case wall_type_3D:
            err = wall_3d_init_offload(&(offload_data->w3d), offload_array,
                                       int_offload_array);
            offload_data->offload_array_length =
                offload_data->w3d.offload_array_length;
            offload_data->int_offload_array_length =
                offload_data->w3d.int_offload_array_length;
            break;

        default:
            /* Unregonized input. Produce error. */
            print_err("Error: Unregonized electric field type.");
            err = 1;
            break;
    }
    if(!err) {
        print_out(VERBOSE_IO, "Estimated memory usage %.1f MB\n",
                  offload_data->offload_array_length
                  * sizeof(real) / (1024.0*1024.0) );
    }

    return err;
}

/**
 * @brief Free offload array and reset parameters
 *
 * This function deallocates the offload_array.
 *
 * This function is host only.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 * @param int_offload_array pointer to pointer to offload array storing integers
 */
void wall_free_offload(wall_offload_data* offload_data, real** offload_array,
                       int** int_offload_array) {
    switch(offload_data->type) {
        case wall_type_2D:
            wall_2d_free_offload(&(offload_data->w2d), offload_array,
                                 int_offload_array);
            break;

        case wall_type_3D:
            wall_3d_free_offload(&(offload_data->w3d), offload_array,
                                 int_offload_array);
>>>>>>> 3365d39a (Added flags to 2D wall and flags are now read from the HDF5 and can be used in a simulation)
            break;
    }
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void wall_offload(wall_data* data) {
    switch(data->type) {
        case wall_type_2D:
<<<<<<< HEAD
            wall_2d_offload(&data->w2d);
=======
            wall_2d_init(&(w->w2d), &(offload_data->w2d), offload_array,
                         int_offload_array);
>>>>>>> 3365d39a (Added flags to 2D wall and flags are now read from the HDF5 and can be used in a simulation)
            break;

        case wall_type_3D:
            wall_3d_offload(&data->w3d);
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

/**
 * @brief Return the flag of a wall element.
 *
 * @param w pointer to wall data struct on target
 * @param idx wall element index
 *
 * @return Flag of the wall element.
 */
int wall_get_flag(wall_data* w, int idx) {
    int flag = 0;
    switch(w->type) {
        case wall_type_2D:
            flag = w->w2d.flag[idx];
            break;

        case wall_type_3D:
            flag = w->w3d.flag[idx];
            break;
    }
    return flag;
}