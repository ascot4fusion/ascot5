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
 * @brief Load wall data and prepare parameters
 *
 * This function fills the relevant wall offload struct with parameters and
 * allocates and fills the offload array. Sets offload array length in the
 * offload struct.
 *
 * The offload data has to have a type when this function is called as it should
 * be set when the offload data is constructed from inputs.
 *
 * This function is host only.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @return zero if initialization succeeded
 */
int wall_init_offload(wall_offload_data* offload_data, real** offload_array) {

    int err = 0;

    switch(offload_data->type) {

        case wall_type_2D:
            err = wall_2d_init_offload(&(offload_data->w2d), offload_array);
            offload_data->offload_array_length =
                offload_data->w2d.offload_array_length;
            break;

        case wall_type_3D:
            err = wall_3d_init_offload(&(offload_data->w3d), offload_array);
            offload_data->offload_array_length =
                offload_data->w3d.offload_array_length;
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
 */
void wall_free_offload(wall_offload_data* offload_data, real** offload_array) {
    switch(offload_data->type) {
        case wall_type_2D:
            wall_2d_free_offload(&(offload_data->w2d), offload_array);
            break;

        case wall_type_3D:
            wall_3d_free_offload(&(offload_data->w3d), offload_array);
            break;
    }
}

/**
 * @brief Initialize wall data struct on target
 *
 * This function copies the wall parameters from the offload struct to the
 * struct on target and sets the wall data pointers to correct offsets in the
 * offload array.
 *
 * @param w pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 *
 * @return zero on success
 */
int wall_init(wall_data* w, wall_offload_data* offload_data,
              real* offload_array) {
    int err = 0;
    switch(offload_data->type) {
        case wall_type_2D:
            wall_2d_init(&(w->w2d), &(offload_data->w2d), offload_array);
            break;

        case wall_type_3D:
            wall_3d_init(&(w->w3d), &(offload_data->w3d), offload_array);
            break;
        default:
            /* Unregonized input. Produce error. */
            err = 1;
            break;
    }
    w->type = offload_data->type;

    return err;
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
 *
 * @return wall element id if hit, zero otherwise
 */
int wall_hit_wall(real r1, real phi1, real z1, real r2, real phi2, real z2,
                  wall_data* w, real* w_coll) {
    int ret = 0;
    switch(w->type) {
        case wall_type_2D:
            ret = wall_2d_hit_wall(r1, phi1, z1, r2, phi2, z2, &(w->w2d));
            break;

        case wall_type_3D:
	    ret = wall_3d_hit_wall(r1, phi1, z1, r2, phi2, z2, &(w->w3d), w_coll);
	    break;
    }
    return ret;
}
