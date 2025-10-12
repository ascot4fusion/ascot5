/**
 * @file wall.h
 * Interface through which wall data is initialized and accessed.
 */
#ifndef WALL_H
#define WALL_H

#include "defines.h"
#include "parallel.h"
#include "wall_data.h"

/**
 * Free allocated resources.
 *
 * @param wall Pointer to the data struct.
 */
void Wall_free(Wall *wall);

/**
 *  Offload data to the accelerator.
 *
 * @param wall pointer to the data struct.
 */
void Wall_offload(Wall *wall);

/**
 * Check if a given directed line segment intersects the wall.
 *
 * This function is intended to be used to check whether a marker collides with
 * the wall. If there is a collision, this function returns an identification
 * number specific to that wall tile. If the marker hits multiple wall elements,
 * only the first one is returned.
 *
 * This is a SIMD function.
 *
 * @param r1 Start point R coordinate [m].
 * @param phi1 Start point phi coordinate [rad].
 * @param z1 Start point z coordinate [rad].
 * @param r2 End point R coordinate [m].
 * @param phi2 End point phi coordinate [rad].
 * @param z2 End point z coordinate [rad].
 * @param wall Pointer to data struct on target.
 * @param w_coll Pointer for storing the parameter in P = P1 + w_coll * (P2-P1),
 *        where P is the point where the collision occurred.
 *
 * @return Wall element id if hit, zero otherwise.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(wall)
size_t Wall_hit_wall(
    real r1, real phi1, real z1, real r2, real phi2, real z2, real *w_coll,
    Wall *wall);
DECLARE_TARGET_END

/**
 * Return the number of wall elements.
 *
 * @param wall pointer to wall data struct on target.
 *
 * @return Number of wall elements or zero on failure.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(wall)
size_t Wall_get_n_elements(Wall *wall);
DECLARE_TARGET_END

/**
 * Return the flag of a wall element.
 *
 * @param wall Pointer to wall data struct on target.
 * @param idx Wall element index.
 *
 * @return Flag of the wall element.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(wall)
int Wall_get_flag(size_t idx, Wall *wall);
DECLARE_TARGET_END

#endif
