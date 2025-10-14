/**
 * @file wall.h
 * Wall interface.
 *
 * Provides functions and datatypes for checking intersections between marker
 * trajectories and wall elements.
 */
#ifndef WALL_H
#define WALL_H

#include "defines.h"
#include "parallel.h"
#include "wall_data.h"

/**
 * Free allocated resources.
 *
 * @param wall The struct whose fields are deallocated.
 */
void Wall_free(Wall *wall);

/**
 * Offload data to the accelerator.
 *
 * @param wall The struct to offload.
 */
void Wall_offload(Wall *wall);

GPU_DECLARE_TARGET_SIMD_UNIFORM(wall)
/**
 * Check if a given directed line segment intersects the wall.
 *
 * This function is intended to be used to check whether a marker intersects
 * with the wall. If there is an intersection, this function returns the index
 * of the wall element (indexing starts from one). If the marker hits multiple
 * wall elements, only the one that is closest to the start point is returned.
 *
 * @param w_coll Parameter indicating the point of intersection.
 *        This parameter is defined by P = P1 + w_coll * (P2-P1). The value is
 *        one if no intersection occurred.
 * @param r1 Start point R coordinate [m].
 * @param phi1 Start point phi coordinate [rad].
 * @param z1 Start point z coordinate [m].
 * @param r2 End point R coordinate [m].
 * @param phi2 End point phi coordinate [rad].
 * @param z2 End point z coordinate [m].
 * @param wall The wall data.
 *
 * @return Wall element index on intersection, zero otherwise.
 */
size_t Wall_eval_intersection(
    real w_coll[1], real r1, real phi1, real z1, real r2, real phi2, real z2,
    Wall *wall);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(wall)
/**
 * Return the number of wall elements.
 *
 * @param wall The wall data.
 *
 * @return Number of wall elements or zero on failure.
 */
size_t Wall_get_n_elements(Wall *wall);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(wall)
/**
 * Return the flag of a wall element.
 *
 * @param idx The index of the wall element (indexing starts from one).
 * @param wall The wall data.
 *
 * @return Flag of the wall element.
 */
int Wall_get_flag(size_t idx, Wall *wall);
DECLARE_TARGET_END

#endif
