/**
 * @file wall_contour2d.h
 * Header file for wall_2d.c
 */
#ifndef WALL_CONTOUR2D_H
#define WALL_CONTOUR2D_H
#include "defines.h"
#include "parallel.h"
#include "wall.h"

/**
 * Initialize 2D contour wall model.
 *
 * @param wall The struct to initialize.
 * @param n Number of vertices in the wall mesh.
 * @param r R coordinates for the wall polygon vertices.
 * @param z z coordinates for the wall polygon vertices.
 * @param flag Integer label for grouping wall elements together.
 *
 * @return Zero if the initialization succeeded.
 */
int WallContour2D_init(
    WallContour2D *wall, size_t n, real r[n], real z[n], int flag[n]);

/**
 * Free allocated resources.
 *
 * @param wall The struct whose fields are deallocated.
 */
void WallContour2D_free(WallContour2D *wall);

/**
 * Offload data to the accelerator.
 *
 * @param wall The struct to offload.
 */
void WallContour2D_offload(WallContour2D *wall);

GPU_DECLARE_TARGET_SIMD_UNIFORM(wall)
/**
 * Check if a given directed line segment intersects the wall.
 *
 * This function is intended to be used to check whether a marker intersects
 * with the wall. If there is an intersection, this function returns the index
 * of the wall element (indexing starts from one). If the marker hits multiple
 * wall elements, only the one that is closest to the start point is returned.
 *
 * Intersections are tested with respect to each wall element.
 *
 * @param w_coll Parameter indicating the point of intersection.
 *        This parameter is defined by P = P1 + w_coll * (P2-P1). The value is
 *        one if no intersection occurred.
 * @param r1 Start point R coordinate [m].
 * @param z1 Start point z coordinate [m].
 * @param r2 End point R coordinate [m].
 * @param z2 End point z coordinate [m].
 * @param wall The wall data.
 *
 * @return Wall element index on intersection, zero otherwise.
 */
size_t WallContour2D_eval_intersection(
    real w_coll[1], real r1, real z1, real r2, real z2, WallContour2D *wall);

#endif
