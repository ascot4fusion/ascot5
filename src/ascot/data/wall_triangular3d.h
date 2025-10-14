/**
 * @file wall_triangular3d.h
 * Header file for wall_3d.c
 */
#ifndef WALL_TRIANGULAR3D_H
#define WALL_TRIANGULAR3D_H
#include "defines.h"
#include "parallel.h"
#include "wall.h"

/**
 * Initialize triangular 3D mesh wall model.
 *
 * Constructs the octree which can be time-consuming process.
 *
 * @param wall The struct to initialize.
 * @param n Number of vertices in the wall mesh.
 * @param vertices Cartesian coordinates of the wall mesh vertices.
 *        The layout is following: [..., x1_i, y1_i, z1_i, x2_i, y2_i, z2_i,
 *        x3_i, y3_i, z3_i, ...], where i is the index of the triangle.
 * @param flag Integer label for grouping wall elements together.
 *
 * @return Zero if the initialization succeeded.
 */
int WallTriangular3D_init(
    WallTriangular3D *wall, size_t n, real *vertices, int *flag);

/**
 * Free allocated resources.
 *
 * @param wall The struct whose fields are deallocated.
 */
void WallTriangular3D_free(WallTriangular3D *wall);

/**
 * Offload data to the accelerator.
 *
 * @param wall The struct to offload.
 */
void WallTriangular3D_offload(WallTriangular3D *wall);

GPU_DECLARE_TARGET_SIMD_UNIFORM(wall)
/**
 * Check if a given directed line segment intersects the wall.
 *
 * This function is intended to be used to check whether a marker intersects
 * with the wall. If there is an intersection, this function returns the index
 * of the wall element (indexing starts from one). If the marker hits multiple
 * wall elements, only the one that is closest to the start point is returned.
 *
 * This function first locates the relevant octree cells that are on the marker
 * trajectory, and then checks for intersections with all elements within those
 * cells.
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
size_t WallTriangular3D_eval_intersection(
    real w_coll[1], real r1, real phi1, real z1, real r2, real phi2, real z2,
    WallTriangular3D *wall);
DECLARE_TARGET_END

#endif
