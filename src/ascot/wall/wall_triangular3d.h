/**
 * @file wall_triangular3d.h
 * Header file for wall_3d.c
 */
#ifndef WALL_TRIANGULAR3D_H
#define WALL_TRIANGULAR3D_H
#include "defines.h"
#include "wall.h"
#include "parallel.h"





/**
 * @brief Initialize 3D wall data and check inputs

 * The default octree depth is defined by macro WALL_OCTREE_DEPTH in wall_3d.h.
 *
 * @param data pointer to the data struct
 * @param nelements number of elements in the wall polygon
 * @param vertices each triangle's vertices' (x,y,z) coordinates in an array as
 *        [..., x1_i, y1_i, z1_i, x2_i, y2_i, z2_i, x3_i, y3_i, z3_i, ...]
 * @param flag integer label for grouping wall elements together.
 *
 * @return zero if initialization succeeded
 */
int WallTriangular3D_init(
    WallTriangular3D* wall, int n, real* vertices, int* flag
);


/**
 * @brief Free allocated resources
 *
 * @param data pointer to data struct
 */
void WallTriangular3D_free(WallTriangular3D* wall);


/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void WallTriangular3D_offload(WallTriangular3D* wall);


/**
 * @brief Check if trajectory from (r1, phi1, z1) to (r2, phi2, z2) intersects
 *        the wall using the octree structure
 *
 * @param r1 start point R coordinate [m]
 * @param phi1 start point phi coordinate [rad]
 * @param z1 start point z coordinate [rad]
 * @param r2 end point R coordinate [m]
 * @param phi2 end point phi coordinate [rad]
 * @param z2 end point z coordinate [rad]
 * @param wdata pointer to data struct on target
 * @param w_coll pointer for storing the parameter in P = P1 + w_coll * (P2-P1),
 *        where P is the point where the collision occurred.
 *
 * @return id, which is the first element id if hit, zero otherwise
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(wall)
int WallTriangular3D_hit_wall(
    real r1, real phi1, real z1, real r2, real phi2, real z2, real* w_coll,
    WallTriangular3D* wall
);
DECLARE_TARGET_END





/**
 * @brief Construct wall octree recursively
 *
 * Constructs the octree array by iterating through all wall triangles and
 * placing them into an octree structure
 *
 * Note that this function allocates extra space at the end of the tree_array
 * to store wall element flags.
 *
 * @param w pointer to wall data
 * @param offload_array the offload array
 * @param tree_array pointer to array storing what octree cells contain
 *        which triangles
 */
void WallTriangular3D_init_octree(WallTriangular3D* wall);


/**
 * @brief Construct wall octree iteratively
 *
 * Constructs the octree array by iterating through all wall triangles and
 * octree grid to identify triangles belonging to each grid cell.
 *
 * Slow, only for testing purposes.
 *
 * @param w pointer to wall data
 * @param offload_array offload array
 */
void WallTriangular3D_init_tree(WallTriangular3D* wall, real* offload_array);


/**
 * @brief Check if trajectory from (r1, phi1, z1) to (r2, phi2, z2) intersects
 *        the wall against all triangles
 *
 * @param r1 start point R coordinate [m]
 * @param phi1 start point phi coordinate [rad]
 * @param z1 start point z coordinate [rad]
 * @param r2 end point R coordinate [m]
 * @param phi2 end point phi coordinate [rad]
 * @param z2 end point z coordinate [rad]
 * @param wdata pointer to data struct on target
 * @param w_coll pointer for storing the parameter in P = P1 + w_coll * (P2-P1),
 *        where P is the point where the collision occurred.
 *
 * @return id is wall element id if hit, zero otherwise*
 */
int WallTriangular3D_hit_wall_full(
    real r1, real phi1, real z1, real r2, real phi2, real z2, real* w_coll,
    WallTriangular3D* wall
);

#endif
